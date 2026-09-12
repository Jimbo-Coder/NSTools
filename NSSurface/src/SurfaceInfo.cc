#include "NSSurface.hxx"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace {
using NSSurface::Point;

// Centred fourth-order radius derivative; initialization requires two ghosts.
CCTK_REAL derivative(const CCTK_REAL *f, int stride, CCTK_REAL dx) {
  return (f[-2*stride] - 8*f[-stride] + 8*f[stride] - f[2*stride]) / (12*dx);
}

void reduce(const cGH *gh, const char *operation, int count, CCTK_REAL *values) {
  std::vector<CCTK_REAL> result(count);
  const int handle = CCTK_ReductionArrayHandle(operation);
  if (handle < 0 || CCTK_ReduceLocArrayToArray1D(
          gh, -1, handle, values, result.data(), count, CCTK_VARIABLE_REAL) < 0)
    CCTK_ERROR("NSSurface: cannot reduce surface information");
  std::copy(result.begin(), result.end(), values);
}
} // namespace

extern "C" void NSSurface_UpdateInfo(CCTK_ARGUMENTS) {
  using namespace NSSurface;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL nan = std::numeric_limits<CCTK_REAL>::quiet_NaN();
  const CCTK_REAL infinity = std::numeric_limits<CCTK_REAL>::infinity();
  const CCTK_INT variables[] = {
      CCTK_VarIndex("ADMBase::gxx"), CCTK_VarIndex("ADMBase::gxy"),
      CCTK_VarIndex("ADMBase::gxz"), CCTK_VarIndex("ADMBase::gyy"),
      CCTK_VarIndex("ADMBase::gyz"), CCTK_VarIndex("ADMBase::gzz")};

  for (int s = 0; s < nsurfaces; ++s) {
    info_valid[s] = 0;
    if (!find_surface[s]) continue;
    CCTK_REAL *const info[] = {
        sf_area, sf_mean_radius, sf_centroid_x, sf_centroid_y, sf_centroid_z,
        sf_quadrupole_xx, sf_quadrupole_xy, sf_quadrupole_xz,
        sf_quadrupole_yy, sf_quadrupole_yz, sf_quadrupole_zz,
        sf_min_radius, sf_max_radius, sf_min_x, sf_min_y, sf_min_z,
        sf_max_x, sf_max_y, sf_max_z};
    for (auto field : info) field[s] = nan;
    if (sf_valid[s] <= 0) {
      info_valid[s] = sf_valid[s];
      continue;
    }

    const int nt = sf_ntheta[s], np = sf_nphi[s];
    const int gt = sf_nghoststheta[s], gp = sf_nghostsphi[s];
    // SphericalSurface allows steering these widths after initialization.
    if (gt < 2 || gp < 2) {
      info_valid[s] = -1;
      if (CCTK_MyProc(cctkGH) == 0)
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "NSSurface[%d]: surface information requires two angular ghost points in theta and phi", s);
      continue;
    }
    const int ni = nt - 2*gt, nj = np - 2*gp, count = ni*nj;
    const int rank = CCTK_MyProc(cctkGH), ranks = CCTK_nProcs(cctkGH);
    const int first = std::size_t(count)*rank/ranks;
    const int last = std::size_t(count)*(rank+1)/ranks;
    const Point origin{sf_origin_x[s], sf_origin_y[s], sf_origin_z[s]};
    const bool reflected[] = {bool(symmetric_x[s]), bool(symmetric_y[s]), bool(symmetric_z[s])};
    const CCTK_REAL dt = sf_delta_theta[s], dp = sf_delta_phi[s];
    const auto *radius = sf_radius + maxntheta*maxnphi*s;
    Coordinates coordinates;
    for (auto &axis : coordinates) axis.resize(last-first);
    for (int ray = first; ray < last; ++ray) {
      const int i = gt + ray%ni, j = gp + ray/ni;
      const CCTK_REAL theta = sf_origin_theta[s] + i*dt;
      const CCTK_REAL phi = sf_origin_phi[s] + j*dp;
      const CCTK_REAL r = radius[i + maxntheta*j];
      coordinates[0][ray-first] = origin[0] + r*std::sin(theta)*std::cos(phi);
      coordinates[1][ray-first] = origin[1] + r*std::sin(theta)*std::sin(phi);
      coordinates[2][ray-first] = origin[2] + r*std::cos(theta);
    }
    const int table = Util_TableCreateFromString("order=1");
    if (table < 0) CCTK_ERROR("NSSurface: cannot create metric interpolation table");
    std::vector<CCTK_REAL> metric[6];
    bool valid = interpolate(cctkGH, table, coordinates, 6, variables, metric);
    Util_TableDestroy(table);

    // The full midpoint theta grid has 2*ni-1 points for a reflected hemisphere.
    // Use QLM's sine-series quadrature; halve endpoints on symmetry planes.
    const int full_nt = symmetric_z[s] ? 2*ni-1 : ni;
    std::vector<CCTK_REAL> weights(ni);
    for (int i = 0; i < ni; ++i) {
      const CCTK_REAL theta = sf_origin_theta[s] + (i+gt)*dt;
      CCTK_REAL sum = 0;
      for (int l = 0; l <= (full_nt-1)/2; ++l)
        sum += std::sin((2*l+1)*theta)/(2*l+1);
      weights[i] = 4*dp*sum/full_nt;
      for (const bool reflection : reflected) if (reflection) weights[i] *= 2;
      if (symmetric_z[s] && i == ni-1) weights[i] /= 2;
    }

    CCTK_REAL sums[11] = {}; // area, radius, x/y/z, xx/xy/xz/yy/yz/zz moments
    CCTK_REAL lower[] = {infinity, infinity, infinity, infinity}; // radius, x/y/z
    CCTK_REAL upper[] = {-infinity, -infinity, -infinity, -infinity};
    for (int ray = first; valid && ray < last; ++ray) {
      const int p = ray-first, i = gt + ray%ni, j = gp + ray/ni;
      const CCTK_REAL theta = sf_origin_theta[s] + i*dt;
      const CCTK_REAL phi = sf_origin_phi[s] + j*dp;
      const CCTK_REAL st = std::sin(theta), ct = std::cos(theta);
      const CCTK_REAL sp = std::sin(phi), cp = std::cos(phi);
      const Point direction{st*cp, st*sp, ct};
      const CCTK_REAL r = radius[i + maxntheta*j];
      const auto *point = radius + i + maxntheta*j;
      const CCTK_REAL rt = derivative(point, 1, dt);
      const CCTK_REAL rp = derivative(point, maxntheta, dp);
      const Point tangent_t{rt*direction[0]+r*ct*cp, rt*direction[1]+r*ct*sp, rt*ct-r*st};
      const Point tangent_p{rp*direction[0]-r*st*sp, rp*direction[1]+r*st*cp, rp*ct};
      const CCTK_REAL gxx = metric[0][p], gxy = metric[1][p], gxz = metric[2][p];
      const CCTK_REAL gyy = metric[3][p], gyz = metric[4][p], gzz = metric[5][p];
      const auto dot = [&](const Point &a, const Point &b) {
        return gxx*a[0]*b[0] + gyy*a[1]*b[1] + gzz*a[2]*b[2]
            + gxy*(a[0]*b[1]+a[1]*b[0]) + gxz*(a[0]*b[2]+a[2]*b[0])
            + gyz*(a[1]*b[2]+a[2]*b[1]);
      };
      const CCTK_REAL determinant = gxx*(gyy*gzz-gyz*gyz)
          - gxy*(gxy*gzz-gxz*gyz) + gxz*(gxy*gyz-gxz*gyy);
      // Pull back the physical 3-metric to obtain dA = sqrt(det(q)) dtheta dphi.
      const CCTK_REAL tt = dot(tangent_t, tangent_t), tp = dot(tangent_t, tangent_p);
      const CCTK_REAL pp = dot(tangent_p, tangent_p), area2 = tt*pp-tp*tp;
      valid = gxx > 0 && gxx*gyy > gxy*gxy && determinant > 0
          && std::isfinite(area2) && area2 > 0;
      if (!valid) break;
      CCTK_REAL weight = std::sqrt(area2)*weights[i-gt];
      if ((symmetric_x[s] || symmetric_y[s]) && (j == gp || j == np-gp-1)) weight /= 2;
      const Point position{r*direction[0], r*direction[1], r*direction[2]};
      sums[0] += weight;
      sums[1] += weight*r;
      for (int a = 0; a < 3; ++a) sums[2+a] += weight*position[a];
      int moment = 5;
      for (int a = 0; a < 3; ++a)
        for (int b = a; b < 3; ++b) sums[moment++] += weight*position[a]*position[b];
      lower[0] = std::min(lower[0], r);
      upper[0] = std::max(upper[0], r);
      for (int a = 0; a < 3; ++a) {
        lower[1+a] = std::min(lower[1+a], reflected[a] ? -std::abs(position[a]) : position[a]);
        upper[1+a] = std::max(upper[1+a], reflected[a] ? std::abs(position[a]) : position[a]);
      }
    }
    for (const auto sum : sums) valid = valid && std::isfinite(sum);
    if (!all_found(cctkGH, valid)) {
      info_valid[s] = -1;
      if (rank == 0)
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "NSSurface[%d]: cannot calculate surface information; check metric and angular resolution", s);
      continue;
    }
    reduce(cctkGH, "sum", 11, sums);
    reduce(cctkGH, "minimum", 4, lower);
    reduce(cctkGH, "maximum", 4, upper);
    sf_area[s] = sums[0];
    sf_mean_radius[s] = sums[1]/sums[0];
    Point centroid;
    for (int a = 0; a < 3; ++a) centroid[a] = reflected[a] ? 0 : sums[2+a]/sums[0];
    sf_centroid_x[s] = origin[0]+centroid[0];
    sf_centroid_y[s] = origin[1]+centroid[1];
    sf_centroid_z[s] = origin[2]+centroid[2];
    int moment = 5;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b) {
        const CCTK_REAL average = a != b && (reflected[a] || reflected[b]) ? 0 : sums[moment]/sums[0];
        sums[moment] = average-centroid[a]*centroid[b];
        ++moment;
      }
    sf_quadrupole_xx[s] = sums[5];
    sf_quadrupole_xy[s] = sums[6];
    sf_quadrupole_xz[s] = sums[7];
    sf_quadrupole_yy[s] = sums[8];
    sf_quadrupole_yz[s] = sums[9];
    sf_quadrupole_zz[s] = sums[10];
    sf_min_radius[s] = lower[0];
    sf_max_radius[s] = upper[0];
    sf_min_x[s] = origin[0]+lower[1]; sf_max_x[s] = origin[0]+upper[1];
    sf_min_y[s] = origin[1]+lower[2]; sf_max_y[s] = origin[1]+upper[2];
    sf_min_z[s] = origin[2]+lower[3]; sf_max_z[s] = origin[2]+upper[3];
    info_valid[s] = 1;
  }
}
