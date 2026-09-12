#include "NSSurface.hxx"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "cctk_Parameter.h"
#include "util_Table.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace NSSurface {

#ifdef CCTK_MPI
MPI_Comm communicator(const cGH *gh) {
  return CCTK_IsFunctionAliased("GetMPICommWorld")
             ? *static_cast<const MPI_Comm *>(GetMPICommWorld(gh))
             : MPI_COMM_WORLD;
}
#endif

bool all_found(const cGH *gh, int found) {
#ifdef CCTK_MPI
  MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MIN, communicator(gh));
#endif
  return found != 0;
}

// Every rank requests its own points; Carpet routes them to the grid owners.
bool interpolate(const cGH *gh, int table, const Coordinates &coordinates,
                 int nvars, const CCTK_INT *variables,
                 std::vector<CCTK_REAL> *values) {
  const int count = coordinates[0].size();
  std::vector<void *> output(nvars);
  std::vector<CCTK_INT> types(nvars, CCTK_VARIABLE_REAL);
  for (int v = 0; v < nvars; ++v) {
    values[v].assign(count, std::numeric_limits<CCTK_REAL>::quiet_NaN());
    output[v] = values[v].data();
  }
  const void *coords[] = {coordinates[0].data(), coordinates[1].data(),
                         coordinates[2].data()};
  const int error = CCTK_InterpGridArrays(
      gh, 3, CCTK_InterpHandle("Lagrange polynomial interpolation"), table,
      CCTK_CoordSystemHandle("cart3d"), count, CCTK_VARIABLE_REAL, coords,
      nvars, variables, nvars, types.data(), output.data());
  int finite = error >= 0;
  for (int v = 0; v < nvars; ++v)
    for (const auto value : values[v]) finite &= std::isfinite(value);
  return all_found(gh, finite);
}

bool find_radii(const cGH *gh, const Point &origin,
                const std::vector<Point> &directions, CCTK_REAL cutoff,
                CCTK_REAL radius, int samples, CCTK_REAL tolerance,
                std::vector<CCTK_REAL> &radii) {
  const int table = Util_TableCreateFromString("order=1");
  if (table < 0) CCTK_ERROR("NSSurface: cannot create interpolation table");
  const int rays = directions.size();
  Coordinates coordinates;
  for (int axis = 0; axis < 3; ++axis) {
    coordinates[axis].resize(1 + samples * rays);
    coordinates[axis][0] = origin[axis];
  }
  for (int ray = 0; ray < rays; ++ray)
    for (int step = 1; step <= samples; ++step) {
      const CCTK_REAL r = radius * step / samples;
      const int p = 1 + ray * samples + step - 1;
      for (int axis = 0; axis < 3; ++axis)
        coordinates[axis][p] = origin[axis] + r * directions[ray][axis];
    }
  std::vector<CCTK_REAL> density;
  const CCTK_INT input = CCTK_VarIndex("HydroBase::rho");
  bool found = interpolate(gh, table, coordinates, 1, &input, &density);
  std::vector<CCTK_REAL> lower(rays), upper(rays);
  if (found) {
    int bracketed = density[0] > cutoff;
    for (int ray = 0; ray < rays; ++ray) {
      int step = 1;
      while (step <= samples && density[1 + ray * samples + step - 1] > cutoff)
        ++step;
      bracketed &= step <= samples;
      lower[ray] = radius * (step - 1) / samples;
      upper[ray] = radius * step / samples;
    }
    found = all_found(gh, bracketed);
  }
  for (auto &axis : coordinates) axis.resize(rays);
  CCTK_REAL width = radius / samples;
  // Keep all ranks in the same collective sequence, including ranks with no
  // assigned rays. MPI and interpolation calls are outside OpenMP regions.
  for (int iteration = 0; found && width > tolerance && iteration < 64; ++iteration) {
    for (int ray = 0; ray < rays; ++ray) {
      const CCTK_REAL r = (lower[ray] + upper[ray]) / 2;
      for (int axis = 0; axis < 3; ++axis)
        coordinates[axis][ray] = origin[axis] + r * directions[ray][axis];
    }
    found = interpolate(gh, table, coordinates, 1, &input, &density);
    if (!found) break;
    for (int ray = 0; ray < rays; ++ray) {
      const CCTK_REAL midpoint = (lower[ray] + upper[ray]) / 2;
      if (density[ray] > cutoff) lower[ray] = midpoint;
      else upper[ray] = midpoint;
    }
    width /= 2;
  }
  found = found && width <= tolerance;
  if (found) {
    radii.resize(rays);
    for (int ray = 0; ray < rays; ++ray)
      radii[ray] = (lower[ray] + upper[ray]) / 2;
  }
  Util_TableDestroy(table);
  return found;
}
} // namespace NSSurface

extern "C" void NSSurface_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const auto *spherical = static_cast<const CCTK_INT *>(
      CCTK_ParameterGet("set_spherical", "SphericalSurface", nullptr));
  const auto *elliptic = static_cast<const CCTK_INT *>(
      CCTK_ParameterGet("set_elliptic", "SphericalSurface", nullptr));
  const auto *ox = static_cast<const CCTK_REAL *>(
      CCTK_ParameterGet("origin_x", "SphericalSurface", nullptr));
  const auto *oy = static_cast<const CCTK_REAL *>(
      CCTK_ParameterGet("origin_y", "SphericalSurface", nullptr));
  const auto *oz = static_cast<const CCTK_REAL *>(
      CCTK_ParameterGet("origin_z", "SphericalSurface", nullptr));

  std::fill_n(info_valid, nsurfaces, 0);
  for (int s = 0; s < 42; ++s) {
    if (!find_surface[s]) continue;
    if (!CCTK_Equals(metric_type, "physical"))
      CCTK_ERROR("NSSurface requires the physical ADMBase metric for surface information");
    if (s >= nsurfaces || spherical[s] || elliptic[s] || auto_res[s])
      CCTK_VERROR("NSSurface[%d]: use an allocated surface with "
                  "set_spherical, set_elliptic and auto_res disabled", s);
    const int nt = sf_ntheta[s] - 2*sf_nghoststheta[s];
    const int np = sf_nphi[s] - 2*sf_nghostsphi[s];
    if (sf_nghoststheta[s] < 2 || sf_nghostsphi[s] < 2)
      CCTK_VERROR("NSSurface[%d]: fourth-order surface derivatives require at least two angular ghost points in theta and phi", s);
    if ((symmetric_z[s] ? 2*nt-1 : nt) < 3 ||
        np < (symmetric_x[s] || symmetric_y[s] ? 2 : 4))
      CCTK_VERROR("NSSurface[%d]: angular grid is too small for surface integration", s);
    if (!std::isfinite(max_radius[s]) || !std::isfinite(rho_surface[s]) ||
        !std::isfinite(radial_tolerance) || radial_tolerance >= max_radius[s])
      CCTK_VERROR("NSSurface[%d]: invalid radius, density or radial tolerance", s);
    if (CCTK_InterpHandle("Lagrange polynomial interpolation") < 0 ||
        CCTK_CoordSystemHandle("cart3d") < 0)
      CCTK_ERROR("NSSurface needs Cartesian grid interpolation and AEILocalInterp");
    sf_origin_x[s] = ox[s];
    sf_origin_y[s] = oy[s];
    sf_origin_z[s] = oz[s];
    sf_active[s] = 1;
    sf_valid[s] = 0;
    std::fill_n(sf_radius + s * maxntheta * maxnphi, maxntheta * maxnphi, 0.0);
  }
}

extern "C" void NSSurface_Find(CCTK_ARGUMENTS) {
  using namespace NSSurface;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  for (int s = 0; s < nsurfaces; ++s) {
    if (!find_surface[s]) continue;
    sf_active[s] = 1;
    sf_valid[s] = 0; // Never let QLM interpret an old contour as current data.
    if (cctk_iteration % compute_every != 0) continue;
    const Point origin{sf_origin_x[s], sf_origin_y[s], sf_origin_z[s]};
    const int nt = sf_ntheta[s], np = sf_nphi[s];
    const int count = nt * np;
    const int rank = CCTK_MyProc(cctkGH), ranks = CCTK_nProcs(cctkGH);
    const int first = std::size_t(count) * rank / ranks;
    const int last = std::size_t(count) * (rank + 1) / ranks;
    std::vector<Point> directions(last - first);
    for (int ray = first; ray < last; ++ray) {
      const CCTK_REAL theta = sf_origin_theta[s] + (ray % nt) * sf_delta_theta[s];
      const CCTK_REAL phi = sf_origin_phi[s] + (ray / nt) * sf_delta_phi[s];
      directions[ray - first] = {std::sin(theta) * std::cos(phi),
                                 std::sin(theta) * std::sin(phi), std::cos(theta)};
    }
    std::vector<CCTK_REAL> local_radii;
    if (!all_found(cctkGH, std::all_of(origin.begin(), origin.end(),
                                      [](CCTK_REAL x) { return std::isfinite(x); })) ||
        !find_radii(cctkGH, origin, directions, rho_surface[s], max_radius[s],
                    radial_samples, radial_tolerance, local_radii)) {
      sf_valid[s] = -1;
      if (CCTK_MyProc(cctkGH) == 0)
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "NSSurface[%d]: no complete rho=%g contour about (%g,%g,%g) "
                 "within r=%g; check origin, search bounds and interpolation", s,
                 double(rho_surface[s]), double(origin[0]), double(origin[1]),
                 double(origin[2]), double(max_radius[s]));
      continue;
    }
    std::vector<CCTK_REAL> radii(count);
#ifdef CCTK_MPI
    std::vector<int> sizes(ranks), offsets(ranks);
    for (int p = 0; p < ranks; ++p) {
      const int begin = std::size_t(count) * p / ranks;
      const int end = std::size_t(count) * (p + 1) / ranks;
      sizes[p] = (end - begin) * sizeof(CCTK_REAL);
      offsets[p] = begin * sizeof(CCTK_REAL);
    }
    MPI_Allgatherv(local_radii.data(), local_radii.size() * sizeof(CCTK_REAL), MPI_BYTE,
                   radii.data(), sizes.data(), offsets.data(), MPI_BYTE, communicator(cctkGH));
#else
    radii.swap(local_radii);
#endif
    for (int j = 0; j < np; ++j)
      for (int i = 0; i < nt; ++i)
        sf_radius[i + maxntheta * (j + maxnphi * s)] = radii[i + nt * j];
    sf_valid[s] = 1;
  }
}
