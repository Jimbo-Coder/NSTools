#include "NSSurface.hxx"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameter.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace NSSurface {
namespace {
constexpr int max_slots = 42;
constexpr CCTK_REAL unset = std::numeric_limits<CCTK_REAL>::quiet_NaN();
constexpr CCTK_REAL lowest = -std::numeric_limits<CCTK_REAL>::max();

// Slot states, kept in NSSurface::slot_state so they survive restarts.
enum SlotState { waiting = 0, following = 1 };

using Candidate = NSTools::Extremum; // value: the local maximum's density

// Per-step scratch filled by NSSurface_ScanDensity: each followed star's peak
// density, the lowest density on the grid (the atmosphere), the finest grid
// spacing, the grid's bounding box and, while a slot waits to be placed, this
// rank's densest local density maxima.
std::array<CCTK_REAL, max_slots> slot_peaks;
CCTK_REAL grid_floor;
CCTK_REAL finest_spacing;
Point grid_min, grid_max;
NSTools::Extrema local_maxima;
std::array<bool, max_slots> warned_floor{}; // atmosphere warning given

// Rays sample a fixed radial lattice: max_radius/64 when max_radius is set,
// else the finest grid spacing (recorded per slot at its first search, so
// restarts and later regridding do not change it). Roots are bisected to
// 2^-12 of that step.
constexpr int lattice_intervals = 64;
constexpr CCTK_REAL tolerance_fraction = 1.0 / 4096;
// Automatic max_radius: this many times the star's largest unclipped radius
// at the previous search (no limit before the first surface).
constexpr CCTK_REAL reach_factor = 2.0;

struct SearchResult {
  bool ok = false;
  Point origin{};
  std::vector<CCTK_REAL> radii; // every ray, replicated on all ranks
  Point shift{};                // contour volume centroid minus origin
  CCTK_REAL cutoff = unset, peak = unset, limit = unset, clipped = unset;
  CCTK_REAL reach = 0; // largest unclipped radius
};

// SphericalSurface's private per-surface parameters.
template <typename T> const T *surface_parameter(const char *name) {
  return static_cast<const T *>(CCTK_ParameterGet(name, "SphericalSurface", nullptr));
}

bool given_origin(const Point &p) {
  return p[0] != 0.0 || p[1] != 0.0 || p[2] != 0.0;
}

std::array<bool, 3> grid_reflections() {
  if (!CCTK_IsThornActive("ReflectionSymmetry"))
    return {false, false, false};
  return {parameter_int("reflection_x", "ReflectionSymmetry", 0) != 0,
          parameter_int("reflection_y", "ReflectionSymmetry", 0) != 0,
          parameter_int("reflection_z", "ReflectionSymmetry", 0) != 0};
}

// Axes whose negative half is covered by a symmetry rather than by the grid:
// reflections, and the half (x) or quadrant (x, y) of rotating symmetries.
std::array<bool, 3> mirrored_axes(const std::array<bool, 3> &reflections) {
  std::array<bool, 3> mirrored = reflections;
  mirrored[0] |= CCTK_IsThornActive("RotatingSymmetry180") ||
                 CCTK_IsThornActive("RotatingSymmetry90");
  mirrored[1] |= CCTK_IsThornActive("RotatingSymmetry90");
  return mirrored;
}

// Distance from origin along dir to the edge of the grid's bounding box,
// including the parts covered by symmetry.
CCTK_REAL distance_to_grid_edge(const Point &origin, const Point &dir) {
  CCTK_REAL r = std::numeric_limits<CCTK_REAL>::infinity();
  for (int a = 0; a < 3; ++a) {
    if (dir[a] > 0)
      r = std::min(r, (grid_max[a] - origin[a]) / dir[a]);
    else if (dir[a] < 0)
      r = std::min(r, (grid_min[a] - origin[a]) / dir[a]);
  }
  return std::max(r, CCTK_REAL(0));
}

// Search limit about slot s's origin.
CCTK_REAL radius_limit(int s, CCTK_REAL outer_radius) {
  DECLARE_CCTK_PARAMETERS;
  if (max_radius[s] > 0)
    return max_radius[s];
  return outer_radius > 0 ? reach_factor * outer_radius
                          : std::numeric_limits<CCTK_REAL>::infinity();
}

CCTK_REAL contour_density(int s, CCTK_REAL peak, CCTK_REAL reference) {
  DECLARE_CCTK_PARAMETERS;
  if (CCTK_Equals(rho_surface_mode, "absolute"))
    return rho_surface[s];
  if (CCTK_Equals(rho_surface_mode, "initial") && reference > 0)
    return rho_surface_fraction[s] * reference;
  return rho_surface_fraction[s] * peak;
}

// Volume centroid of the star-shaped region inside radii(theta, phi),
// relative to the origin the radii are measured from.
Point volume_centroid(int nt, int np, int gt, int gp, CCTK_REAL origin_theta,
                      CCTK_REAL origin_phi, CCTK_REAL dt, CCTK_REAL dp,
                      const std::array<bool, 3> &reflected,
                      const std::vector<CCTK_REAL> &radii) {
  const std::vector<CCTK_REAL> weights =
      theta_weights(nt - 2 * gt, gt, origin_theta, dt, dp, reflected);
  CCTK_REAL volume = 0.0;
  Point moment{0.0, 0.0, 0.0};
  for (int j = gp; j < np - gp; ++j) {
    const CCTK_REAL phi = origin_phi + j * dp;
    const CCTK_REAL half =
        (reflected[0] || reflected[1]) && (j == gp || j == np - gp - 1) ? 0.5 : 1.0;
    for (int i = gt; i < nt - gt; ++i) {
      const CCTK_REAL theta = origin_theta + i * dt;
      const Point dir{std::sin(theta) * std::cos(phi),
                      std::sin(theta) * std::sin(phi), std::cos(theta)};
      // dV = r^2 sin(theta) dr dtheta dphi
      const CCTK_REAL w = half * weights[i - gt] * std::sin(theta);
      const CCTK_REAL r3 = std::pow(radii[i + nt * j], 3);
      volume += w * r3 / 3;
      for (int a = 0; a < 3; ++a)
        moment[a] += w * r3 * radii[i + nt * j] / 4 * dir[a];
    }
  }
  Point shift{0.0, 0.0, 0.0};
  if (volume > 0 && std::isfinite(volume))
    for (int a = 0; a < 3; ++a)
      shift[a] = reflected[a] ? 0.0 : moment[a] / volume;
  return shift;
}

std::vector<CCTK_REAL> interpolate_rho(const cGH *gh, const Coordinates &points,
                                       bool &ok) {
  const int table = Util_TableCreateFromString("order=1");
  if (table < 0)
    CCTK_ERROR("NSSurface: cannot create interpolation table");
  const CCTK_INT rho = CCTK_VarIndex("HydroBase::rho");
  std::vector<CCTK_REAL> density;
  ok = interpolate(gh, table, points, 1, &rho, &density);
  Util_TableDestroy(table);
  return density;
}

// First outward crossing of rho = cutoff on each of this rank's rays.
// Rays are sampled on the fixed lattice r_k = k*step, each up to its own
// limit (the last sample sits on the limit). The first pass stops at `start`
// (about the previous surface); later passes double it for rays still inside
// the star. The brackets, and hence the radii, do not depend on `start`, which
// only saves samples. Rays with no crossing inside their limit are clipped at
// it and flagged in `clipped`. Returns false only if the origin is outside the
// contour or interpolation fails. Every rank makes the same collective calls.
bool find_radii(const cGH *gh, const Point &origin,
                const std::vector<Point> &directions,
                const std::vector<CCTK_REAL> &limits, CCTK_REAL cutoff,
                CCTK_REAL start, CCTK_REAL step, CCTK_REAL tolerance,
                std::vector<CCTK_REAL> &radii, std::vector<char> &clipped) {
  const int rays = directions.size();
  std::vector<CCTK_REAL> lower(rays, 0.0), upper(rays, 0.0);
  std::vector<char> open(rays, 1); // no crossing bracketed yet
  std::vector<long> samples(rays); // lattice points up to each ray's limit
  for (int ray = 0; ray < rays; ++ray)
    samples[ray] = std::max(1L, long(std::ceil(limits[ray] / step - 1e-12)));
  auto radius = [&](int ray, long k) { return std::min(k * step, limits[ray]); };
  auto at = [&](int ray, CCTK_REAL r, Coordinates &points, std::size_t p) {
    for (int a = 0; a < 3; ++a)
      points[a][p] = origin[a] + r * directions[ray][a];
  };

  bool ok;
  std::vector<CCTK_REAL> density =
      interpolate_rho(gh, {{{origin[0]}, {origin[1]}, {origin[2]}}}, ok);
  ok = ok && all_found(gh, density[0] > cutoff);

  for (long inner = 0, outer = std::max(1L, long(std::ceil(start / step))); ok;
       inner = outer, outer *= 2) {
    std::vector<std::pair<int, long>> samples_now; // (ray, k)
    for (int ray = 0; ray < rays; ++ray)
      for (long k = inner + 1; open[ray] && k <= std::min(outer, samples[ray]); ++k)
        samples_now.emplace_back(ray, k);
    Coordinates points;
    for (auto &axis : points)
      axis.resize(samples_now.size());
    for (std::size_t p = 0; p < samples_now.size(); ++p)
      at(samples_now[p].first, radius(samples_now[p].first, samples_now[p].second),
         points, p);
    density = interpolate_rho(gh, points, ok);
    if (!ok)
      break;
    for (std::size_t p = 0; p < samples_now.size(); ++p) {
      const auto [ray, k] = samples_now[p];
      if (open[ray] && density[p] <= cutoff) {
        open[ray] = 0;
        lower[ray] = radius(ray, k - 1);
        upper[ray] = radius(ray, k);
      }
    }
    bool more = false; // open rays with lattice left beyond this pass
    for (int ray = 0; ray < rays; ++ray)
      more |= open[ray] && samples[ray] > outer;
    if (all_found(gh, !more))
      break;
  }
  if (!ok)
    return false;

  clipped = open;
  CCTK_REAL width = 0.0;
  for (int ray = 0; ray < rays; ++ray) {
    if (open[ray])
      lower[ray] = upper[ray] = limits[ray];
    width = std::max(width, upper[ray] - lower[ray]);
  }
  reduce(gh, "maximum", 1, &width);
  // Bisect every bracket together; clipped rays stay at their limit.
  Coordinates points;
  for (auto &axis : points)
    axis.resize(rays);
  for (int iteration = 0; width > tolerance && iteration < 64; ++iteration) {
    for (int ray = 0; ray < rays; ++ray)
      at(ray, (lower[ray] + upper[ray]) / 2, points, ray);
    density = interpolate_rho(gh, points, ok);
    if (!ok)
      return false;
    for (int ray = 0; ray < rays; ++ray)
      (density[ray] > cutoff ? lower : upper)[ray] = (lower[ray] + upper[ray]) / 2;
    width /= 2;
  }
  if (width > tolerance)
    return false;
  radii.resize(rays);
  for (int ray = 0; ray < rays; ++ray)
    radii[ray] = (lower[ray] + upper[ray]) / 2;
  return true;
}

// True when rho stays above cutoff along the straight segment from a to b:
// both points lie in one connected body at this contour level.
bool connected(const cGH *gh, const Point &a, const Point &b, CCTK_REAL cutoff) {
  CCTK_REAL lowest_density;
  return segment_minimum(gh, CCTK_VarIndex("HydroBase::rho"), a, b, segment_samples,
                                  [](CCTK_REAL v) { return v; }, lowest_density) &&
         lowest_density > cutoff;
}

void gather_radii(const cGH *gh, int count, const std::vector<CCTK_REAL> &local,
                  std::vector<CCTK_REAL> &radii) {
  radii.resize(count);
#ifdef CCTK_MPI
  const int ranks = CCTK_nProcs(gh);
  std::vector<int> sizes(ranks), offsets(ranks);
  for (int p = 0; p < ranks; ++p) {
    const int begin = std::size_t(count) * p / ranks;
    const int end = std::size_t(count) * (p + 1) / ranks;
    sizes[p] = (end - begin) * sizeof(CCTK_REAL);
    offsets[p] = begin * sizeof(CCTK_REAL);
  }
  MPI_Allgatherv(local.data(), local.size() * sizeof(CCTK_REAL), MPI_BYTE,
                 radii.data(), sizes.data(), offsets.data(), MPI_BYTE,
                 communicator(gh));
#else
  radii = local;
#endif
}
} // namespace

bool surface_owned(int s) {
  static const auto *find = static_cast<const CCTK_INT *>(
      CCTK_ParameterGet("find_surface", "NSSurface", nullptr));
  return s >= 0 && s < max_slots && find != nullptr && find[s];
}

// QLM's sine-series quadrature in theta for integrals of f(theta) dtheta dphi
// over the ni interior midpoints; a reflected hemisphere uses the full
// 2*ni-1 point rule with its equator halved, and each reflection doubles.
std::vector<CCTK_REAL> theta_weights(int ni, int gt, CCTK_REAL origin_theta,
                                     CCTK_REAL dt, CCTK_REAL dp,
                                     const std::array<bool, 3> &reflected) {
  const int full_nt = reflected[2] ? 2 * ni - 1 : ni;
  std::vector<CCTK_REAL> weights(ni);
  for (int i = 0; i < ni; ++i) {
    const CCTK_REAL theta = origin_theta + (i + gt) * dt;
    CCTK_REAL sum = 0;
    for (int l = 0; l <= (full_nt - 1) / 2; ++l)
      sum += std::sin((2 * l + 1) * theta) / (2 * l + 1);
    weights[i] = 4 * dp * sum / full_nt;
    for (const bool reflection : reflected)
      if (reflection)
        weights[i] *= 2;
    if (reflected[2] && i == ni - 1)
      weights[i] /= 2;
  }
  return weights;
}

} // namespace NSSurface

extern "C" void NSSurface_Init(CCTK_ARGUMENTS) {
  using namespace NSSurface;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  const bool absolute = CCTK_Equals(rho_surface_mode, "absolute");
  const auto *set_spherical = surface_parameter<CCTK_INT>("set_spherical");
  const auto *set_elliptic = surface_parameter<CCTK_INT>("set_elliptic");
  const auto *origin_x = surface_parameter<CCTK_REAL>("origin_x");
  const auto *origin_y = surface_parameter<CCTK_REAL>("origin_y");
  const auto *origin_z = surface_parameter<CCTK_REAL>("origin_z");

  for (int s = 0; s < nsurfaces; ++s) {
    info_valid[s] = 0;
    rho_contour[s] = rho_peak[s] = search_radius[s] = clipped_fraction[s] = unset;
    slot_state[s] = waiting;
    rho_reference[s] = outer_radius[s] = radial_step[s] = 0.0;
    last_origin_x[s] = last_origin_y[s] = last_origin_z[s] = 0.0;
  }
  for (int s = 0; s < max_slots; ++s) {
    if (!surface_owned(s))
      continue;
    if (!CCTK_Equals(metric_type, "physical"))
      CCTK_ERROR("NSSurface requires the physical ADMBase metric for surface information");
    if (s >= nsurfaces || set_spherical[s] || set_elliptic[s] || auto_res[s])
      CCTK_VERROR("NSSurface[%d]: use an allocated surface with "
                  "set_spherical, set_elliptic and auto_res disabled", s);
    const int nt = sf_ntheta[s] - 2 * sf_nghoststheta[s];
    const int np = sf_nphi[s] - 2 * sf_nghostsphi[s];
    if (sf_nghoststheta[s] < 2 || sf_nghostsphi[s] < 2)
      CCTK_VERROR("NSSurface[%d]: fourth-order surface derivatives require at "
                  "least two angular ghost points in theta and phi", s);
    if ((symmetric_z[s] ? 2 * nt - 1 : nt) < 3 ||
        np < (symmetric_x[s] || symmetric_y[s] ? 2 : 4))
      CCTK_VERROR("NSSurface[%d]: angular grid is too small for surface integration", s);
    if (!std::isfinite(max_radius[s]) ||
        !std::isfinite(absolute ? rho_surface[s] : rho_surface_fraction[s]))
      CCTK_VERROR("NSSurface[%d]: invalid radius or density", s);
    if (CCTK_InterpHandle("Lagrange polynomial interpolation") < 0 ||
        CCTK_CoordSystemHandle("cart3d") < 0)
      CCTK_ERROR("NSSurface needs Cartesian grid interpolation and AEILocalInterp");
    sf_origin_x[s] = last_origin_x[s] = origin_x[s];
    sf_origin_y[s] = last_origin_y[s] = origin_y[s];
    sf_origin_z[s] = last_origin_z[s] = origin_z[s];
    // A slot left at the coordinate origin is placed on a star by the search.
    slot_state[s] = given_origin({origin_x[s], origin_y[s], origin_z[s]}) ? following : waiting;
    sf_active[s] = 1;
    sf_valid[s] = 0;
    std::fill_n(sf_radius + s * maxntheta * maxnphi, maxntheta * maxnphi, 0.0);
  }
}

// Checkpoints written before NSSurface kept this state recover it unset.
extern "C" void NSSurface_Recover(CCTK_ARGUMENTS) {
  using namespace NSSurface;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  for (int s = 0; s < nsurfaces && s < max_slots; ++s) {
    if (!surface_owned(s))
      continue;
    const Point origin{sf_origin_x[s], sf_origin_y[s], sf_origin_z[s]};
    if (!(rho_reference[s] >= 0 && std::isfinite(rho_reference[s])))
      rho_reference[s] = 0.0;
    // Unset state (rho_reference also unset) follows the stored origin.
    if ((slot_state[s] != waiting && slot_state[s] != following) ||
        (rho_reference[s] == 0 && given_origin(origin)))
      slot_state[s] = given_origin(origin) ? following : waiting;
    if (!(outer_radius[s] >= 0 && std::isfinite(outer_radius[s])))
      outer_radius[s] = 0.0;
    if (!(radial_step[s] >= 0 && std::isfinite(radial_step[s])))
      radial_step[s] = 0.0; // set again at the next search
    if (!finite({last_origin_x[s], last_origin_y[s], last_origin_z[s]})) {
      last_origin_x[s] = origin[0];
      last_origin_y[s] = origin[1];
      last_origin_z[s] = origin[2];
    }
  }
}

extern "C" void NSSurface_BeginScan(CCTK_ARGUMENTS) {
  using namespace NSSurface;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  const int every = output_cadence(compute_every);
  if (every <= 0 || cctk_iteration % every != 0)
    return;
  slot_peaks.fill(lowest);
  grid_floor = std::numeric_limits<CCTK_REAL>::max();
  finest_spacing = std::numeric_limits<CCTK_REAL>::max();
  local_maxima.clear();
  grid_min.fill(std::numeric_limits<CCTK_REAL>::max());
  grid_max.fill(lowest);
}

// Each point counts toward the bounding box and the peak of the nearest
// followed star within that slot's search limit. While a slot waits to be placed,
// local density maxima (at least all 26 neighbours, above at least one) are
// collected as candidates; flat atmosphere is not a maximum.
extern "C" void NSSurface_ScanDensity(CCTK_ARGUMENTS) {
  using namespace NSSurface;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  const int every = output_cadence(compute_every);
  if (every <= 0 || cctk_iteration % every != 0)
    return;

  bool any_owned = false, any_waiting = false;
  std::vector<int> followed;
  for (int s = 0; s < nsurfaces && s < max_slots; ++s) {
    any_owned |= surface_owned(s);
    any_waiting |= surface_owned(s) && slot_state[s] == waiting;
    if (surface_owned(s) && slot_state[s] == following &&
        finite({sf_origin_x[s], sf_origin_y[s], sf_origin_z[s]}))
      followed.push_back(s);
  }
  if (!any_owned)
    return;

  static const int weight_index = CCTK_VarIndex("CarpetReduce::weight");
  const auto *weight = weight_index >= 0
                           ? static_cast<const CCTK_REAL *>(
                                 CCTK_VarDataPtrI(cctkGH, 0, weight_index))
                           : nullptr;
  const int dj = cctk_ash[0], dk = cctk_ash[0] * cctk_ash[1];
  auto density = [&](int q) { return rho[q]; };
  for (int d = 0; d < 3; ++d)
    finest_spacing = std::min(finest_spacing, cctk_delta_space[d] / cctk_levfac[d]);
  for (int k = cctk_nghostzones[2]; k < cctk_lsh[2] - cctk_nghostzones[2]; ++k)
    for (int j = cctk_nghostzones[1]; j < cctk_lsh[1] - cctk_nghostzones[1]; ++j)
      for (int i = cctk_nghostzones[0]; i < cctk_lsh[0] - cctk_nghostzones[0]; ++i) {
        const int p = CCTK_GFINDEX3D(cctkGH, i, j, k);
        if ((weight != nullptr && weight[p] <= 0) || !std::isfinite(rho[p]))
          continue;
        const Point position{x[p], y[p], z[p]};
        grid_floor = std::min(grid_floor, rho[p]);
        for (int a = 0; a < 3; ++a) {
          grid_min[a] = std::min(grid_min[a], position[a]);
          grid_max[a] = std::max(grid_max[a], position[a]);
        }
        int nearest = -1;
        CCTK_REAL nearest_distance = std::numeric_limits<CCTK_REAL>::infinity();
        for (const int s : followed) {
          const CCTK_REAL d =
              distance(position, {sf_origin_x[s], sf_origin_y[s], sf_origin_z[s]});
          if (d < nearest_distance) {
            nearest = s;
            nearest_distance = d;
          }
        }
        if (nearest >= 0 && nearest_distance <= radius_limit(nearest, outer_radius[nearest])) {
          slot_peaks[nearest] = std::max(slot_peaks[nearest], rho[p]);
        }
        if (any_waiting && rho[p] > 0 && NSTools::local_maximum(density, p, dj, dk))
          local_maxima.add({rho[p], rho[p], position});
      }
}

extern "C" void NSSurface_Find(CCTK_ARGUMENTS) {
  using namespace NSSurface;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int every = output_cadence(compute_every);
  std::vector<int> slots;
  for (int s = 0; s < nsurfaces && s < max_slots; ++s)
    if (surface_owned(s)) {
      slots.push_back(s);
      sf_active[s] = 1;
      sf_valid[s] = 0; // Never let QLM interpret an old contour as current data.
    }
  if (slots.empty() || every <= 0 || cctk_iteration % every != 0)
    return;

  bool any_waiting = false;
  for (const int s : slots)
    any_waiting |= slot_state[s] == waiting;
  const std::vector<Candidate> candidates =
      any_waiting ? local_maxima.gather(cctkGH) : std::vector<Candidate>{};
  reduce(cctkGH, "maximum", max_slots, slot_peaks.data());
  for (auto &peak : slot_peaks)
    if (peak == lowest)
      peak = unset;
  reduce(cctkGH, "minimum", 1, &grid_floor);
  reduce(cctkGH, "minimum", 1, &finest_spacing);
  const auto grid_refl = grid_reflections();
  reduce(cctkGH, "minimum", 3, grid_min.data());
  reduce(cctkGH, "maximum", 3, grid_max.data());
  const auto mirrored = mirrored_axes(grid_refl);
  for (int a = 0; a < 3; ++a)
    if (mirrored[a])
      grid_min[a] = std::min(grid_min[a], -grid_max[a]);
  auto origin_of = [&](int o) -> Point {
    return {sf_origin_x[o], sf_origin_y[o], sf_origin_z[o]};
  };

  // An origin moved since NSSurface last wrote it comes from a centre
  // provider; NSSurface then follows it without moving it. Followed slots are
  // searched first, so waiting slots are placed on the remaining stars.
  std::array<bool, max_slots> external{};
  for (const int s : slots) {
    external[s] = distance(origin_of(s), {last_origin_x[s], last_origin_y[s],
                                          last_origin_z[s]}) > 1e-12;
    if (external[s])
      slot_state[s] = following;
  }
  std::stable_partition(slots.begin(), slots.end(),
                        [&](int s) { return slot_state[s] != waiting; });

  const int rank = CCTK_MyProc(cctkGH), ranks = CCTK_nProcs(cctkGH);
  for (const int s : slots) {
    CCTK_REAL *const radius_s = sf_radius + s * maxntheta * maxnphi;
    const bool had_surface = outer_radius[s] > 0;
    if (!(radial_step[s] > 0))
      radial_step[s] = max_radius[s] > 0 ? max_radius[s] / lattice_intervals : finest_spacing;
    const CCTK_REAL tolerance = tolerance_fraction * radial_step[s];
    const CCTK_REAL limit_s = radius_limit(s, outer_radius[s]);
    auto record = [&](const SearchResult &result) {
      rho_contour[s] = result.cutoff;
      rho_peak[s] = result.peak;
      search_radius[s] = result.limit;
      clipped_fraction[s] = result.clipped;
    };
    auto fail = [&](const SearchResult &result) {
      record(result);
      sf_valid[s] = -1;
      outer_radius[s] = 0.0;
      std::fill_n(radius_s, maxntheta * maxnphi, 0.0);
    };

    const int nt = sf_ntheta[s], np = sf_nphi[s];
    const int gt = sf_nghoststheta[s], gp = sf_nghostsphi[s];
    const int count = nt * np;
    const int first = std::size_t(count) * rank / ranks;
    const int last = std::size_t(count) * (rank + 1) / ranks;
    std::vector<Point> directions;
    for (int ray = first; ray < last; ++ray) {
      const CCTK_REAL theta = sf_origin_theta[s] + (ray % nt) * sf_delta_theta[s];
      const CCTK_REAL phi = sf_origin_phi[s] + (ray / nt) * sf_delta_phi[s];
      directions.push_back({std::sin(theta) * std::cos(phi),
                            std::sin(theta) * std::sin(phi), std::cos(theta)});
    }
    const std::array<bool, 3> reflected{bool(symmetric_x[s]) || grid_refl[0],
                                        bool(symmetric_y[s]) || grid_refl[1],
                                        bool(symmetric_z[s]) || grid_refl[2]};

    // Other followed stars (optionally only those found on this step) that
    // `origin` is, or is not, joined to above `cutoff`.
    auto others = [&](const Point &origin, CCTK_REAL cutoff, bool joined,
                      bool found_now) {
      std::vector<Point> result;
      for (int o = 0; o < nsurfaces && o < max_slots; ++o) {
        if (o == s || !surface_owned(o) || slot_state[o] != following ||
            !finite(origin_of(o)) || (found_now && sf_valid[o] != 1))
          continue;
        const bool same = distance(origin, origin_of(o)) <= tolerance ||
                          connected(cctkGH, origin, origin_of(o), cutoff);
        if (same == joined)
          result.push_back(origin_of(o));
      }
      return result;
    };

    auto search = [&](const Point &origin, CCTK_REAL peak, CCTK_REAL cutoff,
                      CCTK_REAL start) {
      SearchResult result;
      result.origin = origin;
      result.peak = peak;
      result.cutoff = cutoff;
      if (!(result.cutoff > 0))
        return result;
      // Each ray stops at the search limit, at the grid edge, and at the
      // plane halfway to any separate star in its direction.
      const std::vector<Point> companions = others(origin, result.cutoff, false, false);
      std::vector<CCTK_REAL> limits;
      CCTK_REAL nearest = std::numeric_limits<CCTK_REAL>::infinity();
      for (const Point &n : directions) {
        CCTK_REAL limit = std::min(limit_s, distance_to_grid_edge(origin, n));
        for (const Point &c : companions) {
          const CCTK_REAL along = n[0] * (c[0] - origin[0]) + n[1] * (c[1] - origin[1]) +
                                  n[2] * (c[2] - origin[2]);
          if (along > 0)
            limit = std::min(limit, std::pow(distance(c, origin), 2) / (2 * along));
        }
        limits.push_back(limit);
        nearest = std::min(nearest, limit);
      }
      reduce(cctkGH, "minimum", 1, &nearest);
      result.limit = nearest;
      std::vector<CCTK_REAL> local_radii;
      std::vector<char> clipped;
      if (!find_radii(cctkGH, origin, directions, limits, result.cutoff, start,
                      radial_step[s], tolerance, local_radii, clipped))
        return result;
      CCTK_REAL total = std::count(clipped.begin(), clipped.end(), 1);
      reduce(cctkGH, "sum", 1, &total);
      result.reach = 0.0;
      for (std::size_t ray = 0; ray < clipped.size(); ++ray)
        if (!clipped[ray])
          result.reach = std::max(result.reach, local_radii[ray]);
      reduce(cctkGH, "maximum", 1, &result.reach);
      result.clipped = total / count;
      if (result.clipped > max_clipped_fraction)
        return result;
      gather_radii(cctkGH, count, local_radii, result.radii);
      result.shift = volume_centroid(nt, np, gt, gp, sf_origin_theta[s], sf_origin_phi[s],
                                     sf_delta_theta[s], sf_delta_phi[s], reflected,
                                     result.radii);
      result.ok = true;
      return result;
    };

    // An origin NSSurface owns follows the star's core: it moves to the volume
    // centroid of the half-peak contour, which lies deep inside the star and
    // is unaffected by the envelope or by clipping of the surface.
    auto centred = [&](const Point &origin, CCTK_REAL peak) {
      const SearchResult core = search(origin, peak, 0.5 * peak, 0.0);
      if (!core.ok)
        return origin;
      return Point{origin[0] + core.shift[0], origin[1] + core.shift[1],
                   origin[2] + core.shift[2]};
    };

    // The contour must lie well above the atmosphere (otherwise the atmosphere
    // joins the stars); warned once per slot.
    auto level = [&](CCTK_REAL peak) {
      const CCTK_REAL cutoff = contour_density(s, peak, rho_reference[s]);
      if (cutoff > 0 && cutoff < 10 * grid_floor && !warned_floor[s] && rank == 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "NSSurface[%d]: contour density %g is within a factor 10 of the "
                   "lowest density %g on the grid (atmosphere); raise "
                   "rho_surface_fraction", s, double(cutoff), double(grid_floor));
        warned_floor[s] = true;
      }
      return cutoff;
    };

    const Point origin = origin_of(s);
    SearchResult result;
    if (slot_state[s] == following) {
      // Joined to a star found earlier on this step: the stars have merged
      // (or touch at this level). The slot waits to be placed again on a star
      // that separates, keeping its contour level.
      const CCTK_REAL cutoff = level(slot_peaks[s]);
      if (finite(origin) && cutoff > 0 && !others(origin, cutoff, true, true).empty()) {
        if (!external[s])
          slot_state[s] = waiting;
        fail(result);
        continue;
      }
      // Start just beyond the previous surface, when there is one.
      if (finite(origin))
        result = search(external[s] ? origin : centred(origin, slot_peaks[s]),
                        slot_peaks[s], cutoff, 1.25 * outer_radius[s]);
    } else {
      // Place a waiting slot on the densest maximum not joined to a followed
      // star; maxima below placement_fraction of the densest are not stars.
      for (const Candidate &candidate : candidates) {
        if (candidate.value < placement_fraction * candidates.front().value)
          break;
        Point position = candidate.position;
        for (int a = 0; a < 3; ++a)
          if (reflected[a])
            position[a] = 0.0;
        const CCTK_REAL cutoff = level(candidate.value);
        if (candidate.value > cutoff && others(position, cutoff, true, false).empty()) {
          result = search(centred(position, candidate.value), candidate.value, cutoff, 0.0);
          if (result.ok)
            break;
        }
      }
    }

    if (!result.ok) {
      fail(result);
      // Warn when a followed star's surface is lost, or missing at the start.
      if (slot_state[s] == following && (had_surface || cctk_iteration == 0) &&
          rank == 0)
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "NSSurface[%d]: lost the rho=%g contour about (%g,%g,%g); "
                   "nearest limit %g, clipped fraction %g",
                   s, double(result.cutoff), double(origin[0]), double(origin[1]),
                   double(origin[2]), double(result.limit), double(result.clipped));
      continue;
    }

    record(result);
    sf_origin_x[s] = last_origin_x[s] = result.origin[0];
    sf_origin_y[s] = last_origin_y[s] = result.origin[1];
    sf_origin_z[s] = last_origin_z[s] = result.origin[2];
    if (rho_reference[s] <= 0)
      rho_reference[s] = result.peak;
    outer_radius[s] = result.reach;
    slot_state[s] = following;
    for (int j = 0; j < np; ++j)
      for (int i = 0; i < nt; ++i)
        radius_s[i + maxntheta * j] = result.radii[i + nt * j];
    sf_valid[s] = 1;
  }
}
