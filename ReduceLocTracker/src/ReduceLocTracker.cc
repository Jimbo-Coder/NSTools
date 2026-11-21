#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameter.h"
#include "cctk_Parameters.h"

#include "NSTools.hxx"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace {
using namespace NSTools;

enum class Operation { Maximum, Minimum, MaximumAbs, MinimumAbs };

constexpr int max_tracked_slots = 16;

struct Target {
  int slot;
  std::string name;
  Operation operation;
  int varindex;
  Point seed;
  CCTK_REAL target_radius2; // infinite: no limit
  // Seed at the origin: placed by the tracker on a local extremum, and placed
  // again after it is lost.
  bool placed;
};

// Slots tracking one variable with one reduction share the grid.
struct Family {
  int varindex;
  Operation operation;
  std::vector<int> members; // indices into targets, in slot order
};

// Ties at the best value are averaged.
struct Accumulator {
  CCTK_REAL best;
  CCTK_REAL value_sum;
  CCTK_REAL compare_sum;
  CCTK_REAL x_sum;
  CCTK_REAL y_sum;
  CCTK_REAL z_sum;
  CCTK_REAL count;
};

std::vector<Target> targets;
std::vector<Family> families;
std::vector<Accumulator> local_acc;          // per target
std::vector<Extrema> local_extrema; // per family

bool initialised = false;
bool warned_segment = false;

// The value each operation maximises.
CCTK_REAL tracked(Operation operation, CCTK_REAL value) {
  switch (operation) {
  case Operation::Maximum:
    return value;
  case Operation::Minimum:
    return -value;
  case Operation::MaximumAbs:
    return std::fabs(value);
  case Operation::MinimumAbs:
    return -std::fabs(value);
  }
  return value;
}

// The value reported as reduce_loc_comparison_value.
CCTK_REAL comparison(Operation operation, CCTK_REAL value) {
  return operation == Operation::MaximumAbs || operation == Operation::MinimumAbs
             ? std::fabs(value)
             : value;
}

Operation parse_operation(const char *name) {
  if (std::strcmp(name, "maximum") == 0) {
    return Operation::Maximum;
  }
  if (std::strcmp(name, "minimum") == 0) {
    return Operation::Minimum;
  }
  if (std::strcmp(name, "maximum_abs") == 0) {
    return Operation::MaximumAbs;
  }
  if (std::strcmp(name, "minimum_abs") == 0) {
    return Operation::MinimumAbs;
  }
  CCTK_VERROR("Unknown ReduceLocTracker reduction \"%s\".", name);
}

bool should_run(int cctk_iteration, int requested_compute_every) {
  const int every = output_cadence(requested_compute_every);
  return every > 0 && cctk_iteration % every == 0 && !targets.empty();
}

void reset_local() {
  const CCTK_REAL none = -std::numeric_limits<CCTK_REAL>::max();
  for (auto &acc : local_acc) {
    acc = {none, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  }
  for (auto &extrema : local_extrema) {
    extrema.clear();
  }
}

void update_accumulator(Accumulator &acc, CCTK_REAL value, CCTK_REAL compare,
                        CCTK_REAL best, const Point &p) {
  if (acc.count == 0.0 || best > acc.best) {
    acc = {best, value, compare, p[0], p[1], p[2], 1.0};
  } else if (best == acc.best) {
    acc.value_sum += value;
    acc.compare_sum += compare;
    acc.x_sum += p[0];
    acc.y_sum += p[1];
    acc.z_sum += p[2];
    acc.count += 1.0;
  }
}

void add_target(int slot, const char *name, Operation operation, const Point &seed,
                CCTK_REAL radius) {
  const int varindex = CCTK_VarIndex(name);
  if (varindex < 0) {
    CCTK_VERROR("Could not find tracked_var[%d]=\"%s\".", slot, name);
  }
  const int group = CCTK_GroupIndexFromVarI(varindex);
  if (group < 0 || CCTK_GroupTypeI(group) != CCTK_GF || CCTK_GroupDimI(group) != 3 ||
      CCTK_VarTypeI(varindex) != CCTK_VARIABLE_REAL) {
    CCTK_VERROR("tracked_var[%d]=\"%s\" is not a real 3D grid function.", slot, name);
  }
  if (!finite(seed)) {
    CCTK_VERROR("tracked_var[%d] has non-finite target_x0/y0/z0.", slot);
  }
  Target target;
  target.slot = slot;
  target.name = name;
  target.operation = operation;
  target.varindex = varindex;
  target.seed = seed;
  target.target_radius2 =
      radius > 0 ? radius * radius : std::numeric_limits<CCTK_REAL>::infinity();
  target.placed = seed[0] == 0.0 && seed[1] == 0.0 && seed[2] == 0.0;
  targets.push_back(target);
}

void validate_targets() {
  DECLARE_CCTK_PARAMETERS;

  targets.clear();
  for (int n = 0; n < max_tracked_slots; ++n) {
    if (tracked_var[n] != nullptr && std::strlen(tracked_var[n]) > 0) {
      add_target(n, tracked_var[n], parse_operation(reduction[n]),
                 {target_x0[n], target_y0[n], target_z0[n]}, target_radius[n]);
    }
  }
  // Nothing configured: follow up to three density maxima (stars).
  if (targets.empty() && CCTK_VarIndex("HydroBase::rho") >= 0) {
    for (int n = 0; n < 3; ++n) {
      add_target(n, "HydroBase::rho", Operation::Maximum, {0.0, 0.0, 0.0}, 0.0);
    }
  }

  families.clear();
  for (int n = 0; n < static_cast<int>(targets.size()); ++n) {
    auto family = std::find_if(families.begin(), families.end(), [&](const Family &f) {
      return f.varindex == targets[n].varindex && f.operation == targets[n].operation;
    });
    if (family == families.end()) {
      families.push_back({targets[n].varindex, targets[n].operation, {}});
      family = families.end() - 1;
    }
    family->members.push_back(n);
  }

  local_acc.resize(targets.size());
  local_extrema.resize(families.size());
  reset_local();
  initialised = true;
}

void ensure_initialised() {
  if (!initialised) {
    validate_targets();
  }
}

// reduce_loc_count < 0: seeded, not searched yet; 0: lost, or waiting to be
// placed; > 0: valid, with that many tied points.
bool waiting(const Target &target, const CCTK_REAL *count) {
  return target.placed && count[target.slot] == 0.0;
}

// Slots whose nearest points are their own: valid, or seeded and not yet
// searched. A slot that lost its extremum claims none, so it cannot hide a
// merged peak from the slot still following it.
bool claims(const Target &target, const CCTK_REAL *count) {
  return count[target.slot] != 0.0;
}

// Reduce the per-process best values and their tied points.
void reduce_accumulators(const cGH *cctkGH, std::vector<Accumulator> &acc) {
  const int m = static_cast<int>(acc.size());
  if (m == 0) {
    return;
  }
  std::vector<CCTK_REAL> best(m), sums(6 * m, 0.0);
  for (int i = 0; i < m; ++i) {
    best[i] = acc[i].best;
  }
  reduce(cctkGH, "maximum", m, best.data());
  for (int i = 0; i < m; ++i) {
    if (acc[i].count > 0.0 && acc[i].best == best[i]) {
      const CCTK_REAL mine[6] = {acc[i].value_sum, acc[i].compare_sum, acc[i].x_sum,
                                 acc[i].y_sum,     acc[i].z_sum,       acc[i].count};
      for (int c = 0; c < 6; ++c) {
        sums[c * m + i] = mine[c];
      }
    }
  }
  reduce(cctkGH, "sum", 6 * m, sums.data());
  for (int i = 0; i < m; ++i) {
    acc[i] = {best[i],          sums[0 * m + i], sums[1 * m + i], sums[2 * m + i],
              sums[3 * m + i],  sums[4 * m + i], sums[5 * m + i]};
  }
}

// Two extrema of one family are separate when its tracked value dips
// between them, along the straight segment, by at least min_prominence times
// the stronger one's magnitude. Collective.
bool separate(const cGH *cctkGH, const Family &family, const Point &a,
              CCTK_REAL value_a, const Point &b, CCTK_REAL value_b) {
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL lowest;
  const Operation operation = family.operation;
  if (!segment_minimum(cctkGH, family.varindex, a, b, segment_samples,
                                [operation](CCTK_REAL v) { return tracked(operation, v); },
                                lowest)) {
    if (!warned_segment) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Cannot interpolate %s between tracked extrema (needs CarpetInterp "
                 "and AEILocalInterp); treating them as separate.",
                 CCTK_FullVarName(family.varindex));
      warned_segment = true;
    }
    return true;
  }
  const CCTK_REAL stronger = std::max(value_a, value_b);
  const CCTK_REAL weaker = std::min(value_a, value_b);
  return weaker - lowest >= min_prominence * std::fabs(stronger);
}

void publish(const Target &target, const Accumulator &acc, CCTK_REAL *value,
             CCTK_REAL *compare, CCTK_REAL *x, CCTK_REAL *y, CCTK_REAL *z,
             CCTK_REAL *radius, CCTK_REAL *count) {
  const int slot = target.slot;
  value[slot] = acc.value_sum / acc.count;
  compare[slot] = acc.compare_sum / acc.count;
  x[slot] = acc.x_sum / acc.count;
  y[slot] = acc.y_sum / acc.count;
  z[slot] = acc.z_sum / acc.count;
  radius[slot] = std::sqrt(x[slot] * x[slot] + y[slot] * y[slot] + z[slot] * z[slot]);
  count[slot] = acc.count;
}

} // namespace

extern "C" void ReduceLocTracker_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  validate_targets();
  for (int n = 0; n < max_tracked_slots; ++n) {
    reduce_loc_value[n] = 0.0;
    reduce_loc_comparison_value[n] = 0.0;
    reduce_loc_x[n] = target_x0[n];
    reduce_loc_y[n] = target_y0[n];
    reduce_loc_z[n] = target_z0[n];
    reduce_loc_radius[n] = std::sqrt(reduce_loc_x[n] * reduce_loc_x[n] +
                                     reduce_loc_y[n] * reduce_loc_y[n] +
                                     reduce_loc_z[n] * reduce_loc_z[n]);
    reduce_loc_count[n] = 0.0;
  }
  for (const Target &t : targets) {
    if (!t.placed) {
      reduce_loc_count[t.slot] = -1.0; // seeded, not searched yet
    }
  }
}

extern "C" void ReduceLocTracker_Begin(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  ensure_initialised();
  if (!should_run(cctk_iteration, compute_every)) {
    return;
  }
  reset_local();
}

// One pass per family. Each point belongs to the nearest claiming slot
// (within its radius); a slot that lost its extremum searches the points
// nearer to it than to any claiming slot. While a member waits to be placed,
// the family's local extrema are collected.
extern "C" void ReduceLocTracker_Search(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  ensure_initialised();
  if (!should_run(cctk_iteration, compute_every)) {
    return;
  }

  static const int weight_index = CCTK_VarIndex("CarpetReduce::weight");
  const CCTK_REAL *weight =
      weight_index >= 0
          ? static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, weight_index))
          : nullptr;
  const int dj = cctk_ash[0], dk = cctk_ash[0] * cctk_ash[1];

  for (std::size_t f = 0; f < families.size(); ++f) {
    const Family &family = families[f];
    const CCTK_REAL *gf =
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, family.varindex));
    if (gf == nullptr) {
      continue;
    }
    std::vector<int> claiming, lost;
    bool collect = false;
    for (const int n : family.members) {
      const Target &t = targets[n];
      if (waiting(t, reduce_loc_count)) {
        collect = true;
      } else if (claims(t, reduce_loc_count)) {
        claiming.push_back(n);
      } else {
        lost.push_back(n);
      }
    }
    auto location = [&](int n) {
      const int slot = targets[n].slot;
      return Point{reduce_loc_x[slot], reduce_loc_y[slot], reduce_loc_z[slot]};
    };
    std::vector<Point> claim_at, lost_at;
    for (const int n : claiming) {
      claim_at.push_back(location(n));
    }
    for (const int n : lost) {
      lost_at.push_back(location(n));
    }
    auto value = [&](int q) { return tracked(family.operation, gf[q]); };

    for (int k = cctk_nghostzones[2]; k < cctk_lsh[2] - cctk_nghostzones[2]; ++k) {
      for (int j = cctk_nghostzones[1]; j < cctk_lsh[1] - cctk_nghostzones[1]; ++j) {
        for (int i = cctk_nghostzones[0]; i < cctk_lsh[0] - cctk_nghostzones[0]; ++i) {
          const int p = CCTK_GFINDEX3D(cctkGH, i, j, k);
          if ((weight != nullptr && weight[p] <= 0.0) || !std::isfinite(gf[p])) {
            continue;
          }
          const Point position{x[p], y[p], z[p]};
          const CCTK_REAL v = value(p);

          if (collect && i > 0 && j > 0 && k > 0 && i < cctk_lsh[0] - 1 &&
              j < cctk_lsh[1] - 1 && k < cctk_lsh[2] - 1 &&
              local_maximum(value, p, dj, dk)) {
            local_extrema[f].add({v, gf[p], position});
          }

          CCTK_REAL nearest = std::numeric_limits<CCTK_REAL>::infinity();
          int owner = -1;
          for (std::size_t c = 0; c < claiming.size(); ++c) {
            const CCTK_REAL d2 = distance2(position, claim_at[c]);
            if (d2 < nearest) {
              nearest = d2;
              owner = claiming[c];
            }
          }
          if (owner >= 0 && nearest <= targets[owner].target_radius2) {
            update_accumulator(local_acc[owner], gf[p], comparison(family.operation, gf[p]),
                               v, position);
          }
          for (std::size_t l = 0; l < lost.size(); ++l) {
            const CCTK_REAL d2 = distance2(position, lost_at[l]);
            if (d2 <= nearest && d2 <= targets[lost[l]].target_radius2) {
              update_accumulator(local_acc[lost[l]], gf[p],
                                 comparison(family.operation, gf[p]), v, position);
            }
          }
        }
      }
    }
  }
}

extern "C" void ReduceLocTracker_Reduce(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  ensure_initialised();
  if (!should_run(cctk_iteration, compute_every)) {
    return;
  }

  std::vector<Accumulator> acc = local_acc;
  reduce_accumulators(cctkGH, acc);

  for (std::size_t f = 0; f < families.size(); ++f) {
    const Family &family = families[f];
    // Searching slots first, in slot order: one whose extremum is not
    // separate from one found earlier this step has lost it (the extrema
    // merged) and keeps its last location with count 0. Then waiting slots
    // are placed on the strongest separate local extrema.
    std::vector<int> searched, placing;
    for (const int n : family.members) {
      (waiting(targets[n], reduce_loc_count) ? placing : searched).push_back(n);
    }
    std::vector<std::pair<Point, CCTK_REAL>> found;
    auto separate_from_found = [&](const Point &p, CCTK_REAL v) {
      for (const auto &other : found) {
        if (!separate(cctkGH, family, p, v, other.first, other.second)) {
          return false;
        }
      }
      return true;
    };

    for (const int n : searched) {
      const Target &t = targets[n];
      const Accumulator &a = acc[n];
      if (a.count <= 0.0) {
        reduce_loc_count[t.slot] = 0.0;
        continue;
      }
      const Point p{a.x_sum / a.count, a.y_sum / a.count, a.z_sum / a.count};
      if (!separate_from_found(p, a.best)) {
        reduce_loc_count[t.slot] = 0.0;
        continue;
      }
      publish(t, a, reduce_loc_value, reduce_loc_comparison_value, reduce_loc_x,
              reduce_loc_y, reduce_loc_z, reduce_loc_radius, reduce_loc_count);
      found.push_back({p, a.best});
    }

    if (placing.empty()) {
      continue;
    }
    const std::vector<Extremum> candidates = local_extrema[f].gather(cctkGH);
    std::size_t next = 0;
    for (const int n : placing) {
      const Target &t = targets[n];
      for (; next < candidates.size(); ++next) {
        const Extremum &c = candidates[next];
        if (!separate_from_found(c.position, c.value)) {
          continue;
        }
        const Accumulator a{c.value,       c.raw,         comparison(family.operation, c.raw),
                            c.position[0], c.position[1], c.position[2],
                            1.0};
        publish(t, a, reduce_loc_value, reduce_loc_comparison_value, reduce_loc_x,
                reduce_loc_y, reduce_loc_z, reduce_loc_radius, reduce_loc_count);
        found.push_back({c.position, c.value});
        ++next;
        break;
      }
    }
  }
}
