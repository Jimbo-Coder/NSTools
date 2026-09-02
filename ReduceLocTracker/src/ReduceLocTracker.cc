#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Parameter.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace {

enum class Operation { Maximum, Minimum, MaximumAbs, MinimumAbs };

enum class SelectionRegion { Global, NearPoint, FollowPrevious };

struct Target {
  std::string name;
  Operation operation;
  SelectionRegion selection_region;
  int varindex;
  CCTK_REAL target_x0;
  CCTK_REAL target_y0;
  CCTK_REAL target_z0;
  CCTK_REAL target_radius2;
};

std::vector<Target> targets;
std::vector<CCTK_REAL> local_compare;
std::vector<CCTK_REAL> local_value_sum;
std::vector<CCTK_REAL> local_compare_sum;
std::vector<CCTK_REAL> local_x_sum;
std::vector<CCTK_REAL> local_y_sum;
std::vector<CCTK_REAL> local_z_sum;
std::vector<CCTK_REAL> local_count;

bool initialised = false;
int carpet_reduce_weight_varindex = -2;
bool warned_missing_carpet_reduce_weight = false;
bool warned_null_carpet_reduce_weight = false;
constexpr int max_tracked_slots = 16;

bool is_maximum(Operation operation) {
  return operation == Operation::Maximum || operation == Operation::MaximumAbs;
}

bool uses_absolute(Operation operation) {
  return operation == Operation::MaximumAbs || operation == Operation::MinimumAbs;
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

SelectionRegion parse_selection_region(const char *name) {
  if (std::strcmp(name, "global") == 0) {
    return SelectionRegion::Global;
  }
  if (std::strcmp(name, "near_point") == 0) {
    return SelectionRegion::NearPoint;
  }
  if (std::strcmp(name, "follow_previous") == 0) {
    return SelectionRegion::FollowPrevious;
  }
  CCTK_VERROR("Unknown ReduceLocTracker selection_region \"%s\".", name);
}

CCTK_REAL sentinel(Operation operation) {
  const CCTK_REAL large = std::numeric_limits<CCTK_REAL>::max();
  return is_maximum(operation) ? -large : large;
}

bool better(CCTK_REAL value, CCTK_REAL current, Operation operation) {
  return is_maximum(operation) ? value > current : value < current;
}

int parameter_int(const char *name, const char *thorn, int fallback) {
  int type = -1;
  const void *value = CCTK_ParameterGet(name, thorn, &type);
  if (value == nullptr || type != PARAMETER_INT) {
    return fallback;
  }
  return static_cast<int>(*static_cast<const CCTK_INT *>(value));
}

int effective_compute_every(int requested_compute_every) {
  if (requested_compute_every != -2) {
    return requested_compute_every;
  }

  const char *scalar_io_thorns[] = {"CarpetIOScalar", "IOBasic", "IOScalar"};
  for (const char *thorn : scalar_io_thorns) {
    const int scalar_every = parameter_int("outScalar_every", thorn, -2);
    if (scalar_every != -2) {
      return scalar_every;
    }
  }

  const int io_every = parameter_int("out_every", "IOUtil", -1);
  if (io_every != -1) {
    return io_every;
  }
  return parameter_int("out_every", "IO", -1);
}

bool should_run(int cctk_iteration, int requested_compute_every) {
  const int every = effective_compute_every(requested_compute_every);
  return every > 0 && cctk_iteration % every == 0 && !targets.empty();
}

bool has_finite_target(const Target &target) {
  return std::isfinite(target.target_x0) && std::isfinite(target.target_y0) &&
         std::isfinite(target.target_z0);
}

bool inside_selection_region(const Target &target, CCTK_REAL dx, CCTK_REAL dy,
                             CCTK_REAL dz) {
  if (target.selection_region == SelectionRegion::Global) {
    return true;
  }
  return dx * dx + dy * dy + dz * dz <= target.target_radius2;
}

const CCTK_REAL *get_carpet_reduce_weight(const cGH *cctkGH) {
  if (carpet_reduce_weight_varindex == -2) {
    carpet_reduce_weight_varindex = CCTK_VarIndex("CarpetReduce::weight");
  }
  if (carpet_reduce_weight_varindex < 0) {
    if (!warned_missing_carpet_reduce_weight) {
      CCTK_WARN(2, "CarpetReduce::weight is unavailable; ReduceLocTracker is "
                   "continuing without reduction weights.");
      warned_missing_carpet_reduce_weight = true;
    }
    return nullptr;
  }

  const CCTK_REAL *weight = static_cast<const CCTK_REAL *>(
      CCTK_VarDataPtrI(cctkGH, 0, carpet_reduce_weight_varindex));
  if (weight == nullptr && !warned_null_carpet_reduce_weight) {
    CCTK_WARN(2, "CarpetReduce::weight has no data pointer; ReduceLocTracker is "
                 "continuing without reduction weights.");
    warned_null_carpet_reduce_weight = true;
  }
  return weight;
}

void reset_local() {
  for (std::size_t n = 0; n < targets.size(); ++n) {
    local_compare[n] = sentinel(targets[n].operation);
    local_value_sum[n] = 0.0;
    local_compare_sum[n] = 0.0;
    local_x_sum[n] = 0.0;
    local_y_sum[n] = 0.0;
    local_z_sum[n] = 0.0;
    local_count[n] = 0.0;
  }
}

void validate_targets() {
  DECLARE_CCTK_PARAMETERS;

  if (num_tracked_vars > max_tracked_slots) {
    CCTK_VERROR("num_tracked_vars=%d exceeds the compiled ReduceLocTracker "
                "slot count %d.",
                static_cast<int>(num_tracked_vars), max_tracked_slots);
  }

  targets.clear();
  targets.reserve(num_tracked_vars);

  for (int n = 0; n < num_tracked_vars; ++n) {
    if (tracked_var[n] == nullptr || std::strlen(tracked_var[n]) == 0) {
      CCTK_VERROR("tracked_var[%d] is empty but num_tracked_vars=%d.", n,
                  static_cast<int>(num_tracked_vars));
    }
    if (only_positive_x[n] && only_negative_x[n]) {
      CCTK_VERROR("tracked_var[%d] sets both only_positive_x and "
                  "only_negative_x.",
                  n);
    }

    const int varindex = CCTK_VarIndex(tracked_var[n]);
    if (varindex < 0) {
      CCTK_VERROR("Could not find tracked_var[%d]=\"%s\".", n,
                  tracked_var[n]);
    }

    const int group = CCTK_GroupIndexFromVarI(varindex);
    if (group < 0) {
      CCTK_VERROR("Could not find group for tracked_var[%d]=\"%s\".", n,
                  tracked_var[n]);
    }

    if (CCTK_GroupTypeI(group) != CCTK_GF) {
      CCTK_VERROR("tracked_var[%d]=\"%s\" is not a grid function.", n,
                  tracked_var[n]);
    }
    if (CCTK_GroupDimI(group) != 3) {
      CCTK_VERROR("tracked_var[%d]=\"%s\" is not a 3D grid function.", n,
                  tracked_var[n]);
    }
    if (CCTK_VarTypeI(varindex) != CCTK_VARIABLE_REAL) {
      CCTK_VERROR("tracked_var[%d]=\"%s\" is not CCTK_REAL.", n,
                  tracked_var[n]);
    }

    Target target;
    target.name = tracked_var[n];
    target.operation = parse_operation(reduction[n]);
    target.selection_region = parse_selection_region(selection_region[n]);
    target.varindex = varindex;
    target.target_x0 = target_x0[n];
    target.target_y0 = target_y0[n];
    target.target_z0 = target_z0[n];
    target.target_radius2 = target_radius[n] * target_radius[n];
    if (target.selection_region != SelectionRegion::Global) {
      if (target_radius[n] <= 0.0) {
        CCTK_VERROR("tracked_var[%d] uses selection_region=\"%s\" but "
                    "target_radius[%d]=%.17g is not positive.",
                    n, selection_region[n], n,
                    static_cast<double>(target_radius[n]));
      }
      if (!has_finite_target(target)) {
        CCTK_VERROR("tracked_var[%d] uses selection_region=\"%s\" but "
                    "target_x0/y0/z0 are not all finite.",
                    n, selection_region[n]);
      }
    }
    targets.push_back(target);
  }

  local_compare.resize(targets.size());
  local_value_sum.resize(targets.size());
  local_compare_sum.resize(targets.size());
  local_x_sum.resize(targets.size());
  local_y_sum.resize(targets.size());
  local_z_sum.resize(targets.size());
  local_count.resize(targets.size());
  reset_local();
  initialised = true;
}

void ensure_initialised() {
  if (!initialised) {
    validate_targets();
  }
}

} // namespace

extern "C" void ReduceLocTracker_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  validate_targets();
  if (num_tracked_vars == 0) {
    return;
  }
  for (int n = 0; n < max_tracked_slots; ++n) {
    reduce_loc_value[n] = 0.0;
    reduce_loc_comparison_value[n] = 0.0;
    reduce_loc_x[n] = target_x0[n];
    reduce_loc_y[n] = target_y0[n];
    reduce_loc_z[n] = target_z0[n];
    reduce_loc_radius[n] =
        std::sqrt(reduce_loc_x[n] * reduce_loc_x[n] +
                  reduce_loc_y[n] * reduce_loc_y[n] +
                  reduce_loc_z[n] * reduce_loc_z[n]);
    reduce_loc_count[n] = 0.0;
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

extern "C" void ReduceLocTracker_Search(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  ensure_initialised();
  if (!should_run(cctk_iteration, compute_every)) {
    return;
  }

  const int imin = exclude_ghosts ? cctk_nghostzones[0] : 0;
  const int jmin = exclude_ghosts ? cctk_nghostzones[1] : 0;
  const int kmin = exclude_ghosts ? cctk_nghostzones[2] : 0;
  const int imax = exclude_ghosts ? cctk_lsh[0] - cctk_nghostzones[0]
                                  : cctk_lsh[0];
  const int jmax = exclude_ghosts ? cctk_lsh[1] - cctk_nghostzones[1]
                                  : cctk_lsh[1];
  const int kmax = exclude_ghosts ? cctk_lsh[2] - cctk_nghostzones[2]
                                  : cctk_lsh[2];
  const CCTK_REAL *reduction_weight =
      use_carpet_reduce_weight ? get_carpet_reduce_weight(cctkGH) : nullptr;

  for (std::size_t n = 0; n < targets.size(); ++n) {
    const CCTK_REAL *gf = static_cast<const CCTK_REAL *>(
        CCTK_VarDataPtrI(cctkGH, 0, targets[n].varindex));
    if (gf == nullptr) {
      continue;
    }
    const bool follow =
        targets[n].selection_region == SelectionRegion::FollowPrevious;
    const CCTK_REAL center_x = follow ? reduce_loc_x[n] : targets[n].target_x0;
    const CCTK_REAL center_y = follow ? reduce_loc_y[n] : targets[n].target_y0;
    const CCTK_REAL center_z = follow ? reduce_loc_z[n] : targets[n].target_z0;

    for (int k = kmin; k < kmax; ++k) {
      for (int j = jmin; j < jmax; ++j) {
        for (int i = imin; i < imax; ++i) {
          const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
          if (reduction_weight != nullptr && reduction_weight[index] <= 0.0) {
            continue;
          }

          const CCTK_REAL xx = x[index];
          const CCTK_REAL yy = y[index];
          const CCTK_REAL zz = z[index];

          if ((only_positive_x[n] && xx < 0.0) ||
              (only_negative_x[n] && xx > 0.0)) {
            continue;
          }
          if (!inside_selection_region(targets[n], xx - center_x, yy - center_y,
                                       zz - center_z)) {
            continue;
          }

          const CCTK_REAL raw_value = gf[index];
          if (require_finite && !std::isfinite(raw_value)) {
            continue;
          }

          const CCTK_REAL comparison_value =
              uses_absolute(targets[n].operation) ? std::fabs(raw_value)
                                                  : raw_value;
          if (require_finite && !std::isfinite(comparison_value)) {
            continue;
          }

          if (local_count[n] == 0.0 ||
              better(comparison_value, local_compare[n],
                     targets[n].operation)) {
            local_compare[n] = comparison_value;
            local_value_sum[n] = raw_value;
            local_compare_sum[n] = comparison_value;
            local_x_sum[n] = xx;
            local_y_sum[n] = yy;
            local_z_sum[n] = zz;
            local_count[n] = 1.0;
          } else if (comparison_value == local_compare[n]) {
            local_value_sum[n] += raw_value;
            local_compare_sum[n] += comparison_value;
            local_x_sum[n] += xx;
            local_y_sum[n] += yy;
            local_z_sum[n] += zz;
            local_count[n] += 1.0;
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

  const int ntracked = static_cast<int>(targets.size());
  if (ntracked == 0) {
    return;
  }

  std::vector<CCTK_REAL> global_compare(ntracked);
  std::vector<CCTK_REAL> send_sum(6 * ntracked, 0.0);
  std::vector<CCTK_REAL> global_sum(6 * ntracked, 0.0);

  const int max_handle = CCTK_ReductionArrayHandle("maximum");
  const int min_handle = CCTK_ReductionArrayHandle("minimum");
  const int sum_handle = CCTK_ReductionArrayHandle("sum");
  if (max_handle < 0 || min_handle < 0 || sum_handle < 0) {
    CCTK_WARN(0, "Could not obtain required reduction handles.");
  }

  std::vector<CCTK_REAL> max_input(ntracked, 0.0);
  std::vector<CCTK_REAL> max_output(ntracked, 0.0);
  std::vector<CCTK_REAL> min_input(ntracked, 0.0);
  std::vector<CCTK_REAL> min_output(ntracked, 0.0);
  std::vector<CCTK_INT> is_max_op(ntracked, 0);
  std::vector<CCTK_INT> is_min_op(ntracked, 0);

  for (int n = 0; n < ntracked; ++n) {
    if (is_maximum(targets[n].operation)) {
      max_input[n] = local_compare[n];
      min_input[n] = std::numeric_limits<CCTK_REAL>::max();
      is_max_op[n] = 1;
    } else {
      max_input[n] = -std::numeric_limits<CCTK_REAL>::max();
      min_input[n] = local_compare[n];
      is_min_op[n] = 1;
    }
  }

  if (CCTK_ReduceLocArrayToArray1D(cctkGH, -1, max_handle, max_input.data(),
                                   max_output.data(), ntracked,
                                   CCTK_VARIABLE_REAL) ||
      CCTK_ReduceLocArrayToArray1D(cctkGH, -1, min_handle, min_input.data(),
                                   min_output.data(), ntracked,
                                   CCTK_VARIABLE_REAL)) {
    CCTK_WARN(0, "Failed to reduce local comparison extrema.");
  }

  for (int n = 0; n < ntracked; ++n) {
    global_compare[n] = is_max_op[n] ? max_output[n] : min_output[n];

    if (local_count[n] > 0.0 && local_compare[n] == global_compare[n]) {
      send_sum[0 * ntracked + n] = local_value_sum[n];
      send_sum[1 * ntracked + n] = local_compare_sum[n];
      send_sum[2 * ntracked + n] = local_x_sum[n];
      send_sum[3 * ntracked + n] = local_y_sum[n];
      send_sum[4 * ntracked + n] = local_z_sum[n];
      send_sum[5 * ntracked + n] = local_count[n];
    }
  }

  if (CCTK_ReduceLocArrayToArray1D(cctkGH, -1, sum_handle, send_sum.data(),
                                   global_sum.data(), 6 * ntracked,
                                   CCTK_VARIABLE_REAL)) {
    CCTK_WARN(0, "Failed to reduce extrema coordinate sums.");
  }

  for (int n = 0; n < ntracked; ++n) {
    const CCTK_REAL count = global_sum[5 * ntracked + n];
    const bool valid = count > 0.0;
    // An empty search invalidates this sample but retains the last result.
    reduce_loc_count[n] = count;
    if (valid) {
      reduce_loc_value[n] = global_sum[0 * ntracked + n] / count;
      reduce_loc_comparison_value[n] = global_sum[1 * ntracked + n] / count;
      reduce_loc_x[n] = global_sum[2 * ntracked + n] / count;
      reduce_loc_y[n] = global_sum[3 * ntracked + n] / count;
      reduce_loc_z[n] = global_sum[4 * ntracked + n] / count;
      reduce_loc_radius[n] =
          std::sqrt(reduce_loc_x[n] * reduce_loc_x[n] +
                    reduce_loc_y[n] * reduce_loc_y[n] +
                    reduce_loc_z[n] * reduce_loc_z[n]);
      if (verbose) {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "[%d] %s %s value=%+.16e compare=%+.16e "
                   "loc=(%+.16e,%+.16e,%+.16e) radius=%+.16e count=%.0f",
                   n, targets[n].name.c_str(), reduction[n],
                   reduce_loc_value[n], reduce_loc_comparison_value[n],
                   reduce_loc_x[n], reduce_loc_y[n], reduce_loc_z[n],
                   reduce_loc_radius[n], reduce_loc_count[n]);
      }
    }
  }
}
