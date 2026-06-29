# ReduceLocTracker

`ReduceLocTracker` tracks the coordinate location of extrema of selected real
3D grid functions.  Each tracked variable is selected by its Cactus variable
name, so the thorn does not need to inherit from the physics thorn that owns
the target variable.

Example:

```ccl
ActiveThorns = "... ReduceLocTracker ..."

ReduceLocTracker::num_tracked_vars = 2
ReduceLocTracker::tracked_var[0] = "HydroBase::rho"
ReduceLocTracker::tracked_var[1] = "GRHayLHD::rho_star"
ReduceLocTracker::reduction[0] = "maximum"
ReduceLocTracker::reduction[1] = "maximum_abs"
ReduceLocTracker::use_carpet_reduce_weight = "yes"

IOScalar::outScalar_vars = "ReduceLocTracker::reduce_loc_value ReduceLocTracker::reduce_loc_x ReduceLocTracker::reduce_loc_y ReduceLocTracker::reduce_loc_z ReduceLocTracker::reduce_loc_radius ReduceLocTracker::reduce_loc_count ReduceLocTracker::reduce_loc_valid"
```

The thorn outputs scalar arrays indexed in the same order as
`tracked_var[]`.  `reduce_loc_value` is the raw grid-function value at the
selected location, while `reduce_loc_comparison_value` is the value used for the
selection.  These are different for absolute-value reductions.  If more than
one point has the selected extremum, the reported location and value are the
average over tied points and `reduce_loc_count` records the number of ties.
`reduce_loc_radius` is the coordinate radius of the reported location.

By default, `compute_every = -2` updates on the scalar-output cadence:
`IOScalar::outScalar_every` when available, otherwise `IO::out_every`.  Set a
positive `compute_every` explicitly if the tracker should update more often
than scalar output, for example for Trigger logic.

For moving multi-peak data, such as two neutron-star density maxima, each slot
can restrict its search region:

```ccl
ReduceLocTracker::tracked_var[0] = "HydroBase::rho"
ReduceLocTracker::reduction[0] = "maximum"
ReduceLocTracker::selection_region[0] = "follow_previous"
ReduceLocTracker::target_x0[0] = 12.0
ReduceLocTracker::target_y0[0] = 0.0
ReduceLocTracker::target_z0[0] = 0.0
ReduceLocTracker::target_radius[0] = 8.0

ReduceLocTracker::tracked_var[1] = "HydroBase::rho"
ReduceLocTracker::reduction[1] = "maximum"
ReduceLocTracker::selection_region[1] = "follow_previous"
ReduceLocTracker::target_x0[1] = -12.0
ReduceLocTracker::target_y0[1] = 0.0
ReduceLocTracker::target_z0[1] = 0.0
ReduceLocTracker::target_radius[1] = 8.0
```

`selection_region = "near_point"` always searches around the fixed target
point given by `target_x0/y0/z0`.  `selection_region = "follow_previous"` uses
`target_x0/y0/z0` only as the initial seed, then recentres the next search on
the previous valid reduction location.  `target_radius` is constant for both
modes.  During merger, two follow-previous slots can still converge to the same
physical maximum once the peaks are no longer separated.

In Carpet runs, `use_carpet_reduce_weight = "yes"` makes the local search skip
points with `CarpetReduce::weight <= 0`, matching the active region used by
Carpet reductions.  If `CarpetReduce::weight` is not present, the thorn warns
once and continues without this mask.
