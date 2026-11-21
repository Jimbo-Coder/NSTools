# ReduceLocTracker

`ReduceLocTracker` follows the locations of extrema of real 3D grid
functions: density maxima of stars, lapse minima, or any variable named by
its Cactus name (the thorn does not inherit from the thorn that owns it).

## Minimal setup

Activate `ReduceLocTracker` (with `CarpetInterp` and `AEILocalInterp`) and
output its scalars:

```ccl
IOASCII::out0D_vars = "ReduceLocTracker::reduce_loc_reals"
```

With no other parameters, slots 0-2 follow the density maxima of
`HydroBase::rho`: each slot places itself on a separate star, and spare slots
wait with `reduce_loc_count = 0`. To follow something else, name it per slot;
the slots used are those with a `tracked_var`:

```ccl
ReduceLocTracker::tracked_var[0] = "HydroBase::rho"   # two stars
ReduceLocTracker::tracked_var[1] = "HydroBase::rho"
ReduceLocTracker::tracked_var[2] = "ADMBase::alp"     # and their lapse minima
ReduceLocTracker::reduction[2]   = "minimum"
ReduceLocTracker::tracked_var[3] = "ADMBase::alp"
ReduceLocTracker::reduction[3]   = "minimum"
```

`reduction` is `maximum` (default), `minimum`, `maximum_abs` or
`minimum_abs`. Optional per slot:

- `target_x0/y0/z0`: a starting location (for example the stars' centres
  from initial data, which fixes which star is slot 0). Left at 0, the
  tracker places the slot itself.
- `target_radius`: a search radius about the previous location (default 0:
  no limit beyond the points nearest this slot).

`compute_every = -2` (default) updates on the 0D/scalar output cadence
(else `IO::out_every`), as NSSurface does.

## Outputs

`reduce_loc_reals` holds, per slot, `reduce_loc_value` (the variable at the
extremum), `reduce_loc_comparison_value` (the value compared: `|value|` for
the absolute reductions), `reduce_loc_x/y/z`, `reduce_loc_radius` (distance
from the coordinate origin) and `reduce_loc_count`. A count > 0 is a valid
result for the latest search (the number of tied points, whose location and
value are averaged). A count of 0 means the slot has no extremum of its own
at this search (merged, lost, or waiting to be placed); the last valid
location and value are kept. The group is checkpointed: restarts continue
exactly as an uninterrupted run.

## How slots behave

Slots tracking the same variable with the same reduction share the grid:
each grid point belongs to the slot whose last valid location is nearest, so
two slots cannot drift onto the same extremum.

Two extrema are separate when the tracked value (the variable, its negative
for `minimum`, or `+-|value|`) dips along the straight segment between them by
at least `min_prominence` (default 0.05) times the stronger one's magnitude.
Two stars dip to the atmosphere; lapse minima of two neutron stars dip by
10-30% until they nearly touch; grid noise in a stellar core or gauge noise
near the outer boundary dips far less and is never taken for an extremum.
When a slot's extremum is not separate from one found earlier on the same
step (by a lower slot), the extrema have merged: the slot reports count 0
and claims no points until it has a separate extremum again.

A slot without a starting location is placed on the strongest local extremum
(above or below all 26 grid neighbours) separate from every extremum already
tracked, and is placed again after it loses its extremum. Placement does not
depend on where the extrema lie relative to the coordinates. Through a close
encounter or collision, which slot holds which object can change.

The separation test interpolates the variable (linear Lagrange). Without
CarpetInterp and AEILocalInterp the tracker warns once and treats extrema as
separate. The local-extremum and segment code is shared with NSSurface
through `NSTools.hxx`.
