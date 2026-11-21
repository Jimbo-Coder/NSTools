# NSSurface

Find a neutron star's `HydroBase::rho` density contour and publish it to
SphericalSurface. QuasiLocalMeasures consumes the shape and calculates proper
area, areal radius, circumferences and spin. It does not use VolumeIntegrals,
and needs ReduceLocTracker only in the ThornList (for shared code).

The origin is read from `SphericalSurface::sf_origin_*` before each search.
At startup it is initialized from SphericalSurface's `origin_x/y/z` parameters.
An existing centre provider can update it before `NSSurface_Find`; the finder
is ordered after `VI_GRMHD_VolumeIntegralGroup` and before `qlm_calculate`.
NSSurface then follows that origin and never moves it.

Set `find_surface[s] = yes` on one slot per star. A slot whose origin is left
at `(0, 0, 0)` is placed by the finder on the densest local density maximum
(at least all 26 grid neighbours) that is not joined to another slot's star
and is at least `placement_fraction` (0.01) of the densest maximum; from then
on, before each search, its origin is moved to the volume centroid of the
star's core (the contour at half its current peak density).
This does not depend on where the stars are relative to the coordinates.
Several maxima inside one star are all joined to it and are skipped.

Two slots for a binary are enough: when the stars merge, the later slot waits
(`sf_valid = -1`) and the earlier one follows the remnant; if the stars
separate again, as in a close encounter, the waiting slot is placed back on
the separate star. A slot with nothing to follow stays invalid. By default,
`compute_every = -2` matches the 0D/scalar output or `QuasiLocalMeasures`
cadence.

## Minimal setup

Add `NSTools/NSSurface` and `NSTools/ReduceLocTracker` to the ThornList
(ReduceLocTracker need not be active). Activate `NSSurface`,
`SphericalSurface`, `CarpetInterp` and `AEILocalInterp` with the usual grid
and hydro thorns. For a binary:

```ccl
SphericalSurface::nsurfaces = 3       # two stars and a spare
NSSurface::find_surface[0] = yes
NSSurface::find_surface[1] = yes
NSSurface::find_surface[2] = yes

IOASCII::out0D_vars = "SphericalSurface::sf_valid SphericalSurface::sf_origin
                       SphericalSurface::sf_info NSSurface::search_info"
IOASCII::out2D_vars = "SphericalSurface::sf_radius"   # the shapes, if wanted
```

Nothing else needs choosing. The slots place themselves on the stars, the
contour is 1e-4 of each star's initial peak density (a warning is printed
if that is within a factor 10 of the lowest density on the grid, the
atmosphere), the search limit follows each star
(`max_radius = 0`), the sampling step is the finest grid spacing, and the
cadence is the 0D/scalar output cadence. SphericalSurface's default angular
grid (`ntheta = 19`, `nphi = 38` including ghosts) is used unless
`ntheta[s]`/`nphi[s]` are raised.

To have QuasiLocalMeasures compute the area, areal radius and spin of a
surface, activate `QuasiLocalMeasures` and `TGRtensor` and name the slot:

```ccl
SphericalSurface::name[0] = "star 1"
QuasiLocalMeasures::num_surfaces = 1
QuasiLocalMeasures::surface_name[0] = "star 1"
QuasiLocalMeasures::spatial_order = 4
QuasiLocalMeasures::interpolator_options = "order=4"
```

Use the same cadence for any centre provider, HydroBase/ADMBase conversion,
the surface finder and QLM. Up to 42 slots are available, shared with other
SphericalSurface users; several contours can share one origin.

## Symmetry

Cartesian grid symmetries are handled by Carpet's registered symmetry
interpolation, including reflection and rotational symmetry. The origin and
search ball can therefore extend across a symmetry boundary.

NSSurface also accepts SphericalSurface's reduced angular grids. However,
QuasiLocalMeasures itself requires a full angular sphere: keep
`SphericalSurface::symmetric_x/y/z = no` for slots consumed by QLM. This does
not require evolving a full Cartesian domain. Keep `auto_res = no`; explicit
angular resolution makes the surface layout consistent across restarts.

## Contour level

`rho_surface_mode` chooses the contour density:

- `initial` (default): `rho_surface_fraction[s]` times the star's peak density
  at its first search. The level then stays fixed, so the surface moves only
  when the matter does.
- `current`: `rho_surface_fraction[s]` times the star's peak density at each
  search. A falling peak moves this contour outward without any matter moving.
- `absolute`: `rho_surface[s]`.

A star's peak is the largest density on the grid closer to its origin than to
any other followed star's origin, and within its search limit. Both the contour
density and the peak are written out at every search (see Outputs).

## Search and validity

The origin must be inside the contour. The finder samples each angular ray on
a fixed radial lattice, brackets the first downward crossing and bisects it
to 1/4096 of the lattice step using linear density interpolation (no
high-order overshoots near the atmosphere). The step is the finest grid
spacing at the slot's first search, or `max_radius/64` when `max_radius` is
set; it is kept for the rest of the run.

The first pass covers the lattice up to 1.25 times the previous surface's
largest radius; rays still inside the star are extended until they cross or
reach their limit. The lattice is fixed, so the result does not depend on the
previous surface, only the cost does: restarts reproduce an uninterrupted run.

Each ray stops at the search limit, at the edge of the grid (extended across
reflection and rotating symmetries), and, for rays heading toward another
followed star, at the plane halfway between the two origins. The search limit
is `max_radius` when set; by default it is twice the star's largest
unclipped radius at the previous search (none before the first surface), so
it follows a star that expands, and keeps rays into a tidal tail from
running to the grid edge.
A ray with no crossing before its limit is clipped: its radius is set to the
limit. The surface stays valid while the clipped fraction of rays is at most
`max_clipped_fraction` (default 0.25; 0 requires a crossing on every ray).
The core used for re-centring lies well inside these limits, so tidal tails
and clipping do not pull the origin.

Two origins are the same body when the density stays above the contour level
along the straight segment between them. A star joined in this way to a star
found earlier on the same step has merged with it, or touches it at this
level: its slot waits to be placed again on a star that separates, keeping
the contour level from its first placement (so its time series stays on one
level), unless its origin comes from a centre provider, in which case it is
only invalid for that step. The halfway plane applies only between stars
that are not joined.

This representation is intended for star-shaped stars. It does not describe
detached ejecta, overhangs or arbitrary merger topology. Features narrower
than a lattice interval can be missed. Check radial/angular resolution and
contour sensitivity when measuring oscillations.

A warning is printed when a followed star's surface is lost, or cannot be
found at the first iteration; `search_info` shows why.

Only selected surface slots are modified. `sf_valid` is positive only after
a successful search, negative on failure, and zero on skipped iterations.
Failed searches zero the radii. For NSSurface's slots, auxiliary fields have
their own `NSSurface::info_valid` flag:
1 means current information, 0 skipped or unused, and -1 unavailable. Check it
when reading `sf_info`; those fields are NaN on skipped or failed updates.
An auxiliary failure does not change the contour's `sf_valid` or prevent QLM
from processing it. QLM still requires valid metric data for its own calculation.
All active angular points, including angular ghosts, are
searched; padding is initialized separately.

## Restarts

Slot states, the initial peak density, the last surface's largest radius and
the last origin NSSurface wrote are kept in checkpointed grid scalars.
Checkpoints written before these existed are recovered with slots following
their stored origins and the contour level taken again at the first search.

## Parallelism

Angular rays are divided among MPI ranks. Each rank requests interpolation for
its own rays; the completed shape is then replicated for SphericalSurface and
QLM. Search failures are shared before proceeding to the next collective call.

Threading inside interpolation follows the driver's normal Cactus/OpenMP
settings. NSSurface's local point loops are serial; ray work is distributed
through MPI. The finder adds no thread regions or threading parameters. The
auxiliary geometry calculation also distributes its interpolation and point
work among MPI ranks, then reduces the scalar moments and bounds.

## Outputs

Use `QuasiLocalMeasures::qlm_area` and `qlm_radius` for proper area and
`sqrt(area/(4*pi))`. Plot `qlm_radius(t)/qlm_radius(0)` for relative size, and
retain coordinate radii to inspect deformation. The finder supplies
`sf_radius`, `sf_origin`, `sf_active` and `sf_valid`. SphericalSurface supplies
the angular descriptors; QLM derives its diagnostics from this geometry.

`NSSurface::search_info` holds, per slot, the last search's `rho_contour`
(contour density), `rho_peak` (the star's peak density), `search_radius`
(the nearest ray limit) and `clipped_fraction`. They are also set on failed
searches where known, which shows why a search failed.

The separate `NSSurface_UpdateInfo` routine fills all of `sf_info` after each
successful search: proper area, mean coordinate radius, surface centroid,
central quadrupole moments, minimum/maximum coordinate radius and bounding box.
The centroid is a surface-area-weighted coordinate position, not a barycentre;
it does not replace or move the supplied search origin.

Area and moments use the physical `ADMBase::metric`, linearly interpolated at
the surface points. Angular radius derivatives use centred fourth-order
differences throughout. At least two angular ghost points in both theta and phi
are required and checked at initialization; this is SphericalSurface's default.
Integration uses QLM's sine-series quadrature rule, including symmetry-plane
weights. Angular ghosts and padding are excluded from the integrals and extrema.
Reduced angular grids describe reflections of both the shape and metric about
the supplied origin. The mean radius is `integral(r dA)/integral(dA)`, distinct
from a plain angular average and QLM's areal radius. Extrema are sampled bounds;
they converge with angular resolution. `sf_area` and `qlm_area` can differ at
finite resolution because their interpolation and differentiation differ.

Standard text and HDF5 output are both supported; enable either or both:

```ccl
IOASCII::out0D_every = 128
IOASCII::out0D_vars = "
  SphericalSurface::sf_active SphericalSurface::sf_valid
  SphericalSurface::sf_info NSSurface::info_valid NSSurface::search_info
  SphericalSurface::sf_origin SphericalSurface::sf_coordinate_descriptors
  SphericalSurface::sf_shape_descriptors
  QuasiLocalMeasures::qlm_state QuasiLocalMeasures::qlm_scalars
"
IOASCII::out2D_every = 128
IOASCII::out2D_vars = "SphericalSurface::sf_radius"

IOHDF5::out_every = 128
IOHDF5::out_vars = "
  SphericalSurface::sf_radius SphericalSurface::sf_origin
  SphericalSurface::sf_info NSSurface::info_valid
  SphericalSurface::sf_coordinate_descriptors
  SphericalSurface::sf_shape_descriptors SphericalSurface::sf_valid
"
```

An array index is an angular grid index, not a Cartesian coordinate. Reconstruct
the angles using `sf_origin_theta/phi` and `sf_delta_theta/phi`, and respect the
active sizes and ghost widths. Retain the parfile to identify each density
level. Always check validity before plotting a stored contour. QLM's optional
`output_vtk_every` provides a ready-made 3D surface export.

QLM's `qlm_coordspin*` uses the supplied origin as its rotation centre. Keep
the centre convention consistent when comparing spin between surfaces or
runs. QLM's horizon mass and dimensionless Kerr-spin formulas should not be
interpreted as neutron-star masses or used to normalize neutron-star spin.

## Regression

`test/tov_surface.par` is a small standard Cactus regression using TOVSolver,
an octant grid and QLM, in `absolute` mode. It checks the found shape,
auxiliary surface information, search information, validity and areal radius.
It requires the usual Carpet, interpolation and QLM thorns in the configuration.
