# NSSurface

Find a neutron star's `HydroBase::rho = rho_surface` contour and publish it to
SphericalSurface. QuasiLocalMeasures consumes the shape and calculates proper
area, areal radius, circumferences and spin. There is no dependency on
VolumeIntegrals or ReduceLocTracker.

The origin is read from `SphericalSurface::sf_origin_*` before each search.
At startup it is initialized from SphericalSurface's `origin_x/y/z` parameters.
An existing centre provider can update it before `NSSurface_Find`; the finder
is ordered after `VI_GRMHD_VolumeIntegralGroup` and before `qlm_calculate`.
A fixed origin is sufficient for a centred isolated star. The finder does not
move the origin or calculate a mass-weighted centre of mass.

## Minimal setup

Add `NSTools/NSSurface` to the thornlist, and activate `NSSurface`,
`SphericalSurface`, `ADMBase`, `HydroBase`, `CarpetInterp` and `AEILocalInterp` with the
usual grid and evolution thorns. For one surface:

```ccl
SphericalSurface::nsurfaces = 1
SphericalSurface::name[0] = "stellar surface"
SphericalSurface::origin_x[0] = 0
SphericalSurface::origin_y[0] = 0
SphericalSurface::origin_z[0] = 0

# Full angular grid, including two ghost points at each end.
SphericalSurface::maxntheta = 35
SphericalSurface::maxnphi = 68
SphericalSurface::ntheta[0] = 35
SphericalSurface::nphi[0] = 68
SphericalSurface::set_spherical[0] = no
SphericalSurface::set_elliptic[0] = no
SphericalSurface::auto_res[0] = no

NSSurface::find_surface[0] = yes
NSSurface::rho_surface[0] = 1e-8
NSSurface::max_radius[0] = 20
NSSurface::compute_every = 128

QuasiLocalMeasures::num_surfaces = 1
QuasiLocalMeasures::surface_name[0] = "stellar surface"
QuasiLocalMeasures::spatial_order = 4
QuasiLocalMeasures::interpolator_options = "order=4"
```

For the QLM diagnostics in this example, also activate `QuasiLocalMeasures`
and `TGRtensor` and include their dependencies in the thornlist.

The example density and search radius are choices to adapt to the star. The
cutoff must exceed the atmosphere density. The search ball must fit inside
the interpolation domain, with stencil clearance. Use the same synchronized
cadence for the centre update, HydroBase/ADMBase conversion, surface finder and QLM.
For BNS, configure a slot and origin for each star. Up to 42 slots are available,
shared with other SphericalSurface users; several contours can share one origin.

## Symmetry

Cartesian grid symmetries are handled by Carpet's registered symmetry
interpolation, including reflection and rotational symmetry. The origin and
search ball can therefore extend across a symmetry boundary.

NSSurface also accepts SphericalSurface's reduced angular grids. However,
QuasiLocalMeasures itself requires a full angular sphere: keep
`SphericalSurface::symmetric_x/y/z = no` for slots consumed by QLM. This does
not require evolving a full Cartesian domain. Keep `auto_res = no`; explicit
angular resolution makes the surface layout consistent across restarts.

## Search and validity

The origin must be inside the selected density contour. The finder samples
outward along each angular ray, brackets the first downward crossing and
bisects it using linear density interpolation. Linear interpolation avoids
high-order overshoots near the atmosphere. `radial_samples` defaults to 64;
`radial_tolerance` defaults to 1e-4 in coordinate units. These control the
bracketing and root search, not the accuracy of the underlying evolution grid.

This representation is intended for separate, intact, star-shaped stars.
It does not describe detached ejecta, overhangs or arbitrary merger topology.
Features narrower than a radial sampling interval can be missed. Check
radial/angular resolution and cutoff sensitivity when measuring oscillations.

Only selected surface slots are modified. `sf_valid` is positive only after
a complete successful search, negative on failure, and zero on skipped
iterations. Previous radii can remain in storage on failure; consumers must
check validity. For NSSurface's slots, auxiliary fields have their own
`NSSurface::info_valid` flag:
1 means current information, 0 skipped or unused, and -1 unavailable. Check it
when reading `sf_info`; those fields are NaN on skipped or failed updates.
An auxiliary failure does not change the contour's `sf_valid` or prevent QLM
from processing it. QLM still requires valid metric data for its own calculation.
All active angular points, including angular ghosts, are
searched; padding is initialized separately. Origins are never changed after
the radii have been calculated.

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
  SphericalSurface::sf_info NSSurface::info_valid
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
cutoff. Always check validity before plotting a stored contour. QLM's optional
`output_vtk_every` provides a ready-made 3D surface export.

QLM's `qlm_coordspin*` uses the supplied origin as its rotation centre. Keep
the centre convention consistent when comparing spin between surfaces or
runs. QLM's horizon mass and dimensionless Kerr-spin formulas should not be
interpreted as neutron-star masses or used to normalize neutron-star spin.

## Regression

`test/tov_surface.par` is a small standard Cactus regression using TOVSolver,
an octant grid and QLM. It checks the found shape, auxiliary surface information,
validity and areal radius.
It requires the usual Carpet, interpolation and QLM thorns in the configuration.
