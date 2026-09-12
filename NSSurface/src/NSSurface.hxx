#ifndef NSSURFACE_HXX
#define NSSURFACE_HXX

#include "cctk.h"
#include <array>
#include <vector>
#ifdef CCTK_MPI
#include <mpi.h>
#endif

namespace NSSurface {
using Point = std::array<CCTK_REAL, 3>;
using Coordinates = std::array<std::vector<CCTK_REAL>, 3>;

#ifdef CCTK_MPI
MPI_Comm communicator(const cGH *gh);
#endif
bool all_found(const cGH *gh, int found);
bool interpolate(const cGH *gh, int table, const Coordinates &coordinates,
                 int nvars, const CCTK_INT *variables,
                 std::vector<CCTK_REAL> *values);
} // namespace NSSurface
#endif
