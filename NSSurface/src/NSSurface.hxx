#ifndef NSSURFACE_HXX
#define NSSURFACE_HXX

#include "cctk.h"
#include "NSTools.hxx"
#include <array>
#include <vector>

namespace NSSurface {
using namespace NSTools;
bool surface_owned(int s);
std::vector<CCTK_REAL> theta_weights(int ni, int gt, CCTK_REAL origin_theta,
                                     CCTK_REAL dt, CCTK_REAL dp,
                                     const std::array<bool, 3> &reflected);
} // namespace NSSurface
#endif
