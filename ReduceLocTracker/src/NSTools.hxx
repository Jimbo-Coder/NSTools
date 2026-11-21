// Shared by the NSTools thorns (exported as NSTools.hxx): parameter and
// cadence lookup, reductions over processes, local maxima of a grid function
// gathered over all processes, and interpolation along straight segments.
// Thorns including it must USE FUNCTION GetMPICommWorld.
#ifndef NSTOOLS_HXX
#define NSTOOLS_HXX

#include "cctk.h"
#include "cctk_Functions.h"
#include "cctk_Parameter.h"
#include "util_Table.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>
#ifdef CCTK_MPI
#include <mpi.h>
#endif

namespace NSTools {
using Point = std::array<CCTK_REAL, 3>;
using Coordinates = std::array<std::vector<CCTK_REAL>, 3>;

// Points on a segment between two extrema, for joined/separate tests.
constexpr int segment_samples = 64;

inline bool finite(const Point &p) {
  return std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]);
}

inline CCTK_REAL distance2(const Point &a, const Point &b) {
  return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
         (a[2] - b[2]) * (a[2] - b[2]);
}

inline CCTK_REAL distance(const Point &a, const Point &b) {
  return std::sqrt(distance2(a, b));
}

// An integer or boolean parameter of another thorn, or `fallback`.
inline int parameter_int(const char *name, const char *thorn, int fallback) {
  int type = -1;
  const void *value = CCTK_ParameterGet(name, thorn, &type);
  if (value == nullptr || (type != PARAMETER_INT && type != PARAMETER_BOOLEAN))
    return fallback;
  return static_cast<int>(*static_cast<const CCTK_INT *>(value));
}

// compute_every = -2: the 0D/scalar output cadence, else QLM's, else
// IO::out_every, else 128.
inline int output_cadence(int requested) {
  if (requested != -2)
    return requested;
  const std::pair<const char *, const char *> cadences[] = {
      {"out0D_every", "CarpetIOASCII"},      {"out0D_every", "IOASCII"},
      {"outScalar_every", "CarpetIOScalar"}, {"outScalar_every", "IOBasic"},
      {"outScalar_every", "IOScalar"},       {"calculate_every", "QuasiLocalMeasures"},
      {"out_every", "IOUtil"},               {"out_every", "IO"}};
  for (const auto &[name, thorn] : cadences) {
    const int every = parameter_int(name, thorn, -1);
    if (every > 0)
      return every;
  }
  return 128;
}

// In-place reduction ("minimum", "maximum", "sum") over all processes.
inline void reduce(const cGH *gh, const char *operation, int count, CCTK_REAL *values) {
  std::vector<CCTK_REAL> result(count);
  const int handle = CCTK_ReductionArrayHandle(operation);
  if (handle < 0 || CCTK_ReduceLocArrayToArray1D(gh, -1, handle, values, result.data(),
                                                 count, CCTK_VARIABLE_REAL) < 0)
    CCTK_ERROR("NSTools: cannot reduce across processes");
  std::copy(result.begin(), result.end(), values);
}

#ifdef CCTK_MPI
inline MPI_Comm communicator(const cGH *gh) {
  return CCTK_IsFunctionAliased("GetMPICommWorld")
             ? *static_cast<const MPI_Comm *>(GetMPICommWorld(gh))
             : MPI_COMM_WORLD;
}
#endif

// True on every process when `found` is true on every process.
inline bool all_found(const cGH *gh, int found) {
#ifdef CCTK_MPI
  MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MIN, communicator(gh));
#endif
  return found != 0;
}

// True when value(p) is at least each of its 26 neighbours and above at least
// one of them, so flat regions (atmosphere) are not maxima. `dj` and `dk` are
// the index strides; p must not lie on the outermost layer of the array.
template <typename Value>
bool local_maximum(const Value &value, int p, int dj, int dk) {
  const CCTK_REAL centre = value(p);
  bool above = false;
  for (int c = -1; c <= 1; ++c)
    for (int b = -1; b <= 1; ++b)
      for (int a = -1; a <= 1; ++a) {
        const CCTK_REAL neighbour = value(p + a + b * dj + c * dk);
        if (neighbour > centre)
          return false;
        above |= neighbour < centre;
      }
  return above;
}

struct Extremum {
  CCTK_REAL value; // the value being maximised
  CCTK_REAL raw;   // the grid function there
  Point position;
};
inline bool stronger(const Extremum &a, const Extremum &b) { return a.value > b.value; }

// This process's strongest local maxima (a min-heap of at most `capacity`),
// and all processes' maxima, strongest first.
class Extrema {
public:
  static constexpr std::size_t capacity = 32;
  void clear() { heap.clear(); }
  void add(const Extremum &e) {
    if (heap.size() == capacity) {
      if (!(e.value > heap.front().value))
        return;
      std::pop_heap(heap.begin(), heap.end(), stronger);
      heap.pop_back();
    }
    heap.push_back(e);
    std::push_heap(heap.begin(), heap.end(), stronger);
  }
  std::vector<Extremum> gather(const cGH *gh) const {
    constexpr int width = 5;
    const CCTK_REAL none = -std::numeric_limits<CCTK_REAL>::max();
    std::vector<CCTK_REAL> mine(width * capacity, none);
    for (std::size_t c = 0; c < heap.size(); ++c) {
      mine[width * c] = heap[c].value;
      mine[width * c + 1] = heap[c].raw;
      for (int a = 0; a < 3; ++a)
        mine[width * c + 2 + a] = heap[c].position[a];
    }
    std::vector<CCTK_REAL> all = mine;
#ifdef CCTK_MPI
    all.resize(mine.size() * CCTK_nProcs(gh));
    MPI_Allgather(mine.data(), mine.size(), MPI_DOUBLE, all.data(), mine.size(),
                  MPI_DOUBLE, communicator(gh));
#endif
    std::vector<Extremum> extrema;
    for (std::size_t c = 0; c < all.size(); c += width)
      if (all[c] > none)
        extrema.push_back({all[c], all[c + 1], {all[c + 2], all[c + 3], all[c + 4]}});
    std::stable_sort(extrema.begin(), extrema.end(), stronger);
    return extrema;
  }

private:
  std::vector<Extremum> heap;
};

// Interpolates grid functions (linear Lagrange) at this process's points;
// every process requests its own points and Carpet routes them to the grid
// owners. True on every process when every process got finite values.
inline bool interpolate(const cGH *gh, int table, const Coordinates &coordinates,
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
    for (const auto value : values[v])
      finite &= std::isfinite(value);
  return all_found(gh, finite);
}

// Smallest of transform(f) over `samples` + 1 evenly spaced points from a to b
// (both included), on every process; false if interpolation fails.
template <typename Transform>
bool segment_minimum(const cGH *gh, int variable, const Point &a, const Point &b,
                     int samples, const Transform &transform, CCTK_REAL &minimum) {
  const int points = CCTK_MyProc(gh) == 0 ? samples + 1 : 0;
  Coordinates segment;
  for (int axis = 0; axis < 3; ++axis)
    for (int k = 0; k < points; ++k)
      segment[axis].push_back(a[axis] + (b[axis] - a[axis]) * k / samples);
  const int table = Util_TableCreateFromString("order=1");
  if (table < 0)
    CCTK_ERROR("NSTools: cannot create interpolation table");
  const CCTK_INT index = variable;
  std::vector<CCTK_REAL> values;
  const bool ok = interpolate(gh, table, segment, 1, &index, &values);
  Util_TableDestroy(table);
  minimum = std::numeric_limits<CCTK_REAL>::max();
  for (const CCTK_REAL v : values)
    minimum = std::min(minimum, CCTK_REAL(transform(v)));
#ifdef CCTK_MPI
  MPI_Allreduce(MPI_IN_PLACE, &minimum, 1, MPI_DOUBLE, MPI_MIN, communicator(gh));
#endif
  return ok;
}
} // namespace NSTools
#endif
