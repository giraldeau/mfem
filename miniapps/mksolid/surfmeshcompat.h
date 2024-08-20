#pragma once

#include <vector>

// pmp-lib is compiled with cpp-11 and mfem is not yet compatible with this
// standard. This class encapsulate de functionality

namespace pmp {
class SurfaceMesh;
}

class SurfMeshCompat {
public:
  //! Type of curvature to be computed
  //! \ingroup algorithms
  enum class Curvature {
    min,    //!< minimum curvature
    max,    //!< maximum curvature
    mean,   //!< mean curvature
    gauss,  //!< Gauss curvature
    max_abs //!< maximum absolute curvature
  };

  SurfMeshCompat();
  ~SurfMeshCompat();
  void reserve(int nv, int nedge, int ne);
  void add_vertex(double x, double y, double z);
  void add_face(int *index, int n);
  void get_curvature(std::vector<double> &data, Curvature c, int smoothing_step,
                     bool use_tensor, bool use_two_ring);

private:
  pmp::SurfaceMesh *p_imp;
};
