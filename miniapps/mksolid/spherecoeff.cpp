#include "spherecoeff.h"

using namespace mfem;
using namespace std;

real_t ImpSphereCoeff::Eval(ElementTransformation &T,
                            const IntegrationPoint &ip) {

  int dim = T.GetSpaceDim();

  double res = 0.0;

  if (dim == 2) {
    real_t x[2];
    Vector transip(x, 2);
    T.Transform(ip, transip);
    MFEM_ASSERT(transip.Size() == 2, "Dim 2 expected")

    Eigen::Vector2d center(m_center[0], m_center[1]);
    Eigen::Vector2d point(transip[0], transip[1]);
    Eigen::Vector2d cp = point - center;
    double norm = cp.norm();

    if (norm < 1E-6 * m_radius) {
      res = abs(norm - m_radius);
    } else {
      Eigen::Vector2d n = cp.normalized();
      Eigen::Vector2d mid = center + n * (m_radius + m_thick * 0.5);
      res = (point - mid).norm() - m_thick * 0.5;
    }
  } else if (dim == 3) {
    real_t x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    MFEM_ASSERT(transip.Size() == 3, "Dim 3 expected")

    Eigen::Vector3d point(transip[0], transip[1], transip[2]);
    res = Compute(point);
  }
  return res;
}

real_t ImpSphereCoeff::Compute(const Eigen::Vector3d &point) {
  real_t res = 0.0;
  Eigen::Vector3d center(m_center[0], m_center[1], m_center[2]);
  // Shallow sphere with varying thickness
  Eigen::Vector3d cp = point - center;
  double norm = cp.norm();
  if (norm < 1E-6 * m_radius) {
    // FIXME: how to properly handle query points close to the center?
    // The direction of two coincident points is undefined.
    res = abs(norm - m_radius);
  } else {
    Eigen::Vector3d n = cp.normalized();
    Eigen::Vector3d mid = center + n * (m_radius + m_thick * 0.5);
    res = (point - mid).norm() - m_thick * 0.5;
  }
  return res;
}
