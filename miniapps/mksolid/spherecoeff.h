#pragma once

#include "mfem.hpp"
#include <Eigen/Core>

class ImpSphereCoeff : public mfem::Coefficient {
public:
  const mfem::Vector m_center;
  const mfem::real_t m_radius;
  const mfem::real_t m_thick;

  ImpSphereCoeff(const mfem::Vector &center, const mfem::real_t &radius,
                 const mfem::real_t &thick)
      : m_center(center), m_radius(radius), m_thick(thick) {}

  /// Evaluate the coefficient at @a ip.
  virtual mfem::real_t Eval(mfem::ElementTransformation &T,
                            const mfem::IntegrationPoint &ip);

  mfem::real_t Compute(const Eigen::Vector3d &point);
};
