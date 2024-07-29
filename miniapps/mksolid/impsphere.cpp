#include "mfem.hpp"
#include <iostream>
#include <Eigen/Core>

using namespace mfem;
using namespace std;


class ImpSphereCoeff : public Coefficient {
public:
    const Vector m_center;
    const real_t m_radius;

    ImpSphereCoeff(const Vector &center, const real_t &radius):
        m_center(center),
        m_radius(radius) {
    }

    /// Evaluate the coefficient at @a ip.
    virtual real_t Eval(ElementTransformation &T,
                        const IntegrationPoint &ip) {

        real_t x[3];
        Vector transip(x, 3);
        T.Transform(ip, transip);
        real_t thick = 0.1 * transip[0] + 0.05;

        Eigen::Vector3d center(m_center[0], m_center[1], m_center[2]);
        Eigen::Vector3d point(transip[0], transip[1], transip[2]);


        // Shallow sphere with varying thickness
        Eigen::Vector3d cp = point - center;
        double norm = cp.norm();
        if (norm < 1E-6 * m_radius) {
            // FIXME: how to properly handle query points close to the center?
            // The direction of two coincident points is undefined.
            return abs(norm - m_radius);
        }
        Eigen::Vector3d n = cp.normalized();
        Eigen::Vector3d mid = center + n * (m_radius + thick * 0.5);
        double res = (point - mid).norm() - thick * 0.5;

        return res;
    }

};

int main() {
    cout << "BEGIN Implicit Sphere" << endl;

    ImpSphereCoeff sphere(Vector({1, 1, 1}), 0.75);
    Mesh mesh = Mesh::MakeCartesian3D(30, 30, 30, Element::HEXAHEDRON, 2.0, 2.0, 2.0);

    H1_FECollection fec(2, mesh.Dimension());
    FiniteElementSpace fespace(&mesh, &fec);

    GridFunction x(&fespace);
    x.ProjectCoefficient(sphere);

    {
        mfem::ParaViewDataCollection paraview_dc("Implicit", &mesh);
        paraview_dc.SetPrefixPath("ParaView");
        paraview_dc.SetLevelsOfDetail(2);
        paraview_dc.SetDataFormat(VTKFormat::BINARY);
        paraview_dc.RegisterField("function",&x);
        paraview_dc.Save();
    }

    cout << "END Implicit Sphere" << endl;
    return 0;
}
