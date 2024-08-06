#include "mfem.hpp"
#include <iostream>
#include <Eigen/Core>

#include "../meshing/mesh-fitting.hpp"

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
        MFEM_ASSERT(transip.Size() == 3, "Dim 3 expected")

        real_t thick = m_radius * 0.1;
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

void SetMaterial(Mesh &mesh, GridFunction &mat, const GridFunction &surf_fit_gf0) {
    // Set material gridfunction
    for (int i = 0; i < mesh.GetNE(); i++)
    {
        mat(i) = material_id(i, surf_fit_gf0);
        mesh.SetAttribute(i, mat(i) + 1);
    }
}

void SelectDofsFitting(ParMesh &pmesh,
                       ParGridFunction &mat,
                       ParGridFunction &surf_fit_gf0,
                       ParGridFunction &surf_fit_mat_gf,
                       Array<bool> &surf_fit_marker)
{
    mat.ExchangeFaceNbrData();
    const Vector &FaceNbrData = mat.FaceNbrData();
    for (int j = 0; j < surf_fit_marker.Size(); j++)
    {
        surf_fit_marker[j] = false;
    }
    surf_fit_mat_gf = 0.0;

    Array<int> dof_list;
    Array<int> dofs;
    for (int i = 0; i < pmesh.GetNumFaces(); i++)
    {
        auto tr = pmesh.GetInteriorFaceTransformations(i);
        if (tr != NULL)
        {
            int mat1 = mat(tr->Elem1No);
            int mat2 = mat(tr->Elem2No);
            if (mat1 != mat2)
            {
                surf_fit_gf0.ParFESpace()->GetFaceDofs(i, dofs);
                dof_list.Append(dofs);
            }
        }
    }
    for (int i = 0; i < pmesh.GetNSharedFaces(); i++)
    {
        auto tr = pmesh.GetSharedFaceTransformations(i);
        if (tr != NULL)
        {
            int faceno = pmesh.GetSharedFace(i);
            int mat1 = mat(tr->Elem1No);
            int mat2 = FaceNbrData(tr->Elem2No - pmesh.GetNE());
            if (mat1 != mat2)
            {
                surf_fit_gf0.ParFESpace()->GetFaceDofs(faceno, dofs);
                dof_list.Append(dofs);
            }
        }
    }
    for (int i = 0; i < dof_list.Size(); i++)
    {
        surf_fit_marker[dof_list[i]] = true;
        surf_fit_mat_gf(dof_list[i]) = 1.0;
    }
}

double ComputeMinDetJ(ParMesh &pmesh,
                      ParFiniteElementSpace *pfespace,
                      IntegrationRules &irules,
                      int quad_order)
{
    real_t min_detJ = infinity();
    const int NE = pmesh.GetNE();
    for (int i = 0; i < NE; i++)
    {
        const IntegrationRule &ir =
            irules.Get(pfespace->GetFE(i)->GetGeomType(), quad_order);
        ElementTransformation *transf = pmesh.GetElementTransformation(i);
        for (int j = 0; j < ir.GetNPoints(); j++)
        {
            transf->SetIntPoint(&ir.IntPoint(j));
            min_detJ = min(min_detJ, transf->Jacobian().Det());
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &min_detJ, 1,
                  MPITypeMap<real_t>::mpi_type, MPI_MIN, MPI_COMM_WORLD);

    return min_detJ;
}

int main(int argc, char **argv) {
    cout << "BEGIN Implicit Sphere" << endl;

    Mpi::Init(argc, argv);
    int myid = Mpi::WorldRank();
    Hypre::Init();

    int mesh_poly_deg = 1;
    int quad_order = 8;
    real_t surface_fit_const = 100.0;

    ImpSphereCoeff ls_coeff(Vector({0.5, 0.5, 0.5}), 0.45);
    std::unique_ptr<ParMesh> pmesh;
    {
        std::unique_ptr<Mesh> mesh(new Mesh(Mesh::MakeCartesian3D(20, 20, 20, Element::HEXAHEDRON, 1.0, 1.0, 1.0)));
        mesh->EnsureNCMesh();
        pmesh.reset(new ParMesh(MPI_COMM_WORLD, *mesh));
    }

    if (true) {
        // Example of refinement
        cout << "mesh size before: " << pmesh->GetNE() << endl;
        Array<int> el;
        el.Append(0);
        pmesh->GeneralRefinement(el);
        cout << "mesh size after: " << pmesh->GetNE() << endl;
    }

    // Finite element space on the level-set mesh
    H1_FECollection fec(mesh_poly_deg, pmesh->Dimension());
    ParFiniteElementSpace pfespace(pmesh.get(), &fec, pmesh->Dimension());

    // === TMOP ===
    ParGridFunction x(&pfespace);
    pmesh->SetNodalGridFunction(&x);
    x.SetTrueVector();

    // Initial nodes position
    ParGridFunction x0(&pfespace);
    x0 = x;

    // TMOP Metric for 3D mesh
    TMOP_Metric_303 metric;
    TargetConstructor::TargetType target_t = TargetConstructor::IDEAL_SHAPE_UNIT_SIZE;
    TargetConstructor target_c(target_t, MPI_COMM_WORLD);
    target_c.SetNodes(x0);

    // Must be allocated, owned by the NonlinearForm
    TMOP_Integrator *tmop_integ = new TMOP_Integrator(&metric, &target_c);
    tmop_integ->SetIntegrationRules(IntRulesLo, quad_order);

    // Q: what is this doing?
    pmesh->ExchangeFaceNbrData();

    // Surface fitting
    L2_FECollection mat_coll(0, pmesh->Dimension());
    H1_FECollection surf_fit_fec(mesh_poly_deg, pmesh->Dimension());

    ParFiniteElementSpace surf_fit_fes(pmesh.get(), &surf_fit_fec);
    ParFiniteElementSpace mat_fes(pmesh.get(), &mat_coll);

    ParGridFunction mat(&mat_fes);
    ParGridFunction surf_fit_mat_gf(&surf_fit_fes);
    ParGridFunction surf_fit_mat_gf_interface(&surf_fit_fes);
    ParGridFunction surf_fit_gf0(&surf_fit_fes);
    Array<bool> surf_fit_marker(surf_fit_gf0.Size());
    ConstantCoefficient surf_fit_coeff(surface_fit_const);

    // Evaluate the level-set
    surf_fit_gf0.ProjectCoefficient(ls_coeff);
    SetMaterial(*pmesh, mat, surf_fit_gf0);

    GridFunctionCoefficient coeff_mat(&mat);
    surf_fit_mat_gf.ProjectDiscCoefficient(coeff_mat, GridFunction::ARITHMETIC);
    // FIXME: I don't understand why we do this shananigan, seems to yield identity
    surf_fit_mat_gf.SetTrueVector(); // set internal dof vector using the restriction matrix
    surf_fit_mat_gf.SetFromTrueVector(); // apply the conforming prolongation matrix

    SelectDofsFitting(*pmesh.get(), mat, surf_fit_gf0,
                      surf_fit_mat_gf_interface, surf_fit_marker);

    AdvectorCG adapt_surface;

    cout << "Setup surface fitting... " << flush;
    tmop_integ->EnableSurfaceFitting(surf_fit_gf0, surf_fit_marker,
                                     surf_fit_coeff, adapt_surface);

    cout << "Done! " << endl;

    pmesh->SetAttributes();

    // NonlinearForm
    ParNonlinearForm a(&pfespace);
    ConstantCoefficient *metric_coeff1 = NULL;
    a.AddDomainIntegrator(tmop_integ);

    double min_detJ = ComputeMinDetJ(*pmesh.get(), &pfespace, IntRulesLo, quad_order);
    if (myid == 0) {
        cout << "Minimum det(J) of the original mesh is " << min_detJ << endl;
    }

    {
        ParaViewDataCollection paraview_dc("LevelSet", pmesh.get());
        paraview_dc.SetPrefixPath("ParaView");
        paraview_dc.SetLevelsOfDetail(2);
        paraview_dc.SetDataFormat(VTKFormat::BINARY);
        paraview_dc.RegisterField("levelset",&surf_fit_gf0);
        paraview_dc.RegisterField("material",&surf_fit_mat_gf);
        paraview_dc.RegisterField("material_interface",&surf_fit_mat_gf_interface);
        paraview_dc.Save();
    }

    cout << "END Implicit Sphere" << endl;
    return 0;
}
