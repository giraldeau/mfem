#include "mfem.hpp"
#include <iostream>
#include <Eigen/Core>

#include "../meshing/mesh-fitting.hpp"
#include "spherecoeff.h"

using namespace mfem;

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

    int mesh_poly_deg = 3;
    int quad_order = 8;
    real_t surface_fit_const = 100.0;

    ImpSphereCoeff ls_coeff(Vector({0.5, 0.5, 0.5}), 0.25, 0.1);
    std::unique_ptr<ParMesh> pmesh;
    {
        //std::unique_ptr<Mesh> mesh(new Mesh(Mesh::MakeCartesian3D(20, 20, 20, Element::HEXAHEDRON, 1.0, 1.0, 1.0)));
        std::unique_ptr<Mesh> mesh(new Mesh(Mesh::MakeCartesian2D(20, 20, Element::TRIANGLE, 1.0, 1.0)));
        mesh->EnsureNCMesh();
        pmesh.reset(new ParMesh(MPI_COMM_WORLD, *mesh));
    }

    const int dim = pmesh->Dimension();

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
    std::unique_ptr<TMOP_QualityMetric> metric;
    if (dim == 2) {
        metric.reset(new TMOP_Metric_058);
    } else if (dim == 3) {
        metric.reset(new TMOP_Metric_303);
    }

    TargetConstructor::TargetType target_t = TargetConstructor::IDEAL_SHAPE_UNIT_SIZE;
    TargetConstructor target_c(target_t, MPI_COMM_WORLD);
    target_c.SetNodes(x0);

    IntegrationRules *irules = &IntRulesLo;

    // Must be allocated, owned by the NonlinearForm
    TMOP_Integrator *tmop_integ = new TMOP_Integrator(metric.get(), &target_c);
    tmop_integ->SetIntegrationRules(*irules, quad_order);

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

    {
        ParaViewDataCollection paraview_dc("Debug", pmesh.get());
        paraview_dc.SetPrefixPath("ParaView");
        paraview_dc.SetLevelsOfDetail(2);
        paraview_dc.SetDataFormat(VTKFormat::BINARY);
        paraview_dc.RegisterField("levelset",&surf_fit_gf0);
        paraview_dc.RegisterField("material",&surf_fit_mat_gf);
        paraview_dc.RegisterField("material_interface",&surf_fit_mat_gf_interface);
        paraview_dc.Save();
    }

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

    double min_detJ = ComputeMinDetJ(*pmesh.get(), &pfespace, *irules, quad_order);
    if (myid == 0) {
        cout << "Minimum det(J) of the original mesh is " << min_detJ << endl;
    }

    const real_t init_energy = a.GetParGridFunctionEnergy(x);
    real_t init_metric_energy = init_energy;

    // Note: the surface fit coefficient is part of the function energy. By
    // setting it to zero, we effectively exclude it. But why do we take it
    // into account in the first place?

    if (surface_fit_const > 0.0)
    {
        surf_fit_coeff.constant   = 0.0;
        init_metric_energy = a.GetParGridFunctionEnergy(x);
        surf_fit_coeff.constant  = surface_fit_const;
    }

    // Fix boundaries
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;
    a.SetEssentialBC(ess_bdr);

    // Use the MINRES + Jacobi
    MINRESSolver *minres = new MINRESSolver(MPI_COMM_WORLD);
    minres->SetMaxIter(100);
    minres->SetRelTol(1e-12);
    minres->SetAbsTol(0.0);
    minres->SetPrintLevel(1);

    // auto hs = new HypreSmoother;
    // hs->SetType(HypreSmoother::l1Jacobi);
    // hs->SetPositiveDiagonal(true);
    // minres->SetPreconditioner(*hs);

    // Assume homogenous mesh... but later we supply all integration rules for mixed mesh?
    Geometry::Type geom_type = pfespace.GetFE(0)->GetGeomType();
    const IntegrationRule &ir = irules->Get(geom_type, quad_order);
    TMOPNewtonSolver solver(pfespace.GetComm(), ir, 0); /* 0 = Newton */
    solver.SetIntegrationRules(*irules, quad_order);
    solver.SetPreconditioner(*minres);
    solver.SetMaxIter(200);
    solver.SetRelTol(1e-5);
    solver.SetAbsTol(0.0);
    solver.SetMinimumDeterminantThreshold(0.001*min_detJ); // 1.5e-5
    solver.SetPrintLevel(1);
    solver.SetOperator(a);

    // Solve it!
    Vector b(0);
    solver.Mult(b, x.GetTrueVector());
    x.SetFromTrueVector();

    // Displacement
    x0 -= x;

    // Compute the final energy of the functional.
    const real_t fin_energy = a.GetParGridFunctionEnergy(x);
    real_t fin_metric_energy = fin_energy;
    if (surface_fit_const > 0.0)
    {
        surf_fit_coeff.constant  = 0.0;
        fin_metric_energy  = a.GetParGridFunctionEnergy(x);
        surf_fit_coeff.constant  = surface_fit_const;
    }

    if (myid == 0)
    {
        std::cout << std::scientific << std::setprecision(4);
        cout << "Initial strain energy: " << init_energy
             << " = metrics: " << init_metric_energy
             << " + extra terms: " << init_energy - init_metric_energy << endl;
        cout << "  Final strain energy: " << fin_energy
             << " = metrics: " << fin_metric_energy
             << " + extra terms: " << fin_energy - fin_metric_energy << endl;
        cout << "The strain energy decreased by: "
             << (init_energy - fin_energy) * 100.0 / init_energy << " %." << endl;

        // FIXME: this is blocking
        if (false) {
            real_t err_avg, err_max;
            tmop_integ->GetSurfaceFittingErrors(x, err_avg, err_max);
            if (myid == 0)
            {
                std::cout << "Avg fitting error: " << err_avg << std::endl
                          << "Max fitting error: " << err_max << std::endl;
            }
        }
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
