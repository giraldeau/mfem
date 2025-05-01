//                                MFEM Example 2
//
// Compile with: make ex2
//
// Sample runs:  ex2 -m ../data/beam-tri.mesh
//               ex2 -m ../data/beam-quad.mesh
//               ex2 -m ../data/beam-tet.mesh
//               ex2 -m ../data/beam-hex.mesh
//               ex2 -m ../data/beam-wedge.mesh
//               ex2 -m ../data/beam-quad.mesh -o 3 -sc
//               ex2 -m ../data/beam-quad-nurbs.mesh
//               ex2 -m ../data/beam-hex-nurbs.mesh
//
// Test the linear elasticity with a surface mesh in the XY plane. The membrane is stretched in the +X direction.
//                   3
//    FIXED 3 |-----------| 2 Displacement >>>
//            |           |
//            |           |
//            |-----------|
//                  0

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

static int cnt = 0;

void InitDisplacement(const Vector &x, Vector &u)
{
   std::cout << "InitDisplacement: " << (cnt++) << " " << x(0) << " " << x(1) << std::endl;
   u = 0.0;
   u(0) = 0.1;
   // twish sheet
   if (u.Size() == 3) {
      u(2) = 0.2 * x(1) - 0.1;
   }
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../data/inline-quad.mesh";
   int order = 1;
   bool static_cond = false;
   bool visualization = 1;
   int refinements = 0;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&refinements, "-rf", "--refinements", "Number of refinements.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral or hexahedral elements with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();
   int spaceDim = mesh->SpaceDimension();

   // convert to 3d space
   if (spaceDim == 2) {
      mesh->SetCurvature(1, false, 3);
      ofstream mesh_ofs("initial3d.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);
      exit(0);
   }

   mesh->EnsureNodes();


   // 4. Refine the mesh to increase the resolution.
   for (int l = 0; l < refinements; l++)
   {
      mesh->UniformRefinement();
   }

   std::cout << "dim: " << dim << std::endl;
   std::cout << "spaceDim: " << spaceDim << std::endl;
   std::cout << "NV: " << mesh->GetNV() << std::endl;
   std::cout << "NE: " << mesh->GetNE() << std::endl;

   // FIXME: If we have a 2D surface in 3D space and want to use the ElasticityIntegrator, I think we need
   // to use FiniteElementCollection of dim (integration on the surface), but a FiniteElementSpace of vdim=3,
   // because we have 3 degrees of freedom (x, y, z) per node.
   //
   // As such, we will need a new integrator for surface meshes in 3D implementing BST for instance.
   //
   // The other solution is to create a volume mesh by extrusion and do the calculation in 3D.

   // 5. Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
   //    associated with the mesh nodes.
   FiniteElementCollection *fec = new H1_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec, dim);
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl << "Assembling: " << flush;


#if FALSE
   for (int i = 0; i < mesh->bdr_attributes.Size(); i++) {
      Array<int> ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
      ess_bdr = 0;
      ess_bdr[i] = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

      std::cout << "Boundary " << i << std::endl;
      for (int j = 0; j < ess_tdof_list.Size(); ++j) {
         std::cout << "   TDOF " << i << " " << ess_tdof_list[j] << std::endl;
         
      }

   }
   exit(0);
#endif


   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking only
   //    boundary attribute 1 from the mesh as essential and converting it to a
   //    list of true dofs.
   Array<int> ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[1] = 1;
   //ess_bdr[3] = 1;
   fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system. In this case, there is no force.

   LinearForm *b = new LinearForm(fespace);
   cout << "r.h.s. ... " << flush;
   b->Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0.0;

   // Impose the displacement on the boundary 1
   Vector disp(spaceDim);
   disp = 0.0;
   disp(0) = 0.1;
   disp(1) = 0.1;
   disp(2) = 0.1;
   VectorConstantCoefficient init_xc(disp);
   VectorFunctionCoefficient init_xf(spaceDim, InitDisplacement);

   Array<int> dc_bdr(mesh->bdr_attributes.Max());
   dc_bdr = 0;
   dc_bdr[1] = 1;
   x.ProjectBdrCoefficient(init_xc, dc_bdr);

   // FIXME: mesh nodes and the grid functions are incompatible
   // What does SetNodalFESpace exactly? What if we want to keep a linear mesh?
   mesh->SetNodalFESpace(fespace);
   // mesh->SetCurvature(order);

   // Save the initial solution
   {
      GridFunction *nodes = mesh->GetNodes();
      *nodes += x;
      std::cout << "Nodes: " << std::endl;
      nodes->Print();
      std::cout << "Displacements:" << std::endl;
      x.Print();

      ofstream mesh_ofs("initial.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);
      *nodes -= x;
   }

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.
   ConstantCoefficient lambda_coeff(50.0);
   ConstantCoefficient mu_coeff(1.0);
   BilinearForm *a = new BilinearForm(fespace);
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff,mu_coeff));

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   cout << "matrix ... " << flush;
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();

   SparseMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   cout << "done." << endl;

   cout << "Size of linear system: " << A.Height() << endl;

#ifndef MFEM_USE_SUITESPARSE
   // 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system Ax=b with PCG.
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 500, 1e-8, 0.0);
#else
   // 11. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(A);
   umf_solver.Mult(B, X);
#endif

   // 12. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, *b, x);

   // 14. Save the displaced mesh and the inverted solution (which gives the
   //     backward displacements to the original grid). This output can be
   //     viewed later using GLVis: "glvis -m displaced.mesh -g sol.gf".
   {
      GridFunction *nodes = mesh->GetNodes();
      *nodes += x;
      x *= -1;

      ParaViewDataCollection dc("ex2s");
      dc.SetPrefixPath("ParaView");
      dc.SetHighOrderOutput(true);
      dc.SetLevelsOfDetail(4);
      dc.SetMesh(mesh);
      dc.RegisterField("x", &x);
      dc.Save();

      ofstream mesh_ofs("displaced.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);
      ofstream sol_ofs("sol.gf");
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }


   // 16. Free the used memory.
   delete a;
   delete b;
   if (fec)
   {
      delete fespace;
      delete fec;
   }
   delete mesh;

   return 0;
}
