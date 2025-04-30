//                                MFEM Example 0
//
// Compile with: make ex0
//
// Sample runs:  ex0
//               ex0 -m ../data/fichera.mesh
//               ex0 -m ../data/square-disc.mesh -o 2
//
// Description: This example code demonstrates the most basic usage of MFEM to
//              define a simple finite element discretization of the Laplace
//              problem -Delta u = 1 with zero Dirichlet boundary conditions.
//              General 2D/3D mesh files and finite element polynomial degrees
//              can be specified by command line options.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   // 1. Parse command line options.
   string mesh_file = "../data/star.mesh";
   int order = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");
   args.ParseCheck();

   // 2. Read the mesh from the given mesh file, and refine once uniformly.
   Mesh mesh(mesh_file);
   mesh.UniformRefinement();
   mesh.UniformRefinement();
   //mesh = Mesh::MakeSimplicial(mesh);

   // 3. Define a finite element space on the mesh. Here we use H1 continuous
   //    high-order Lagrange finite elements of the given order.
   H1_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fec);
   cout << "Number of unknowns: " << fespace.GetTrueVSize() << endl;

   // 4. Extract the list of all the boundary DOFs. These will be marked as
   //    Dirichlet in order to enforce zero boundary conditions.
   Array<int> boundary_dofs;
   fespace.GetBoundaryTrueDofs(boundary_dofs);

   // 5. Define the solution x as a finite element grid function in fespace. Set
   //    the initial guess to zero, which also sets the boundary conditions.
   GridFunction x(&fespace);
   x = 0.0;

   // 6. Set up the linear form b(.) corresponding to the right-hand side.
   ConstantCoefficient one(0.0);
   LinearForm b(&fespace);
   b.AddDomainIntegrator(new DomainLFIntegrator(one));
   b.Assemble();

   // 7. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new DiffusionIntegrator);
   a.Assemble();
   a.Finalize();
   {
      std::ofstream outb("b.mat");
      b.Print(outb);

      std::ofstream outa("a.mat");
      a.PrintMatlab(outa);
   }

   // 8. Form the linear system A X = B. This includes eliminating boundary
   //    conditions, applying AMR constraints, and other transformations.
   SparseMatrix A;
   Vector B, X;
   a.FormLinearSystem(boundary_dofs, x, b, A, X, B);

   {
      std::ofstream outb("B.mat");
      b.Print(outb);

      std::ofstream outa("A.mat");
      A.PrintMatlab(outa);
   }

   // Lagrangian Multiplicator Bloc
   Array<int> block_offsets(3); // number of blocs + 1
   block_offsets[0] = 0;
   block_offsets[1] = fespace.GetVSize(); // system unknowns
   block_offsets[2] = 2; // number of constraints  DOF
   block_offsets.PartialSum();

   BlockVector x2(block_offsets);
   BlockVector rhs(block_offsets);

   // Lagrange Multipliers are new unkowns in the system
   Vector LM(2);
   LM = 0.0;
   x2.GetBlock(0) = X;
   x2.GetBlock(1) = LM;

   // Define the constraints in the RHS
   Vector fix(2);
   fix(0) = 0.2;
   fix(1) = -0.2;
   rhs.GetBlock(0) = B;
   rhs.GetBlock(1) = fix;

   // Define the constrained DOFs bloc.
   // Hardcoded DOFs 11 and 13
   BlockMatrix blockMatrix(block_offsets);
   SparseMatrix G(2, a.Width());
   G.Add(0, 11, 1.0);
   G.Add(1, 13, 1.0);
   G.Finalize();
   std::unique_ptr<SparseMatrix> Gt(Transpose(G));

   blockMatrix.SetBlock(0,0, &A);
   blockMatrix.SetBlock(0,1, Gt.get());
   blockMatrix.SetBlock(1, 0, &G);

   // The UMFPackSolver direct solver does not support BlocMatrix. We have to convert
   // the bloc matrix to a monolitic sparse matrix. A more efficient way would be to
   // create the monolitic sparse matrix to begin with.
   std::unique_ptr<SparseMatrix> monolithic(blockMatrix.CreateMonolithic());
   {
      std::ofstream out("system.m");
      monolithic->PrintMatlab(out);
   }

   // 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
#ifdef MFEM_USE_SUITESPARSE
   UMFPackSolver solver(*monolithic);
   solver.Mult(rhs, x2);
#else
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 200, 1e-12, 0.0);
#endif

   // 10. Recover the solution x as a grid function and save to file. The output
   //     can be viewed using GLVis as follows: "glvis -m mesh.mesh -g sol.gf"
   GridFunction u;
   u.MakeRef(&fespace, x2.GetBlock(0), 0);
   x2.GetBlock(1);

   u.Save("sol.gf");
   mesh.Save("mesh.mesh");

   {
      // ... but I prefer to use ParaView
      ParaViewDataCollection dc("ex0LM");
      dc.SetPrefixPath("ParaView");
      dc.SetMesh(&mesh);
      dc.RegisterField("u", &u);
      dc.Save();
   }

   return 0;
}
