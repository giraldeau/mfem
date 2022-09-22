//                                MFEM Tutorial 3
//
// Compile with: make t3
//
// Sample runs:  t3
//
// Description: This tutorial is about non-conforming meshes and how it affects
// the system and the solution. We show the difference between conforming and
// non-conforming refinement. We solve the Laplace equation as a model problem.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   // 1. Parse command line options.
   const char *mesh_file = "mesh.mesh";
   const char *sol_file = "sol.gf";
   int order = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");
   args.ParseCheck();

   // 2. Create a square 2D mesh with 4 quad elements. In this mesh, there is
   // only one unkown value at the center. We refine 3 out of 4 elements with
   // non-conforming refinement. This creates two hanging nodes adjacent to the
   // unrefined element.
   Mesh mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL);
   Array<int> el_list({0, 1, 2});
   mesh.GeneralRefinement(el_list);
   mesh.PrintInfo();

   // 3. Define a finite element space on the mesh. Here we use H1 continuous
   //    high-order Lagrange finite elements of the given order.
   H1_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fec);

   int vdofs = fespace.GetVSize();
   int tdofs = fespace.GetTrueVSize();
   int cdofs = vdofs - tdofs;

   cout << "Number of dofs             : " << vdofs << endl;
   cout << "Number of unknowns         : " << tdofs << endl;
   cout << "Number of constrained dofs : " << cdofs << endl;

   // 5. Define the solution x as a finite element grid function in fespace.
   // Set the initial guess to zero. This also sets all boundaries to zero.
   // Then, we set the left side to 1. We use a boolean marker array to apply
   // the value to the these boundary dofs. The function MakeCartesian2D()
   // creates one attribute for each sides.
   GridFunction x(&fespace);
   x = 0.0;

   {
      Array<int> bdr_marker(mesh.bdr_attributes.Max());
      bdr_marker = 0;
      bdr_marker[1] = 1;

      ConstantCoefficient one(1.0);
      x.ProjectBdrCoefficient(one, bdr_marker);

      std::ofstream out("init.gf");
      x.Save(out);
   }

   // 6. Set up the linear form b(.) corresponding to the right-hand side. Here
   // the rhs is simply zero.
   LinearForm b(&fespace);
   b.Assemble();

   // 7. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new DiffusionIntegrator);
   a.Assemble();
   a.Finalize();

   // The call to FormLinearSystem() modifies the system matrix and setup the
   // rhs. We will inspect the bilinear form matrix before the call and notice
   // how it is modified in practice.
   {
      const SparseMatrix s = a.SpMat();
      std::printf("System size before         : %d x %d\n",
                  s.Height(), s.Width());
      std::ofstream out("a1.mat");
      s.PrintMatlab(out);
   }

   // 8. Form the linear system A X = B. This includes eliminating boundary
   //    conditions, applying AMR constraints, and other transformations.
   Array<int> boundary_dofs;
   fespace.GetBoundaryTrueDofs(boundary_dofs);
   SparseMatrix A;
   Vector B, X;
   a.FormLinearSystem(boundary_dofs, x, b, A, X, B);

   // Inspect the modified system matrix and rhs. The resulting matrix is
   // smaller, that is, rows and columns of constrained dofs are deleted. The
   // solution vector X and the rhs B are also smaller to match the system.
   //
   // We also dump the prolongation and restriction matrices. The restriction
   // matrix R has size (tdofs,vdofs) with ones on the diagonal for conforming
   // dofs. But for dofs to remove from the linear system, the column is set to
   // zero. We obtain X (the vector without constrainted dofs) with the product
   // R.x.
   //
   // The prolongation matrix P works in the reverse direction. It has size
   // (vdofs,tdofs) and has mainly ones on the diagonal for conforming dofs.
   // The rows corresponding to constrained dofs have zero on the diagonal and
   // weights representing the linear interpolation of the neighbors values.
   // The result is that we get the full solution x using the product P.X.
   {
      const SparseMatrix s = a.SpMat();
      std::printf("System size after          : %d x %d\n",
                  s.Height(), s.Width());
      std::ofstream out("a2.mat");
      s.PrintMatlab(out);

      const SparseMatrix *R = fespace.GetConformingRestriction();
      std::printf("Restriction matrix size    : %d x %d\n",
                  R->Height(), R->Width());
      std::ofstream rmat("rmat.mat");
      R->PrintMatlab(rmat);

      const SparseMatrix *P = fespace.GetConformingProlongation();
      std::printf("Prolongation matrix size   : %d x %d\n",
                  P->Height(), P->Width());
      std::ofstream pmat("pmat.mat");
      P->PrintMatlab(pmat);
   }

   // 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
   cout << "Solving...\n";
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 200, 1e-12, 0.0);

   // 10. Recover the solution x as a grid function and save to file. The output
   //     can be viewed using GLVis as follows: "glvis -m mesh.mesh -g sol.gf"
   a.RecoverFEMSolution(X, b, x);
   x.Save(sol_file);
   mesh.Save(mesh_file);

   return 0;
}
