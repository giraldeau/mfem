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
// Description:  This example code solves a simple linear elasticity problem
//               describing a multi-material cantilever beam.
//
//               Specifically, we approximate the weak form of -div(sigma(u))=0
//               where sigma(u)=lambda*div(u)*I+mu*(grad*u+u*grad) is the stress
//               tensor corresponding to displacement field u, and lambda and mu
//               are the material Lame constants. The boundary conditions are
//               u=0 on the fixed part of the boundary with attribute 1, and
//               sigma(u).n=f on the remainder with f being a constant pull down
//               vector on boundary elements with attribute 2, and zero
//               otherwise. The geometry of the domain is assumed to be as
//               follows:
//
//                                 +----------+----------+
//                    boundary --->| material | material |<--- boundary
//                    attribute 1  |    1     |    2     |     attribute 2
//                    (fixed)      +----------+----------+     (pull down)
//
//               The example demonstrates the use of high-order and NURBS vector
//               finite element spaces with the linear elasticity bilinear form,
//               meshes with curved elements, and the definition of piece-wise
//               constant and vector coefficient objects. Static condensation is
//               also illustrated.
//
//               We recommend viewing Example 1 before viewing this example.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

class CircleDistance
{
public:
   CircleDistance(real_t centerX, real_t centerY, real_t radius)
      : centerX(centerX), centerY(centerY), radius(radius)
   {
   }

   // Method to compute the implicit distance to the circle
   real_t operator()(const Vector& x) const
   {
      double dx = x(0) - centerX;
      double dy = x(1) - centerY;
      double d = std::sqrt(dx * dx + dy * dy);
      return d - radius;
   }

   // Method to compute the direction vector (normalized) at a given location
   void operator()(const Vector &x, Vector &dir) const
   {
      real_t dx = x(0) - centerX;
      real_t dy = x(1) - centerY;
      real_t magnitude = std::sqrt(dx * dx + dy * dy);
      if (magnitude > 1e-9)
      {
         dir(0) = -dx / magnitude;
         dir(1) = -dy / magnitude;
      }
   }

private:
   real_t centerX;
   real_t centerY;
   real_t radius;
};

SparseMatrix* BuildContactMatrix(FiniteElementSpace& fespace, CircleDistance &dist)
{
   int dim = fespace.GetVDim();

   // Count the number of contact constraints
   int n_rows = 0;
   for (int i = 0; i < fespace.GetNBE(); ++i)
   {
      Array<int> dofs;
      fespace.GetBdrElementDofs(i, dofs);
      std::cout << i << " ";
      dofs.Print();
      std::cout << std::endl;
   }


   SparseMatrix * mout = new SparseMatrix(n_rows, fespace.GetTrueVSize());
   mout->Finalize();

   return mout;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int order = 1;
   bool static_cond = false;
   bool visualization = 1;
   real_t lambda = 1.0;
   real_t mu = 1.0;

   OptionsParser args(argc, argv);
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&lambda, "-l", "--lambda","Lambda parameter");
   args.AddOption(&mu, "-m", "--mu", "Mu parameter");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
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
   Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(10, 10, Element::TRIANGLE, false, 10.0, 5.0));
   int dim = mesh->Dimension();

   FiniteElementCollection *fec = new H1_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec, dim);
   FiniteElementSpace *fespace_ls = new FiniteElementSpace(mesh, fec);
   mesh->EnsureNodes(); // make sure the node array exists
   mesh->SetNodalFESpace(fespace); // required to move the nodes
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl << "Assembling: " << std::endl;

   // BC
   Array<int> ess_bdr(mesh->bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[0] = 1;
   std::cout << "Essential Boundaries: ";
   ess_bdr.Print();
   std::cout << std::endl;

   Array<int> ess_tdof_list;
   fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   std::cout << "Essential TDof: " << std::endl;
   ess_tdof_list.Print();

   // There is no forces on the RHS
   LinearForm *b = new LinearForm(fespace);
   cout << "r.h.s. ... " << std::endl;
   b->Assemble();

   GridFunction x(fespace);
   x = 0.0;

   ConstantCoefficient lambda_func(lambda);
   ConstantCoefficient mu_func(mu);

   BilinearForm *a = new BilinearForm(fespace);
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func,mu_func));

   cout << "matrix ... " << std::endl;
   a->Assemble();

   SparseMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   cout << "done." << endl;

   cout << "Size of linear system: " << A.Height() << endl;


   ParaViewDataCollection dc("Contact01");
   dc.SetPrefixPath("ParaView");
   dc.SetMesh(mesh);

   /*
    * Contact
    */
   CircleDistance circle(5, 8, 4);
   FunctionCoefficient ls_coefficient(circle);
   VectorFunctionCoefficient lsv_coefficient(dim, circle);

   GridFunction x_ls(fespace_ls);
   x_ls.ProjectCoefficient(ls_coefficient);

   GridFunction x_dir(fespace);
   x_dir.ProjectCoefficient(lsv_coefficient);

   dc.RegisterField("ls", &x_ls);
   dc.RegisterField("dir", &x_dir);
   dc.Save();

   //SparseMatrix *C = BuildContactMatrix(*fespace);
   exit(0);

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
