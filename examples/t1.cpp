//                                MFEM Tutorial 0
//
// Compile with: make t1
//
// Sample runs:  t1
//
// Description: In this tutorial, we are attaching a finite element collection
// to a mesh to enumerate the degrees of freedom (dof) of the domain. We use a
// GridFunction to evaluate a function at each dof and save the result for
// visualization. We look at the sparsity pattern of the resulting linear
// system.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

double func(const Vector &v)
{
  double val = 0.0;
  if (v.Size() == 1) {
    val = std::sin(2 * M_PI * v[0]);
  } else {
    val = std::sin(2 * M_PI * v[0]) + std::cos(2 * M_PI * v[1]);
  }
  return val;
}

int main(int argc, char *argv[])
{
   // 1. Parse command line options.
   const char *mesh_file = "mesh.mesh";
   const char *sol_file = "sol.gf";
   int order = 1;
   int refinements = 0;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file input.");
   args.AddOption(&sol_file, "-os", "--output-sol", "Solution file output."); // FIXME: standardize options
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");
   args.AddOption(&refinements, "-r", "--refinements", "Mesh refinement steps.");
   args.ParseCheck();

   // 2. Load the mesh, perform refinement and show caracteristics
   Mesh mesh(mesh_file);
   for (int i = 0; i < refinements; i++) {
     mesh.UniformRefinement();
   }
   mesh.PrintInfo();

   // 3. Create finite element space for the mesh. Here we use Lagrange
   // elements with of specified order.
   H1_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fec);

   // 4. Display finite element caracteristics. Here, GetVSize() returns the
   // number of all dofs, whereas GetTrueVSize() returns the number of true
   // dofs. For plain mesh, the two numbers are equal. However, for
   // non-conforming mesh, true dofs is smaller than the the overall dofs.
   // Basically, the dofs at hanging nodes are interpolated from the true dofs.
   // These dofs are not part of the linear system when solving an equation.

   cout << "Number of dofs: " << fespace.GetVSize() << endl;
   cout << "Number of true dofs: " << fespace.GetTrueVSize() << endl;

   // 5. Next we get the true dofs index and print them.
   Array<int> bdr_tdof;
   fespace.GetBoundaryTrueDofs(bdr_tdof);
   cout << "Number of boundary true dofs: " << bdr_tdof.Size() << endl;
   bdr_tdof.Print();

   // 5. Create a grid function, project a function on this space and and save
   // it to file. The function will be called for each dof coordinate. In
   // practice, ProjectCoefficient() loops over each dofs of each elements
   // independently. Therefore, the function will be evaluated multiple times
   // with the same position for elements with a common nodes.
   GridFunction gf(&fespace);
   FunctionCoefficient coef(func);
   gf.ProjectCoefficient(coef);
   gf.Save(sol_file);

   // Here we extract values for each true dofs only. Ee set the true dofs
   // vector to 1 and values at the boundary with 0. Then, we set the grid
   // function using the modified true dofs vector. If the mesh is
   // non-conforming, SetFromTrueDofs() interpolates the dependent values at
   // hanging nodes.
   Vector u;
   gf.GetTrueDofs(u);

   u.SetSubVectorComplement(bdr_tdof, 1.0);
   u.SetSubVector(bdr_tdof, 0.0);
   gf.SetFromTrueDofs(u);
   gf.Save("u.gf");

   // 6. Now we generate the sparsity pattern of the linear system
   // corresponding to our finite element space. Here we use a MassIntegrator
   // term phi(i)*phi(j)) to enumerate the dofs of the system. We do not
   // actually compute the matrix entries, we only identify where entries are
   // non-zero.
   BilinearForm system(&fespace);
   system.UsePrecomputedSparsity();
   system.AddDomainIntegrator(new MassIntegrator());
   system.AllocateMatrix();

   // Access the core matrix and set all entries to 1. The assignment
   // operator=(double) implicitely loops over all entries.
   SparseMatrix mat = system.SpMat();
   mat = 1.0;

   // Now we write this matrix in a file. You can load this file in matlab or
   // octave to see the sparsity pattern. For instance:
   //
   //   > spy(spconvert(load("system.mat")))
   //
   mat.PrintInfo(cout);
   std::ofstream out("system.mat");
   mat.PrintMatlab(out);

   return 0;
}
