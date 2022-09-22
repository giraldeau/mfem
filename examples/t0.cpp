//                                MFEM Tutorial 0
//
// Compile with: make t0
//
// Sample runs: ../t0 -d 1 -o 1 -out mesh-1d-o1
//              ../t0 -d 1 -o 2 -out mesh-1d-o2
//              ../t0 -d 2 -t 0 -out mesh-2d-tri-o1
//              ../t0 -d 2 -t 1 -o 1 -out mesh-3d-quad-o1
//              ../t0 -d 2 -t 1 -o 2 -out mesh-2d-quad-o2
//              ../t0 -d 2 -t 1 -o 3 -out mesh-2d-quad-o3
//              ../t0 -d 3 -t 0 -out mesh-3d-tet-o1
//              ../t0 -d 3 -t 1 -out mesh-3d-hex-o1
//              ../t0 -d 3 -t 1 -s 10 -out mesh-3d-hex-fine-o1
//
// Description: This tutorial explores the creation of basic cartesian unit
// meshes. We observe the effect of the dimensions (1, 2, or 3d), the element
// type (simplex or quad), and uniform refinement. Then, we attach a continuous
// piece-wise polynomial finite element collection (H1) to it and enumerate the
// degrees of freedom (dof) of the domain. We use a GridFunction to evaluate a
// function at each dof and save the result for visualization. We look at the
// sparsity pattern of the resulting bilinear system.

#include <fstream>
#include <iostream>
#include <sstream>

#include "mfem.hpp"

using namespace std;
using namespace mfem;

// Function returning the field value (scalar) at the given coordinate. This
// function is evaluated at each each dof when projecting the function coefficient.
double func(const Vector &v)
{
   double val = 0.0;
   if (v.Size() == 1)
   {
      val = std::sin(2 * M_PI * v[0]);
   }
   else if (v.Size() == 2)
   {
      val = std::sin(2 * M_PI * v[0]) + std::cos(2 * M_PI * v[1]);
   }
   else
   {
      val = std::sin(2 * M_PI * v[0]) + std::cos(2 * M_PI * v[1]) +
            std::sin(2 * M_PI * v[2] + M_PI / 2);
   }
   return val;
}

int main(int argc, char *argv[])
{
   // 1. Parse command line options.
   const char *out_prefix = "t1";
   int subdiv = 3;
   int dim = 2;
   int etype = 1;
   int order = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&out_prefix, "-out", "--output", "Output prefix.");
   args.AddOption(&dim, "-d", "--dim", "Dimension.");
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree.");
   args.AddOption(&etype, "-t", "--element-type",
                  "Element type: 0 - simplex, 1 - quad.");
   args.AddOption(&subdiv, "-s", "--subdivision",
                  "Number of subdivision in each dimensions.");
   args.ParseCheck();

   // 2. Generate the mesh according to settings
   Mesh mesh;
   if (dim == 1)
   {
      mesh = Mesh::MakeCartesian1D(subdiv);
   }
   else if (dim == 2)
   {
      // FIXME: What does the option generate_edges is doing exactly?
      Element::Type type = (etype == 0) ? Element::TRIANGLE: Element::QUADRILATERAL;
      mesh = Mesh::MakeCartesian2D(subdiv, subdiv, type, false);
   }
   else if (dim == 3)
   {
      Element::Type type = (etype == 0) ? Element::TETRAHEDRON: Element::HEXAHEDRON;
      mesh = Mesh::MakeCartesian3D(subdiv, subdiv, subdiv, type);
   }
   else
   {
      MFEM_ABORT("invalid dimension");
   }

   // 3. Show mesh characteristics and save the mesh
   mesh.PrintInfo();
   {
      stringstream ss;
      ss << out_prefix << ".mesh";
      mesh.Save(ss.str().c_str());
   }

   // 4. The attributes represent sub-domain of the mesh. For example, if there
   // are multiple materials in the domain, we can use the attribute to select
   // the material properties for that sub-domain. Here we print all attributes
   // of the domain, which is always 1 in that case.
   cout << "Mesh attributes:\n";
   mesh.attributes.Print(cout);

   // 5. The boundary attributes allow to apply boundary conditions on specific
   // group of points (1D), line (2D), or faces (3D). By default, one attribute
   // is created per boundary of the geometry. (2 points in 1D, 4 lines in 2D
   // and 6 faces in 3D).
   cout << "Mesh boundary attributes:\n";
   mesh.bdr_attributes.Print(cout);

   // 6. Create finite element space for the mesh. Here we use continuous
   // piece-wise polynomial elements with of specified order.
   H1_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fec);

   // 7. Display finite element caracteristics. Here, GetVSize() returns the
   // number of all dofs, whereas GetTrueVSize() returns the number of true
   // dofs. For plain mesh, the two numbers are equal. However, for
   // non-conforming mesh, true dofs is smaller than the total number of dofs.
   // In this situation, dofs at hanging nodes are interpolated from the true
   // dofs. These dofs are therefore not part of the linear system when solving
   // an equation.
   cout << "Number of dofs: " << fespace.GetVSize() << endl;
   cout << "Number of true dofs: " << fespace.GetTrueVSize() << endl;

   // 8. Next we get the true dofs index and print them. Enumerating these dof
   // is useful to apply boundary conditions.
   Array<int> bdr_tdof;
   fespace.GetBoundaryTrueDofs(bdr_tdof);
   cout << "Number of boundary true dofs: " << bdr_tdof.Size() << endl;
   bdr_tdof.Print();

   // 9. Create a grid function, project a function on this space and and save
   // it to file. The function will be called for each dof coordinate. In
   // practice, ProjectCoefficient() loops over each dofs of each elements
   // independently. Therefore, the function will be evaluated multiple times
   // with the same position for elements with a common nodes.
   GridFunction gf(&fespace);
   FunctionCoefficient coef(func);
   gf.ProjectCoefficient(coef);
   {
      stringstream ss;
      ss << out_prefix << ".gf";
      gf.Save(ss.str().c_str());
   }

   // 10. Now we generate the sparsity pattern of the linear system
   // corresponding to our finite element space. Here we use a MassIntegrator
   // term phi_i*phi_j to enumerate the dofs of the system. We do not
   // actually compute the matrix entries, we only identify where entries are
   // non-zero.
   BilinearForm system(&fespace);
   system.UsePrecomputedSparsity();
   system.AddDomainIntegrator(new MassIntegrator());
   system.AllocateMatrix();

   // 11. Access the core matrix and set all entries to 1. The assignment
   // operator=(double) implicitely loops over all non-zero entries of the
   // sparse matrix.
   SparseMatrix mat = system.SpMat();
   mat = 1.0;

   // 12. Now we write this matrix in a file. You can load this file in matlab or
   // octave to see the sparsity pattern. For instance:
   //
   //   > spy(spconvert(load("t1.mat")))
   //
   mat.PrintInfo(cout);
   {
      stringstream ss;
      ss << out_prefix << ".mat";
      std::ofstream out(ss.str());
      mat.PrintMatlab(out);
   }

   return 0;
}
