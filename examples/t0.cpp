//                                MFEM Tutorial 0
//
// Compile with: make t0
//
// Sample runs:  t0 -om mesh-1d.mesh
//               t0 -d 2 -t 0 -om mesh-2d-tri.mesh
//               t0 -d 2 -t 1 -om mesh-3d-quad.mesh
//               t0 -d 3 -t 0 -om mesh-3d-tet.mesh
//               t0 -d 3 -t 1 -om mesh-3d-hex.mesh
//               t0 -d 3 -t 1 -r 2 -om mesh-3d-hex-fine.mesh
//
// Description: This tutorial explores the creation of basic cartesian meshes.
// We explore changing the dimensions and the element type. The mesh can be
// uniformly refined.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
   // 1. Parse command line options.
   const char *mesh_file = "mesh.mesh";
   int subdiv = 3;
   int dim = 2;
   int etype = 1;
   double length = 1.0;
   int refinements = 0;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-om", "--output-mesh", "Mesh file output.");
   args.AddOption(&dim, "-d", "--dim", "Dimension.");
   args.AddOption(&etype, "-t", "--element-type", "Element type: 0 - simplex, 1 - quad.");
   args.AddOption(&subdiv, "-s", "--subdivision", "Number of subdivision in each dimensions.");
   args.AddOption(&length, "-l", "--length", "Length of the mesh in each dimension.");
   args.AddOption(&refinements, "-r", "--refinements", "Post-meshing refinement steps.");
   args.ParseCheck();

   // 2. Generate the mesh according to settings
   Mesh mesh;
   if (dim == 1) {
      mesh = Mesh::MakeCartesian1D(subdiv, length);
   } else if (dim == 2) {
     // FIXME: What does the option generate_edges is doing exactly?
      Element::Type type = (etype == 0) ? Element::TRIANGLE: Element::QUADRILATERAL;
      mesh = Mesh::MakeCartesian2D(subdiv, subdiv, type, false, length, length);
   } else if (dim == 3) {
     Element::Type type = (etype == 0) ? Element::TETRAHEDRON: Element::HEXAHEDRON;
     mesh = Mesh::MakeCartesian3D(subdiv, subdiv, subdiv, type, length, length, length);
   } else {
     MFEM_ABORT("invalid dimension");
   }

   // 3. Perform uniform refinement
   for (int i = 0; i < refinements; i++) {
     mesh.UniformRefinement();
   }

   // 4. Show mesh characteristics
   mesh.PrintInfo();

   // The attributes represent sub-domain of the mesh. For example, if there
   // are multiple materials in the domain, we can use the attribute to select
   // the material properties for that sub-domain. Here we print all attributes
   // of the domain, which is always 1 in that case.

   cout << "Mesh attributes:\n";
   mesh.attributes.Print(cout);

   // The boundary attributes allow to apply boundary conditions on specific
   // group of faces (3D), line (2D) or points (1D). By default, one attribute
   // is created per boundary geometry (6 in 3D, 4 in 2D and 2 in 1D).

   cout << "Mesh boundary attributes:\n";
   mesh.bdr_attributes.Print(cout);

   mesh.Save(mesh_file);
   return 0;
}
