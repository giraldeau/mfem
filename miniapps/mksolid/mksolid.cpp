#include "mfem.hpp"
#include <iostream>

#include "tmd/TriangleMeshDistance.h"
#include <Eigen/Geometry>

/*
 * Compute the distance field from a surface mesh.
 */

using namespace mfem;
using namespace std;

typedef Eigen::AlignedBox<double, 3> AABB;

int main(int argc, char *argv[]) {

  int refinement = 2;
  std::string outname = "solid";
  std::string mesh_name;
  int simplex = 0;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_name, "-m", "--mesh", "Input mesh (surface in 3D)");
  args.AddOption(&outname, "-o", "--output", "Output file basename");

  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return 1;
  }
  args.PrintOptions(std::cout);

  Mesh surf(mesh_name);

  // The distance calculation requires triangles
  if (surf.HasGeometry(Geometry::SQUARE)) {
    surf = Mesh::MakeSimplicial(surf);
  }

  // Compute the AABB for the mesh
  AABB surf_aabb;
  ComputeMeshAABB(surf, surf_aabb);
  std::cout << "Input mesh AABB min:\n"
            << surf_aabb.min() << "\n"
            << surf_aabb.max() << "\n";

  // Build the triangle distance index
  tmd::TriangleMeshDistance tmd;
  tmd.construct(surf);

  // Create a mesh to evaluate the distance field
  // FIXME: uniform grid is not scalable. Implement AMR
  Eigen::Vector3d box = 3.0 * surf_aabb.sizes();
  Mesh mesh(Mesh::MakeCartesian3D(100, 5, 50, Element::QUADRILATERAL, box.x(),
                                  box.y(), box.z()));

  // Align both meshes
  {
    AABB mesh_aabb(Eigen::Vector3d(0, 0, 0), box);
    Eigen::Vector3d translate = surf_aabb.center() - mesh_aabb.center();
    std::cout << "mesh translation: " << translate << std::endl;
    for (int i = 0; i < mesh.GetNV(); i++) {
      real_t *v = mesh.GetVertex(i);
      v[0] += translate[0];
      v[1] += translate[1];
      v[2] += translate[2];
    }
  }

  int order = 2;
  H1_FECollection fec(order, 3);

  FiniteElementSpace h1_fespace(&mesh, &fec);
  FiniteElementSpace h1_fespace_vdim3(&mesh, &fec, mesh.SpaceDimension());

  GridFunction dist(&h1_fespace);
  GridFunction dist_grad(&h1_fespace_vdim3);

  std::cout << "NDofs: " << h1_fespace.GetNDofs() << std::endl;

  FunctionCoefficient tmd_fc([&](const Vector &coord) {
    // FIXME: use the feature and barycentric coordinate to determine the actual
    // thickness
    auto res = tmd.signed_distance(coord);
    return res.distance;
  });

  dist.ProjectCoefficient(tmd_fc);

  // Gradient of the distance function
  GradientGridFunctionCoefficient dist_grad_coeff(&dist);
  dist_grad.ProjectCoefficient(dist_grad_coeff);

  {
    // Paraview
    surf.Save("surf.vtk");

    ParaViewDataCollection dc("Solid", &mesh);
    dc.SetPrefixPath("ParaView");
    dc.SetLevelsOfDetail(order);
    dc.SetHighOrderOutput(true);
    dc.SetDataFormat(VTKFormat::BINARY);
    dc.RegisterField("dist", &dist);
    dc.RegisterField("dist_grad", &dist_grad);
    dc.Save();
  }

  {
    // glvis
    mesh.Save("mksolid.mesh");
    dist.Save("mksolid-dist.gf");
    dist_grad.Save("mksolid-dist_grad.gf");
  }

  return 0;
}

