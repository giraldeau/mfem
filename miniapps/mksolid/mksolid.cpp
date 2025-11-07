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

void ComputeMeshAABB(const Mesh &mesh, AABB &aabb) {
  for (int i = 0; i < mesh.GetNV(); i++) {
    aabb.extend(Eigen::Vector3d(mesh.GetVertex(i)));
  }
}

int main(int argc, char *argv[]) {

  int refinement = 0;
  std::string outname = "solid";
  std::string mesh_name;
  int simplex = 0;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_name, "-m", "--mesh", "Input mesh (surface in 3D)");
  args.AddOption(&refinement, "-r", "--refinements", "Number of uniform refinements");
  args.AddOption(&simplex, "-s", "--simplex", "Output grid as simplexes");
  args.AddOption(&outname, "-o", "--output", "Output file basename");

  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return 1;
  }
  args.PrintOptions(std::cout);

  Mesh surf(mesh_name);

  // Refine the surface mesh. This is useful only for curved meshes. If we subdivide triangles,
  // the subdivided triangles lies in the original triangle plane, therefore it has no effect
  // on the level-set definition.
  if (refinement > 0) {
    for (int i = 0; i < refinement; i++) {
      surf.UniformRefinement();
    }
  }

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
  double margin = surf_aabb.diagonal().norm() * 0.25;
  surf_aabb.min().array() -= margin;
  surf_aabb.max().array() += margin;
  Eigen::Vector3d box = surf_aabb.sizes();
  Mesh mesh(Mesh::MakeCartesian3D(1, 1, 1, Element::QUADRILATERAL,
                                  box.x(), box.y(), box.z()));

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

  if (simplex == 1) {
    mesh = Mesh::MakeSimplicial(mesh);
    mesh.Finalize(true);
  }

  int order = 1;
  H1_FECollection fec(order, 3);
  FiniteElementSpace h1_fespace(&mesh, &fec);
  GridFunction dist(&h1_fespace);
  GridFunction sigmoid(&h1_fespace);

  std::cout << "NDofs: " << h1_fespace.GetNDofs() << std::endl;

  double thick = 1.0;
  double half_thick = thick * 0.5;
  FunctionCoefficient tmd_fc([&](const Vector &coord) {
    // FIXME: use the feature and barycentric coordinate
    // to determine the actual thickness
    auto res = tmd.signed_distance(coord);

    // simple distance to level-set
    //return res.distance;

    // constant mid-plane thickness
    return std::abs(res.distance) - half_thick;
  });

  FunctionCoefficient tmd_fc_sigmoid([&](const Vector &coord) {
    // FIXME: use the feature and barycentric coordinate to determine the actual
    // thickness
    auto res = tmd.signed_distance(coord);

    // simple distance to level-set
    //return res.distance;

    // constant mid-plane thickness
    real_t v1 = std::abs(res.distance) - half_thick;
    real_t v2 = (cosh(10*v1));
    real_t v3 = 1.0 / (v2*v2);
    return v3;
  });

  // Initial state
  dist.ProjectCoefficient(tmd_fc);
  sigmoid.ProjectCoefficient(tmd_fc_sigmoid);

  // Refinement step: detect elements that contains the sigmoid of the level-set
  LpErrorEstimator estimator(2, tmd_fc_sigmoid, sigmoid);
  ThresholdRefiner refiner(estimator);
  refiner.SetMaxElements(200000);

  // Prepare ParaView Data Collection
  ParaViewDataCollection dc("Solid", &mesh);
  dc.SetPrefixPath("ParaView");
  dc.SetDataFormat(VTKFormat::BINARY);
  dc.RegisterField("dist", &dist);
  dc.RegisterField("sigmoid", &sigmoid);

  // Also save the surface mesh to VTK to superimpose to the domain
  {
    std::ofstream out("surf.vtk");
    surf.PrintVTK(out);
  }

  // AMR Loop
  for (int it = 0; ; it++) {
    dc.SetTime(it);
    dc.SetCycle(it);
    dc.Save();

    int cdofs = h1_fespace.GetTrueVSize();
    cout << "AMR iteration " << it << " ndofs: " << cdofs << endl;

    refiner.Apply(mesh);
    if (refiner.Stop()) {
      cout << "Stopping criterion satisfied. Stop." << endl;
      break;
    }

    h1_fespace.Update();
    sigmoid.Update();
    dist.Update();

    dist.ProjectCoefficient(tmd_fc);
    sigmoid.ProjectCoefficient(tmd_fc_sigmoid);
  }

  {
    // glvis
    mesh.Save("mksolid.mesh");
    dist.Save("mksolid-dist.gf");
  }

  return 0;
}

