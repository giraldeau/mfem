#include "mfem.hpp"
#include <iostream>

#include <Eigen/Geometry>
#include "tmd/TriangleMeshDistance.h"

/*
 * Compute the distance field from a surface mesh.
 */

using namespace mfem;
using namespace std;

typedef Eigen::AlignedBox<double, 3> AABB;

// Unused at the moment
void DiffuseField(GridFunction &field, int smooth_steps)
{
   // Setup the Laplacian operator
   BilinearForm *Lap = new BilinearForm(field.FESpace());
   Lap->AddDomainIntegrator(new DiffusionIntegrator());
   Lap->Assemble();
   Lap->Finalize();

   // Setup the smoothing operator
   DSmoother *S = new DSmoother(0,1.0,smooth_steps);
   S->iterative_mode = true;
   S->SetOperator(Lap->SpMat());

   Vector tmp(field.Size());
   tmp = 0.0;
   S->Mult(tmp, field);

   delete S;
   delete Lap;
}

void ComputeMeshAABB(const Mesh& mesh, AABB& aabb)
{
  for (int i = 0; i < mesh.GetNV(); i++) {
    aabb.extend(Eigen::Vector3d(mesh.GetVertex(i)));
  }
}
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
  std::cout << "Input mesh AABB min:\n" << surf_aabb.min() << "\n" << surf_aabb.max() << "\n";

  // Build the triangle distance index
  tmd::TriangleMeshDistance tmd;
  tmd.construct(surf);

  // Create a mesh to evaluate the distance field
  // FIXME: uniform grid is not scalable. Implement AMR
  Eigen::Vector3d box = 2.0 * (surf_aabb.sizes() + Eigen::Vector3d(0, 0, 10));
  Mesh mesh(Mesh::MakeCartesian3D(10, 10, 10, Element::QUADRILATERAL,
                                  box.x(), box.y(), box.z()));

  {
    AABB mesh_aabb;
    ComputeMeshAABB(mesh, mesh_aabb);
    std::cout << "BEFORE Grid mesh AABB min:\n" << mesh_aabb.min() << "\n" << mesh_aabb.max() << "\n";
  }

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

  // Check bounding box
  {
    AABB mesh_aabb;
    ComputeMeshAABB(mesh, mesh_aabb);
    std::cout << "AFTER Grid mesh AABB min:\n" << mesh_aabb.min() << "\n" << mesh_aabb.max() << "\n";
  }

  int order = 3;
  H1_FECollection fec(order, 3);
  FiniteElementSpace fespace(&mesh, &fec);
  GridFunction dist(&fespace);

  std::cout << "NDofs: " << fespace.GetNDofs() << std::endl;

  FunctionCoefficient tmd_fc([&](const Vector &coord){
    auto res = tmd.signed_distance(coord);
    return res.distance;
  });

  dist.ProjectCoefficient(tmd_fc);

  {
    // Paraview
    ParaViewDataCollection dc("Solid", &mesh);
    dc.SetPrefixPath("ParaView");
    dc.SetLevelsOfDetail(order);
    dc.SetHighOrderOutput(true);
    dc.SetDataFormat(VTKFormat::BINARY);
    dc.RegisterField("distance", &dist);
    dc.Save();
  }


  {
    // glvis
    mesh.Save("mksolid.mesh");
    dist.Save("mksolid.gf");
  }

  return 0;
}


/*
  //DiffuseField(ls, 2);

  ConstantCoefficient one(1.0);
  DiffusionIntegrator integ(one);
  FiniteElementSpace flux_fespace(&mesh, &fec, dim);
  ZienkiewiczZhuEstimator estimator(integ, ls, flux_fespace);
  estimator.SetAnisotropic();

  double max_elem_error = 5.0e-3;
  ThresholdRefiner refiner(estimator);
  //refiner.SetTotalErrorFraction(0.7);
  refiner.SetTotalErrorFraction(0.0);
  refiner.SetLocalErrorGoal(max_elem_error);

  int max_dofs = 1000;
  int cdofs = fespace.GetTrueVSize();
  int it = 0;
  while (cdofs < max_dofs && it < 10)
  {
     cdofs = fespace.GetTrueVSize();
     cout << "\nAMR iteration " << it << endl;
     cout << "Number of dofs: " << cdofs << endl;

     ls.ProjectCoefficient(coeff2);
     //DiffuseField(ls, 2);

     refiner.Reset();
     refiner.Apply(mesh);

     fespace.Update();
     ls.SetSize(fespace.GetVSize());
     ls.ProjectCoefficient(coeff2);

     if (refiner.Stop())
     {
        cout << "Stopping criterion satisfied. Stop." << endl;
        break;
     }

     it++;
  }
  */
