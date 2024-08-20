#include "mfem.hpp"
#include <iostream>

#include <Eigen/Geometry>
#include "tmd/TriangleMeshDistance.h"

/*
 * Compute the distance field from a surface mesh.
 */

using namespace mfem;
using namespace std;

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
  Eigen::AlignedBox<double, 3> aabb;
  for (int i = 0; i < surf.GetNV(); i++) {
    real_t *v = surf.GetVertex(i);
    aabb.extend(Eigen::Vector3d(v));
  }

  std::cout << "AABB min:\n" << aabb.min() << "\n";
  std::cout << "AABB max:\n" << aabb.max() << "\n";

  // Build the triangle distance index
  tmd::TriangleMeshDistance tmd;
  tmd.construct(surf);

  // Create a mesh to evaluate the distance field

  Eigen::Vector3d box = aabb.sizes() * 1.2;
  Mesh mesh(Mesh::MakeCartesian3D(30, 3, 30, Element::QUADRILATERAL,
                                  box.x(), box.y(), box.z() * 10));

  int order = 3;
  H1_FECollection fec(order, 3);
  FiniteElementSpace fespace(&surf, &fec);
  GridFunction dist(&fespace);

  FunctionCoefficient tmd_fc([&](const Vector &coord){
    auto res = tmd.signed_distance(coord);
    return res.distance;
  });

  dist.ProjectCoefficient(tmd_fc);

  ParaViewDataCollection dc("Solid", &mesh);
  dc.SetPrefixPath("ParaView");
  dc.SetLevelsOfDetail(order);
  dc.SetHighOrderOutput(true);
  dc.SetDataFormat(VTKFormat::BINARY);
  dc.RegisterField("distance", &dist);
  dc.Save();

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
