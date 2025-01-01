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

// Unused at the moment
void DiffuseField(GridFunction &field, int smooth_steps) {
  // Setup the Laplacian operator
  BilinearForm *Lap = new BilinearForm(field.FESpace());
  Lap->AddDomainIntegrator(new DiffusionIntegrator());
  Lap->Assemble();
  Lap->Finalize();

  // Setup the smoothing operator
  DSmoother *S = new DSmoother(0, 1.0, smooth_steps);
  S->iterative_mode = true;
  S->SetOperator(Lap->SpMat());

  Vector tmp(field.Size());
  tmp = 0.0;
  S->Mult(tmp, field);

  delete S;
  delete Lap;
}

void ComputeMeshAABB(const Mesh &mesh, AABB &aabb) {
  for (int i = 0; i < mesh.GetNV(); i++) {
    aabb.extend(Eigen::Vector3d(mesh.GetVertex(i)));
  }
}

// Update point positions using Runge-Kutta method
void AdvecPoint(Mesh &mesh, const GridFunction &field, const Eigen::Vector3d &pos,
               const double &step, Eigen::Vector3d &res) {

  // Clunky for one point
  Array<int> elem_ids;
  Array<IntegrationPoint> ips;
  DenseMatrix mat(3, 1);
  mat(0, 0) = pos.x();
  mat(1, 0) = pos.y();
  mat(2, 0) = pos.z();

  int found = mesh.FindPoints(mat, elem_ids, ips);

  cout << found << endl;
  if (found == 0) {
    cout << "YUCK";
    return;
  }

  int elem_id = elem_ids[0];
  IntegrationPoint ip = ips[0];
  Element *el = mesh.GetElement(elem_id);
  ElementTransformation *T = mesh.GetElementTransformation(elem_id);

  Vector gradient;
  T->SetIntPoint(&ip);
  field.GetGradient(*T, gradient);

  Eigen::Vector3d grad_vector(gradient[0], gradient[1], gradient[2]);
  res = pos + step * grad_vector;

  // move the integration point
  Eigen::Vector3d ipv_curr(ip.x, ip.y, ip.z);
  Eigen::Vector3d ipv_next;
  ipv_next = ipv_curr + step * grad_vector;
  ip.Set3(ipv_next.x(), ipv_next.y(), ipv_next.z());

  // DEBUG
  if (0) {
    cout << "MovePoint grad: " << grad_vector << endl;
    cout << "MovePoint pos: " << pos << endl;
    cout << "MovePoint res: " << res << endl;
    cout << "MovePoint ip1: " << ipv_curr << endl;
    cout << "MovePoint ip2: " << ipv_next << endl;
  }

  // FIXME: it is invalid to apply the physical gradient to the integration point
  // We must apply the value of CalcDShape directly?
  // Must change the element if the integration point is outside

  // k1 = gradientField(point);
  // k2 = gradientField(point + stepSize * k1 / 2.0);
  // k3 = gradientField(point + stepSize * k2 / 2.0);
  // k4 = gradientField(point + stepSize * k3);
  // point += stepSize * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;


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

  // Experiment: take a particle and move it in the gradient field
  int n = 100;
  std::vector<Eigen::Vector3d> line(n);

  // starting point
  Eigen::Vector3d p0(0, 0, 0);
  Eigen::Vector3d p1;

  for (int i = 0; i < 10; i++) {
    AdvecPoint(mesh, dist, p0, 1e-3, p1);
    cout << "(" << p0.x() << "," << p0.y() << "," << p0.z() << ") "
         << "(" << p1.x() << "," << p1.y() << "," << p1.z() << ")" << endl;
    p0 = p1;
  }

  // Create the volume mesh from the level-set

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
