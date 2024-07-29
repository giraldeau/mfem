#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace mfem;
using namespace std;

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
  Mesh mesh = Mesh::MakeCartesian2D(100, 100, Element::QUADRILATERAL, 1);
  mesh.EnsureNCMesh(true);

  int order = 1;
  int dim = mesh.Dimension();
  int sdim = mesh.SpaceDimension();
  H1_FECollection fec(order, dim);
  FiniteElementSpace fespace(&mesh, &fec);

  Vector c1({0.5, 0.5});
  double r1 = 0.25;
  FunctionCoefficient coeff2([&](const Vector &p) -> double {
    const double xc = p(0) - 0.5, yc = p(1) - 0.5;
    const double r = sqrt(xc*xc + yc*yc);
    double len = (r >= 0.2 && r <= 0.4) ? -1.0 : 1.0;
    cout << p(0) << " " << p(1) << " " << len << "\n";
    return len;
  });
  FunctionCoefficient coeff3([&](const Vector &p) -> double {
    const double xc = p(0) - 0.5, yc = p(1) - 0.5;
    const double r = sqrt(xc*xc + yc*yc);
    return std::tanh(2.0*(r-0.3));
  });

  GridFunction ls(&fespace);
  ls.ProjectCoefficient(coeff3);

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

  ParaViewDataCollection dc("MkSolid", &mesh);
  dc.SetPrefixPath("ParaView");
  dc.SetLevelsOfDetail(order);
  dc.SetHighOrderOutput(true);
  dc.SetDataFormat(VTKFormat::BINARY);
  dc.RegisterField("ls", &ls);
  dc.Save();

  return 0;
}
