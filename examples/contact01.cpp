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
#include <unordered_set>

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

   // Method to compute the direction vector at a given location
   void operator()(const Vector &x, Vector &dir) const
   {
      double dx = x(0) - centerX;
      double dy = x(1) - centerY;
      double d = std::sqrt(dx * dx + dy * dy);
      // what to use for undefined direction?
      dir = 0.0;
      if (d > 1e-12)
      {
         dir(0) = -1.0 * dx / d;
         dir(1) = -1.0 * dy / d;
      }
   }

private:
   real_t centerX;
   real_t centerY;
   real_t radius;
};

struct ContactData {
   SparseMatrix *C;
   Vector *gap;
};

void BuildContactMatrix(struct ContactData& state, FiniteElementSpace& fespace, FunctionCoefficient &dist, VectorFunctionCoefficient &dir)
{
   int dim = fespace.GetVDim();

   // Count the number of contact constraints
   // The constraint is active if the level-set distance is negative
   std::unordered_map<int, int> row_map;
   std::vector<int> row_dofs;
   for (int i = 0; i < fespace.GetNBE(); ++i)
   {
      int attr = fespace.GetBdrAttribute(i);

      // FIXME: restrict to attribute 3 (top edge)
      if (attr == 3)
      {
         Array<int> dofs;
         Array<int> vdofs;
         fespace.GetBdrElementDofs(i, dofs);

         ElementTransformation *Tr = fespace.GetBdrElementTransformation(i);
         const FiniteElement *fe = fespace.GetBE(i);
         const IntegrationRule& nodes = fe->GetNodes();
         // FIXME: difference between integration point and the DOF
         // I think here IntegrationRule is used by convenience, because in
         // general nodes and quadrature points are not related.
         for (int j = 0; j < dofs.Size(); ++j)
         {
            const IntegrationPoint &ip = nodes.IntPoint(j);
            real_t d = dist.Eval(*Tr, ip);
            bool sign = (d < 0);
            if (d < 0)
            {
               const auto &iter = row_map.find(dofs[j]);
               if (iter == row_map.end())
               {
                  row_map[dofs[j]] = row_dofs.size();
                  row_dofs.push_back(dofs[j]);
                  std::cout << "dof: " << dofs[j] << " " << d << std::endl;
               }
            }
         }
      }
   }

   // Display the DOFs in order
   std::cout << "n_rows: " << row_dofs.size() << std::endl;
   for (const auto& dof : row_dofs)
   {
      const auto &row = row_map[dof];
      std::cout << row << " " << dof << std::endl;
   }

   // Compute the actual constraint matrix and RHS
   state.C = new SparseMatrix(row_dofs.size(), fespace.GetTrueVSize());
   state.gap = new Vector(row_dofs.size());
   *state.gap = 0.0;
   for (int i = 0; i < fespace.GetNBE(); ++i)
   {
      int attr = fespace.GetBdrAttribute(i);
      if (attr == 3)
      {
         ElementTransformation *Tr = fespace.GetBdrElementTransformation(i);
         const FiniteElement *fe = fespace.GetBE(i);
         const IntegrationRule& nodes = fe->GetNodes();

         Array<int> dofs;
         fespace.GetBdrElementDofs(i, dofs);
         for (int j = 0; j < dofs.Size(); ++j)
         {
            // DOF is the node number, but there are 2 dof per node here
            // With Ordering::byNODES, the VDofs of a node are not consecutive
            int node = dofs[j];


            const IntegrationPoint &ip = nodes.IntPoint(j);
            Tr->SetIntPoint(&ip);
            real_t val = dist.Eval(*Tr, ip);
            if (val < 0)
            {
               Vector nor(3);
               dir.Eval(nor, *Tr, ip);

               int row = row_map.at(node);
               (*state.gap)[row] += val;
               std::cout << "row " << row << " " << val << " norm: " << nor[0] << "," << nor[1] << std::endl;

               for (int d = 0; d < dim; ++d)
               {
                  int vdof = fespace.DofToVDof(node, d);
                  std::cout << "   " << node << " " << d << " " << vdof << std::endl;
                  state.C->Add(row, vdof, nor[d]);
               }
            }
         }
      }
   }

   state.C->Finalize();
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int order = 1;
   bool static_cond = false;
   bool visualization = 1;
   real_t lambda = 1.0;
   real_t mu = 1.0;
   int refinements = 0;

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
   args.AddOption(&refinements, "-r", "--refinements", "Number of refinement levels");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);
   Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(10, 10, Element::TRIANGLE, false, 10.0, 5.0));
   int dim = mesh->Dimension();
   for (int r = 0; r < refinements; ++r) {
      mesh->UniformRefinement();
   }

   FiniteElementCollection *fec = new H1_FECollection(order, dim);
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec, dim);
   FiniteElementSpace *fespace_ls = new FiniteElementSpace(mesh, fec);
   mesh->EnsureNodes(); // make sure the node array exists
   mesh->SetNodalFESpace(fespace); // required to move the nodes
   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl << "Assembling: " << std::endl;
   {
      std::ofstream ofs("contact.mesh");
      mesh->Print(ofs);
   }

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
   dc.SetHighOrderOutput(false);
   dc.SetLevelsOfDetail(0);
   dc.SetMesh(mesh);

   /*
    * Contact
    */
   real_t dy = 2.0;
   CircleDistance circle(5, 9 - dy, 4);
   FunctionCoefficient ls_coefficient(circle);
   VectorFunctionCoefficient lsv_coefficient(dim, circle);

   GridFunction x_ls(fespace_ls);
   x_ls.ProjectCoefficient(ls_coefficient);

   GridFunction x_dir(fespace);
   x_dir.ProjectCoefficient(lsv_coefficient);

   struct ContactData state{};
   BuildContactMatrix(state, *fespace, ls_coefficient, lsv_coefficient);
   {
      std::ofstream mat("C.mat");
      state.C->PrintMatlab(mat);
      std::ofstream gap("gap.mat");
      state.gap->Print(gap);
   }
   GSSmoother M(A);
   SchurConstrainedSolver * solver = new SchurConstrainedSolver(A, *state.C, M);
   solver->SetConstraintRHS(*state.gap);
   solver->SetRelTol(1e-5);
   solver->SetMaxIter(2000);
   solver->SetPrintLevel(1);
   solver->Mult(B, X);

   Vector lm;
   solver->GetMultiplierSolution(lm);
   lm.Print();

   a->RecoverFEMSolution(X, *b, x);
   GridFunction *nodes = mesh->GetNodes();
   *nodes += x;
   fespace_ls->Update();

   {
      ofstream mesh_ofs("displaced.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      ofstream sol_ofs("sol.mesh");
      GridFunction x_ls_end(fespace_ls);
      x_ls_end.ProjectCoefficient(ls_coefficient);
      x_ls_end.Save(sol_ofs);

      dc.RegisterField("x_ls_end", &x_ls_end);
      dc.Save();
   }

   // 16. Free the used memory.
   delete a;
   delete b;
   if (fec)
   {
      delete fespace;
      delete fec;
   }
   delete mesh;
   delete solver;

   return 0;
}
