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
   SparseMatrix C;
   Vector gap;
};

std::string generateFilename(const std::string& templateStr, int number) {
   std::size_t bufferSize = templateStr.size() + 10; // Extra space for the number
   std::vector<char> buffer(bufferSize);
   std::snprintf(buffer.data(), buffer.size(), templateStr.c_str(), number);
   return std::string(buffer.data());
}

void BuildContactMatrix(struct ContactData& state, FiniteElementSpace& fespace,
   FunctionCoefficient &dist, VectorFunctionCoefficient &dir,
   std::vector<int> &row_dofs, Vector &lm_prev)
{
   int dim = fespace.GetVDim();

   // Count the number of contact constraints
   // The constraint is active if the level-set distance is negative
   row_dofs.clear();
   std::unordered_map<int, int> row_map;
   for (int i = 0; i < fespace.GetNBE(); ++i)
   {
      int attr = fespace.GetBdrAttribute(i);
      // FIXME: restrict to attribute 3 (top edge)
      if (attr == 3)
      {
         Array<int> dofs;
         fespace.GetBdrElementDofs(i, dofs);
         ElementTransformation *Tr = fespace.GetBdrElementTransformation(i);
         const FiniteElement *fe = fespace.GetBE(i);
         const IntegrationRule& nodes = fe->GetNodes();
         // FIXME: difference between integration point and the DOF
         // I think here IntegrationRule is used by convenience, because in
         // general nodes and quadrature points are not related.
         for (int j = 0; j < dofs.Size(); ++j)
         {
            int node = dofs[j];
            const IntegrationPoint &ip = nodes.IntPoint(j);
            real_t gap_val = dist.Eval(*Tr, ip);
            real_t lm_val = lm_prev[node];
            if (gap_val < 0 && lm_val >= 0.0)
            {
               const auto &iter = row_map.find(node);
               if (iter == row_map.end())
               {
                  row_map[node] = row_dofs.size();
                  row_dofs.push_back(node);
                  std::cout << "dof: " << dofs[j] << " " << gap_val << " " << lm_val << std::endl;
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
   state.C = SparseMatrix(row_dofs.size(), fespace.GetTrueVSize());
   state.gap = Vector(row_dofs.size());
   state.gap = 0.0;
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

            real_t lm_val = lm_prev[node];

            const IntegrationPoint &ip = nodes.IntPoint(j);
            Tr->SetIntPoint(&ip);
            real_t gap_val = dist.Eval(*Tr, ip);
            if (gap_val < 0 && lm_val >= 0.0)
            {
               Vector nor(3);
               dir.Eval(nor, *Tr, ip);

               int row = row_map.at(node);
               state.gap[row] += gap_val;
               std::cout << "row " << row << " " << gap_val << " norm: " << nor[0] << "," << nor[1] << std::endl;

               for (int d = 0; d < dim; ++d)
               {
                  int vdof = fespace.DofToVDof(node, d);
                  std::cout << "   " << node << " " << d << " " << vdof << std::endl;
                  state.C.Add(row, vdof, nor[d]);
               }
            }
         }
      }
   }

   state.C.Finalize();
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

   ParaViewDataCollection dc("Contact03");
   dc.SetPrefixPath("ParaView");
   dc.SetMesh(mesh);

   int max_iter = 4;
   struct ContactData state{};
   Vector lm_prev(fespace->GetNV());
   // load step
   int nstep = 10;
   for (int step = 0; step < nstep; ++step) {

      std::cout << "Step " << step << std::endl;

      /*
       * Contact
       */

      real_t center_y = 10.0 + (-3.0 * step) / (nstep - 1);
      CircleDistance circle(5, center_y, 4);
      FunctionCoefficient ls_coefficient(circle);
      VectorFunctionCoefficient lsv_coefficient(dim, circle);

      GridFunction x_ls(fespace_ls);
      x_ls.ProjectCoefficient(ls_coefficient);

      GridFunction x_dir(fespace);
      x_dir.ProjectCoefficient(lsv_coefficient);

      std::vector<int> row_dofs;
      std::vector<int> row_dofs_prev;
      GridFunction x(fespace);

      lm_prev = 0.0;
      for (int iter = 0; iter < max_iter; ++iter) {
         std::cout << "ITER " << iter << std::endl;
         // There is no forces on the RHS
         LinearForm b(fespace);
         cout << "r.h.s. ... " << std::endl;
         b.Assemble();

         // reset solution
         x = 0.0;

         ConstantCoefficient lambda_func(lambda);
         ConstantCoefficient mu_func(mu);

         BilinearForm a(fespace);
         a.AddDomainIntegrator(new ElasticityIntegrator(lambda_func,mu_func));

         cout << "matrix ... " << std::endl;
         a.Assemble();

         SparseMatrix A;
         Vector B, X;
         a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
         cout << "done." << endl;

         cout << "Size of linear system: " << A.Height() << endl;

         BuildContactMatrix(state, *fespace, ls_coefficient, lsv_coefficient, row_dofs, lm_prev);
         {
            std::ofstream mat(generateFilename("C_%03d.mat", iter));
            state.C.PrintMatlab(mat);
            std::ofstream gap(generateFilename("gap_%03d.dat", iter));
            state.gap.Print(gap);
         }

         // Solver
         GSSmoother M(A);
         SchurConstrainedSolver solver(A, state.C, M);
         solver.SetConstraintRHS(state.gap);
         solver.SetRelTol(1e-5);
         solver.SetMaxIter(2000);
         solver.SetPrintLevel(1);
         solver.Mult(B, X);

         Vector lm;
         //lm_prev = 0.0;
         solver.GetMultiplierSolution(lm);
         for (int i = 0; i < row_dofs.size(); ++i) {
            std::cout << "LM " << row_dofs[i] << " " << lm[i] << std::endl;
            lm_prev[row_dofs[i]] = lm[i];
         }

         a.RecoverFEMSolution(X, b, x);

         // end criteria: no change in active set
         if (row_dofs == row_dofs_prev) {
            break;
         }
         row_dofs_prev = row_dofs;
      }

      GridFunction *nodes = mesh->GetNodes();
      *nodes += x;

      // save solution only at the end of the iteration
      dc.SetCycle(step);
      dc.SetTime((real_t)step);
      dc.Save();

      // Move nodes to their original position
      *nodes -= x;

   }

   delete fespace;
   delete fec;
   delete mesh;

   return 0;
}
