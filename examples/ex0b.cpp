//                                MFEM Example 0
//
// Compile with: make ex0
//
// Sample runs:  ex0
//               ex0 -m ../data/fichera.mesh
//               ex0 -m ../data/square-disc.mesh -o 2
//
// Description: This example code demonstrates the most basic usage of MFEM to
//              define a simple finite element discretization of the Laplace
//              problem -Delta u = 1 with zero Dirichlet boundary conditions.
//              General 2D/3D mesh files and finite element polynomial degrees
//              can be specified by command line options.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

class MyIntegrator : public BilinearFormIntegrator {
public:
   // The coefficient must be a pointer. A reference cannot be null, so if we have
   // mutually exclusive options, we would need to define them all.
   Coefficient *Q;
   MyIntegrator(const IntegrationRule *ir) : BilinearFormIntegrator(ir) { }
   MyIntegrator(Coefficient &q, const IntegrationRule *ir) : BilinearFormIntegrator(ir) {
      Q = &q;
   }
   ~MyIntegrator() override {}
   void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Trans, DenseMatrix &elmat) override {

      int ndof = el.GetDof();
      int dim = el.GetDim(); // Segment, Surface, Volume
      int spaceDim = Trans.GetSpaceDim(); // 1D, 2D, 3D
      bool square = (dim == spaceDim);
      // Possible combinations
      // Segment: 1D, 2D, 3D
      // Surface: 2D, 3D
      // Volume: 3D
      // In general, dim <= spaceDim must be true

      DenseMatrix dshape(ndof, dim);
      DenseMatrix dshapedxt(ndof, spaceDim);
      elmat.SetSize(ndof);

      const IntegrationRule *ir = GetIntegrationRule(el, Trans);
      if (ir == NULL) {
         ir = &IntRules.Get(el.GetGeomType(), el.GetOrder());
      }

      std::cout << "el.GetGeomType: " << el.GetGeomType() << std::endl;
      std::cout << "el.GetOrder: " << el.GetOrder() << std::endl;
      std::cout << "ir->GetNPoints: " << ir->GetNPoints() << std::endl;

#if false
      for (int i = 0; i < ir->GetNPoints(); i++) {
         const IntegrationPoint &ip = ir->IntPoint(i);
         std::cout << "ip " << i << " " << ip.weight << " " << ip.x << " " << ip.y << " " << ip.z << std::endl;
      }
#endif

      for (int i = 0; i < ir->GetNPoints(); i++) {
         const IntegrationPoint &ip = ir->IntPoint(i);
         el.CalcDShape(ip, dshape); // Gradient of the shape function

         Trans.SetIntPoint(&ip);

         // Weigth to scale from reference to physical quantity (Jacobian determinant)
         // Actually, this function returns $ \sqrt{\lvert J^T J \rvert} $
         // For square matrices, it simply returns det(J), but in the case of non-square matrices,
         // it returns sqrt(det(J^T J)). (J^T J) produce a square matrix where the determinant is
         // defined, but causes the area to be squared, and that is why the sqrt is required.
         real_t w = Trans.Weight(); // element length/surface/volume

         // Integration point weight
         real_t gauss_w = ip.weight;

         // Jacobian = matrix of the element transformation
         // Adjugate = generalize the inverse to non-square matrices (AKA adjoint)
         // inv(J) = (1/det(J)) * adj(J)
         // det(J)*inv(J) = adj(J)

         // AdjugateJacobian = / adj(J),         if J is square
         //                    \ adj(J^t.J).J^t, otherwise

         // FIXME: Combined weight: why do we devide here by the matrix weight?
         real_t JxW = gauss_w / (square ? w : w*w*w);

         // Include in the weight the value of the diffusion coefficient (i.e. conductivity)
         // The conductivity can vary accord to space and time
         // If Q is null, then it is equivalent of Q = 1.0
         if (Q) {
            JxW *= Q->Eval(Trans, ip);
         }

         const DenseMatrix& jac = Trans.Jacobian();
         const DenseMatrix& adj = Trans.AdjugateJacobian();

         std::cout << "ip " << i << " " << ip.weight << " " << ip.x << " " << ip.y << " " << ip.z << std::endl;
         std::cout << "dim:" << dim << std::endl;
         std::cout << "spaceDim:" << spaceDim << std::endl;
         std::cout << "jac:" << jac.Height() << ", " << jac.Width() << std::endl;
         std::cout << "adj:" << adj.Height() << ", " << adj.Width() << std::endl;
         std::cout << "ndof:" << ndof << std::endl;
         std::cout << "trans.Weight():" << Trans.Weight() << std::endl;
         std::cout << "gauss_w:" << Trans.Weight() << std::endl;

         // dshapedxt = dshape * AdjugateJacobian : transform dshape to physical space
         // To transform a local gradient to physical, we have to multiply by J^-T (which is the transposed inverse of the Jacobian)
         // Because the inverse is undefined for non-square matrix, the adjugate is used instead. The adjugate contains the factor det(J).
         // The result is dshapedxt has t prefix for transposed, because we multiply by the non-transposed adjugate.
         // To get dshapedx itself, I guess we would need to do call
         // Mult(dshape, Trans.TransposeAdjugateJacobian(), dshapedx)
         Mult(dshape, Trans.AdjugateJacobian(), dshapedxt);

         // AAt += a * A * A^t
         // Compact equivalent to the product of grad_phi_i X grad_phi_j
         // Replaces two nested loops over the integration points (i, j)
         // Because dshapedxt includes the factor det(J), the product dshapedxt*dshapedx has det(J)^2 factor built-in.
         // We need to cancel one det(J), that is why we have w = ip.weight / w.
         // For non-square matrices, we have a factor of (w^3) to cancel, because (adj(J^t.J).J^t)^2
         AddMult_a_AAt(w, dshapedxt, elmat);

      }

   }
};

int main(int argc, char *argv[])
{

   // 1. Parse command line options.
   string mesh_file = "../data/star.mesh";
   int order = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");
   args.ParseCheck();

   // 2. Read the mesh from the given mesh file, and refine once uniformly.
   Mesh mesh(mesh_file);
   mesh.UniformRefinement();

   std::cout << "Mesh Dim: " << mesh.Dimension() << std::endl;
   std::cout << "Mesh SpaceDim: " << mesh.SpaceDimension() << std::endl;

   // 3. Define a finite element space on the mesh. Here we use H1 continuous
   //    high-order Lagrange finite elements of the given order.
   H1_FECollection fec(order, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fec);
   cout << "Number of unknowns: " << fespace.GetTrueVSize() << endl;

   // 4. Extract the list of all the boundary DOFs. These will be marked as
   //    Dirichlet in order to enforce zero boundary conditions.
   Array<int> boundary_dofs;
   fespace.GetBoundaryTrueDofs(boundary_dofs);

   // 5. Define the solution x as a finite element grid function in fespace. Set
   //    the initial guess to zero, which also sets the boundary conditions.
   GridFunction x(&fespace);
   x = 0.0;

   // 6. Set up the linear form b(.) corresponding to the right-hand side.
   ConstantCoefficient one(1.0);
   LinearForm b(&fespace);
   b.AddDomainIntegrator(new DomainLFIntegrator(one));
   b.Assemble();

   // 7. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new MyIntegrator(one, nullptr));
   a.Assemble();

   // 8. Form the linear system A X = B. This includes eliminating boundary
   //    conditions, applying AMR constraints, and other transformations.
   SparseMatrix A;
   Vector B, X;
   a.FormLinearSystem(boundary_dofs, x, b, A, X, B);

   // 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 200, 1e-12, 0.0);

   // 10. Recover the solution x as a grid function and save to file. The output
   //     can be viewed using GLVis as follows: "glvis -m mesh.mesh -g sol.gf"
   a.RecoverFEMSolution(X, b, x);
   x.Save("sol.gf");
   mesh.Save("mesh.mesh");

   return 0;
}
