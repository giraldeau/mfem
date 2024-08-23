#include "mfem.hpp"
#include <iostream>

using namespace mfem;
using namespace std;

void Barycentric(Vector a, Vector b, Vector c, Vector p, double &u, double &v, double &w)
{
    Vector v0(a.Size()), v1(a.Size()), v2(a.Size());
    subtract(b, a, v0);
    subtract(c, a, v1);
    subtract(p, a, v2);
    double d00 = v0 * v0;
    double d01 = v0 * v1;
    double d11 = v1 * v1;
    double d20 = v2 * v0;
    double d21 = v2 * v1;
    double denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0 - v - w;
}

double TriArea2D(double x1, double y1,
                       double x2, double y2,
                       double x3, double y3) {

    return (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
}

int TestPointTriangle(Vector p, Vector a, Vector b, Vector c) {
    double u, v, w;
    Barycentric(a, b, c, p, u, v, w);
    return v >= 0.0 && w >= 0.0 && (v + w) <= 1.0;
}

int IsConvexQuad(Vector a, Vector b, Vector c, Vector d)
{
    Vector bd(a.Size()), ba(a.Size()), bc(a.Size());
    Vector bda(a.Size()), bdc(a.Size());

    subtract(d, b, bd);
    subtract(a, b, ba);
    subtract(c, b, bc);

    bd.cross3D(ba, bda);
    bd.cross3D(bc, bdc);
    if ((bda * bdc) >= 0.0)
        return 0;

    Vector ac(a.Size()), ad(a.Size()), ab(a.Size());
    Vector acd(a.Size()), acb(a.Size());

    subtract(c, a, ac);
    subtract(d, a, ad);
    subtract(b, a, ab);

    ac.cross3D(ad, acd);
    ac.cross3D(ab, acb);

    return acd * acb < 0.0;
}

real_t simple_coeff(const Vector &x)
{
    return x[0] + x[1]*x[1] + x[2]*x[2]*x[2] - 1;
}

real_t scalar_field(const Vector &x)
{
  const int dim = x.Size();
  if (dim == 2)
  {
    return x[0] + x[1]*x[1];
  }
  else
  {
    return x[0] + x[1]*x[1];
  }
}

int main() {

    // how to model a simple point in 2D/3D ?
    {
        Mesh m; // points is an array of vertex
        Point p; // reference to the mesh entity
        Vector u({1, 2, 3}); // a point as a small vector?
        Vertex v({2, 2}); // dim 3 mem by default
    }

    {
        Vector u({2, 2});
        Vector v({4, -4});

        cout << "norm1  : " << u.Norml1() << endl; // sum
        cout << "norm2  : " << u.Norml2() << endl; // length
        cout << "norminf: " << u.Normlinf() << endl; // max
        cout << "normp: " << u.Normlp(3) << endl; // from norm1 to norminf
        for (int i = 1; i < 1E6; i *= 2) {
            cout << "normp: " << i << " " << u.Normlp(i) << endl;
        }

        cout << "dot : " << u * v << endl;
        cout << "add: " ;
        Vector res(u.Size());
        add(u, v, res);
        res.Print();
        cout << endl;
        subtract(u, v, res);
        cout << "sub: ";
        res.Print();
        cout << endl;
    }

    {
        Vector p1({2, 4});
        Vector p2({6, 10});
        Vector p3({10, 2});
        Vector p4({4, 5});
        Vector p5({10, 7});
        printf("isInside %d\n", TestPointTriangle(p4, p1, p2, p3));
        printf("isInside %d\n", TestPointTriangle(p5, p1, p2, p3));
    }

    {
        Vector v1({1, 2, 3});
        Vector v2({4, 5, 6});
        Vector res(v1.Size());
        v1.cross3D(v2, res);
        res.Print();
    }

    {
        Vector a1({1, 5, 1}); // convex
        Vector a2({3, 4, 1}); // concave
        Vector b({3, 8, 1});
        Vector c({6, 3, 1});
        Vector d({2, 2, 1});

        cout << "IsConvexQuad: " << IsConvexQuad(a1, b, c, d) << endl;
        cout << "IsConvexQuad: " << IsConvexQuad(a2, b, c, d) << endl;
    }


    {
        // Create a mesh with basic elements

        Mesh mesh(2, 0, 0, 0, 3);
        int v1 = mesh.AddVertex(0, 0, 0);
        int v2 = mesh.AddVertex(0, 1, 0);
        int v3 = mesh.AddVertex(1, 1, 0);
        int v4 = mesh.AddVertex(1, 0, 0);

        int v5 = mesh.AddVertex(0, 0, 0);
        int v6 = mesh.AddVertex(0, 1, 0);
        int v7 = mesh.AddVertex(1, 1, 0);
        int v8 = mesh.AddVertex(1, 0, 0);

        int q1 = mesh.AddQuad(v1, v2, v3, v4);
        mesh.AddBdrSegment(v1, v2);
        mesh.AddBdrSegment(v2, v3);
        mesh.AddBdrSegment(v3, v4);
        mesh.AddBdrSegment(v4, v1);
        mesh.Finalize(false, true);
        //mesh.UniformRefinement();

        Mesh mesh2(Mesh::MakeSimplicial(mesh));
        //mesh2.UniformRefinement();

        Array<int> el;
        el.Append(1);
        mesh2.GeneralRefinement(el);

        ofstream ofs("quad.mesh");
        mesh2.Print(ofs);
    }


    {
        // Experiment with finite element evaluation
        Mesh mesh(Mesh::MakeCartesian2D(1, 1, Element::QUADRILATERAL));

        // mesh.SetVertices(Vector({
        //     -1.0, 2.0, 0.0, 3.0, // x
        //     -1.0, 0.0, 2.0, 3.0, // y
        // }));

        H1_FECollection fec(1, 2);
        FiniteElementSpace fes(&mesh, &fec);
        GridFunction gf(&fes);
        gf.SetFromTrueDofs(Vector({
            10, 10, 10, 9
        }));


        {
            mesh.Save("simplequad.mesh");
            gf.Save("simplequad.gf");
        }

        // Evaluate at point
        Element *el = mesh.GetElement(0);
        const FiniteElement *fe = fes.GetFE(0);
        ElementTransformation *T = fes.GetElementTransformation(0);

        // obtain the global coordinates from the parametric coordinates on the element
        Vector p1(fe->GetDim());
        IntegrationPoint ip1;
        ip1.Set2(1, 1);
        T->Transform(ip1, p1);
        p1.Print();

        IntegrationPoint ip2;
        T->TransformBack(p1, ip2);
        cout << ip2.x << " " << ip2.y << endl;

        IntegrationPoint ip3;
        ip3.Set2(1.0, 1.0);
        // GetValue uses parametric element space, not global space!
        real_t v1 = gf.GetValue(0, ip3);
        cout << v1 << endl;


        Vector grad(fe->GetDim());
        T->SetIntPoint(&ip3);
        gf.GetGradient(*T, grad);
        cout << "gradient" << endl;
        grad.Print();

    }

    {
        // Experiment with finite element evaluation
        Mesh mesh(Mesh::MakeCartesian3D(10, 10, 10, Element::QUADRILATERAL));

        H1_FECollection fec(1, 3);
        FiniteElementSpace fes(&mesh, &fec);
        FunctionCoefficient simple(simple_coeff);
        GridFunction gf(&fes);
        gf.ProjectCoefficient(simple);

        {
            mesh.Save("simplecube.mesh");
            gf.Save("simplecube.gf");
        }

        GeometryRefiner refiner;

    }

    {
      // Compute the gradient of a scalar field
      // https://github.com/mfem/mfem/issues/865
      int order = 2;
      Mesh mesh(Mesh::MakeCartesian3D(10, 10, 10, Element::QUADRILATERAL));

      H1_FECollection h1_fec(order, mesh.Dimension());

      // Using Nedelec element we obtain implicitely a vector solution
      ND_FECollection nd_fec(order, mesh.Dimension());

      // Q: Can we use H1 elements to compute the gradient using GradientInterpolator?
      // A: Yes, but I think the vectors have limited continuity accros element boundaries.
      // We also need to specify vector dim vdim=2 to the FiniteElementSpace
      // Weird: we have to use at least 2nd order (no matter the order of the
      // finite element scalar field), otherwise the recovered solution is wrong.
      //H1_FECollection nd_fec(order, mesh.Dimension());

      // Weird: it looks like the solution obtained with RT element is inverted?!
      //RT_FECollection nd_fec(order, mesh.Dimension());

      FiniteElementSpace h1_fes(&mesh, &h1_fec);
      FiniteElementSpace nd_fes(&mesh, &nd_fec);

      FunctionCoefficient ls(scalar_field);

      GridFunction x(&h1_fes); // field
      GridFunction dx(&nd_fes); // field gradient for DiscreteLinearOperator
      GridFunction dx2(&nd_fes); // field gradient for GradientGridFunctionCoefficient

      cout << "scalar size: " << x.Size() << endl;
      cout << "gradient size: " << dx.Size() << endl;

      x.ProjectCoefficient(ls);

      // Method 1: use DiscreteLinearOperator to calculate the gradient
      DiscreteLinearOperator grad(&h1_fes, &nd_fes);
      grad.AddDomainInterpolator(new GradientInterpolator);
      grad.Assemble();
      grad.Finalize();
      grad.Mult(x, dx);

      mesh.Save("field.mesh");
      x.Save("x.gf");
      dx.Save("dx.gf");

      // Method 2: Project the gradient on another grid
      GradientGridFunctionCoefficient grad_coeff(&x);
      dx2.ProjectCoefficient(grad_coeff);

      // Method 3: compute the gradient inside the element using GetGradient().
      // Q: is this equivalent to do the interpolation of the derivative of the shape function?
      // A: yes, using CalcDShape(), we get the partial derivatives of each shape function
      Vector gradient;
      int elem_id = 64;
      Element *elem = mesh.GetElement(elem_id);
      ElementTransformation *T = mesh.GetElementTransformation(elem_id);

      Vector center;
      //Vector barycenter2;
      IntegrationPoint barycenter2;
      // barycentric coordinates to physical center
      const IntegrationPoint &barycenter = Geometries.GetCenter(elem->GetGeometryType());
      T->Transform(barycenter, center);

      // TransformBack is not trivial, uses newton iteration with tolerance
      // But why is this necessary?
      T->TransformBack(center, barycenter2);

      cout << "Element center barycentric (orig): " << barycenter.x << " "
           << barycenter.y << " " << barycenter.z << endl;
      cout << "Element center physical          : " << center[0] << " " << center[1] << " " << center[2] << endl;
      cout << "Element center barycentric (back): " << barycenter2.x << " "
           << barycenter2.y << " " << barycenter2.z << endl;

      T->SetIntPoint(&barycenter);
      x.GetGradient(*T, gradient);

      // Method 4: Compute the gradient manually
      double eps = 1e-6;
      Vector p1 = center;
      Vector p2({center[0]+eps, center[1]});
      Vector p3({center[0], center[1]+eps});
      double d1 = scalar_field(p1);
      double d2 = scalar_field(p2);
      double d3 = scalar_field(p3);

      double slope_x = (d2 - d1) / eps;
      double slope_y = (d3 - d1) / eps;
      Vector manual_gradient({slope_x, slope_y});

      cout << "DiscreteLinearOperator: ";
      Vector val_grad2;
      dx.GetVectorValue(elem_id, barycenter, val_grad2);
      val_grad2.Print();

      cout << "GradientGridFunctionCoefficient: ";
      Vector val_grad1;
      dx2.GetVectorValue(elem_id, barycenter, val_grad1);
      val_grad1.Print();

      cout << "x.GetGradient(): ";
      gradient.Print();

      cout << "Gradient calculated manually (finite difference): ";
      manual_gradient.Print();

      {
        // Q: how to actually output a vector solution?
        // A: NumperOfComponents is set automatically and we have nothing special to do.
        //
        // NOTE: paraview consider the data array as vector ONLY IF there are exactly
        // 3 components. For 2D data, then we need to use the Calculator filter with
        // grad_X*iHat+grad_Y*jHat
        ParaViewDataCollection dc("Gradient3D", &mesh);
        dc.SetPrefixPath("ParaView");
        dc.SetDataFormat(VTKFormat::BINARY);
        dc.RegisterField("x", &x);
        dc.RegisterField("dx", &dx);
        dc.RegisterField("dx2", &dx2);
        dc.Save();
      }

    }


    return 0;
}
