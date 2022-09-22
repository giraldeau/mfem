//                                MFEM Tutorial 1
//
// Compile with: make t1
//
// Sample runs: ./t1
//
// Description: This tutorial is about the operations at the core of the finite
// element method, namely transformations to and from the reference element and
// numerical integration. We focus on a 2d triangular element to make
// visualization easier.

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct Point2 {
   double x;
   double y;
};

// We will integrate this field over the element.
double func(const Vector &v)
{
   MFEM_ASSERT(v.Size() == 2, "unsupported dim");
   double val = 0.0;
   return val;
}

Mesh *make_mesh()
{
   Mesh * mesh = new Mesh(2, 3, 1, 3);

   // vertices
   double vc[2];
   vc[0] = 3.0; vc[1] = 4.0;
   mesh->AddVertex(vc);
   vc[0] = 8.0; vc[1] = 3.0;
   mesh->AddVertex(vc);
   vc[0] = 6.0; vc[1] = 7.0;
   mesh->AddVertex(vc);

   // element
   Array<int> vert(3);
   vert[0] = 0; vert[1] = 1; vert[2] = 2;
   mesh->AddTri(vert, 1);

   // boundary
   Array<int> sv(2);
   sv[0] = 0; sv[1] = 1;
   mesh->AddBdrSegment(sv, 1);
   sv[0] = 1; sv[1] = 2;
   mesh->AddBdrSegment(sv, 2);
   sv[0] = 2; sv[1] = 0;
   mesh->AddBdrSegment(sv, 3);

   mesh->FinalizeTriMesh(1, 0, true);

   return mesh;
}

void write_points(ostream &f, const vector<Point2> &v)
{
   for (const auto &p : v) {
      f << p.x << " " << p.y << "\n";
   }
}

void write_tri(ostream &f, const vector<Point2> &v)
{
   int n = 0;
   for (const auto &p : v) {
      f << p.x << " " << p.y << "\n";
      // insert blank line after 3 points
      if (++n % 3 == 0) {
         f << "\n";
      }
   }
}

static const char *gnuplot_commands = R"(
set term png
set output 't1.png'
set xlabel 'X'
set ylabel 'Y'
set title 'transforms'
set grid
set key off
set size ratio -1
set style data lines
plot 't1_tri.txt' using 1:2 with filledcurves closed, \
     't1_ref.txt' with points pt 7 lt 0, \
     't1_trans.txt' with points pt 7 lt 0
)";

int main(int argc, char *argv[])
{
   // 1. Parse command line options.
   int order = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree.");
   args.ParseCheck();

   // 2. Create the mesh. The mesh contains one 2d triangle with some arbitrary
   // vertices.
   std::unique_ptr<Mesh> mesh(make_mesh());

   // 3. Show mesh characteristics and save the mesh
   mesh->PrintInfo();
   mesh->Save("t1.mesh");

   // 4. Transform points in the reference element to the deformed element.
   vector<Point2> ref_pts = {
      { 0.2, 0.6 },
      { 0.2, 0.4 },
      { 0.2, 0.2 },
      { 0.4, 0.4 },
      { 0.4, 0.2 },
      { 0.6, 0.2 }
   };
   vector<Point2> trans_pts(ref_pts.size());
   IsoparametricTransformation Tr;
   mesh->GetElementTransformation(0, &Tr);

   for (size_t i = 0; i < ref_pts.size(); i++) {
      const Point2 &p = ref_pts[i];
      IntegrationPoint ip;
      Vector trans;
      ip.Set2(p.x, p.y);
      // transform the point from the reference element to physical element
      Tr.Transform(ip, trans);
      trans_pts[i] = { trans[0], trans[1] };
   }

   {
      vector<Point2> tri_coords = {
         { 0, 0 },
         { 0, 1 },
         { 1, 0 },

         { 3, 4 },
         { 8, 3 },
         { 6, 7 },
      };
      ofstream out_pref("t1_ref.txt");
      ofstream out_ptrans("t1_trans.txt");
      ofstream out_tri("t1_tri.txt");
      write_points(out_pref, ref_pts);
      write_points(out_ptrans, trans_pts);
      write_tri(out_tri, tri_coords);
      ofstream out_plot("t1.plt");
      out_plot << gnuplot_commands;
   }

   return 0;
}
