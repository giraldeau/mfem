
#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;


int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../data/inline-quad.mesh";
   bool vis = true;
   OptionsParser args(argc, argv);
   args.AddOption(&vis, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   Mesh mesh(mesh_file, 1, 1);

   Array<int> ref_elems(1);
   ref_elems = 2;
   mesh.GeneralRefinement(ref_elems,-1,0);
   ref_elems = 5;
   mesh.GeneralRefinement(ref_elems,-1,0);

   <<<<<<< HEAD
   Array<int> faces({0,1,10,26,28,29});
   =======
      Array<int> faces({0,10,26,28,29});
   >>>>>>> zz_test

   for (int i = 0; i<faces.Size(); i++)
   {
      Array<int> elems;
      <<<<<<< HEAD
      int type = mesh.GetFaceElements2(faces[i],elems);
      cout << "type = " << type << endl;
      cout << "elems = " ; elems.Print();
      =======
         mesh.GetFaceElements(faces[i],elems);

      >>>>>>> zz_test
      if (vis)
      {
         char vishost[] = "localhost";
         int  visport   = 19916;
         socketstream patch_sock(vishost, visport);
         L2_FECollection fec(1,mesh.Dimension());
         FiniteElementSpace fes(&mesh,&fec);
         GridFunction vis_gf(&fes);
         vis_gf = 0.0;
         Array<int> dofs;
         for (int i=0; i<elems.Size(); i++)
         {
            int el = elems[i];
            fes.GetElementDofs(el, dofs);
            vis_gf.SetSubVector(dofs,1.0);
         }
         patch_sock.precision(8);
         patch_sock << "solution\n" << mesh << vis_gf <<
                    "keys rRmjnppppp\n" << flush;
      }
   }

   return 0;
}