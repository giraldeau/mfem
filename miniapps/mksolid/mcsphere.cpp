#define MC_IMPLEM_ENABLE
#include "MC.h"
#include "spherecoeff.h"
#include "tetgen.h"

#include <fstream>
#include <iostream>

void CalculateNormal(const MC::mcVec3f &p1, const MC::mcVec3f &p2,
                     const MC::mcVec3f &p3, MC::mcVec3f &normal) {
  MC::mcVec3f u = p2 - p1;
  MC::mcVec3f v = p3 - p1;
  normal.x = (u.y * v.z) - (u.z * v.y);
  normal.y = (u.z * v.x) - (u.x * v.z);
  normal.z = (u.x * v.y) - (u.y * v.x);
}

int VectorsDirectionSimilar(const MC::mcVec3f &v1, const MC::mcVec3f &v2) {
  const MC::mcVec3f n1 = MC::mc_internalNormalize(v1);
  const MC::mcVec3f n2 = MC::mc_internalNormalize(v2);
  float dot = n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
  return (dot > 0.9); // close enough
}

int main() {

  MC::mcMesh mesh;

  {
    const int n = 50;
    float *field = new float[n * n * n];

    ImpSphereCoeff ls(mfem::Vector({0.5, 0.5, 0.5}), 0.25, 0.1);
    float dx = 1.0 / (n - 1);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
          field[(k * n + j) * n + i] =
              ls.Compute(Eigen::Vector3d(i * dx, j * dx, k * dx));
        }
      }
    }

    // Compute isosurface using marching cube
    MC::marching_cube(field, n, n, n, mesh);
  }

  {
    // Check normals
    int n = mesh.indices.size() / 3;
    for (size_t i = 0; i < n; i++) {
      int idx = i * 3;

      int i1 = mesh.indices[idx + 0];
      int i2 = mesh.indices[idx + 1];
      int i3 = mesh.indices[idx + 2];
      const MC::mcVec3f &p1 = mesh.vertices[i1];
      const MC::mcVec3f &p2 = mesh.vertices[i2];
      const MC::mcVec3f &p3 = mesh.vertices[i3];
      MC::mcVec3f n1, n2 = {0, 0, 0};
      CalculateNormal(p1, p2, p3, n1);

      n2 += mesh.normals[i1];
      n2 += mesh.normals[i2];
      n2 += mesh.normals[i3];
      n2.x /= 3;
      n2.y /= 3;
      n2.z /= 3;

      int res = VectorsDirectionSimilar(n1, n2);
      if (res == 0) {
        std::cout << "NORMAL CHECK FAILED: " << i << " " << res << std::endl;
      }
    }
  }

  {
    // Export the result as an .obj file
    std::ofstream out;
    out.open("test.obj");
    if (out.is_open() == false)
      return 1;
    out << "g Obj\n";
    for (size_t i = 0; i < mesh.vertices.size(); i++) {
      out << "v " << mesh.vertices.at(i).x << " " << mesh.vertices.at(i).y
          << " " << mesh.vertices.at(i).z << '\n';
    }
    for (size_t i = 0; i < mesh.vertices.size(); i++) {
      out << "vn " << mesh.normals.at(i).x << " " << mesh.normals.at(i).y << " "
          << mesh.normals.at(i).z << '\n';
    }
    for (size_t i = 0; i < mesh.indices.size(); i += 3) {
      out << "f " << mesh.indices.at(i) + 1 << "//" << mesh.indices.at(i) + 1
          << " " << mesh.indices.at(i + 1) + 1 << "//"
          << mesh.indices.at(i + 1) + 1 << " " << mesh.indices.at(i + 2) + 1
          << "//" << mesh.indices.at(i + 2) + 1 << '\n';
    }
    out.close();
  }

  {
    tetgenio in, out;

    in.numberofpoints = mesh.vertices.size();
    REAL *pointlist = new REAL[in.numberofpoints * 3];
    for (int i = 0; i < mesh.vertices.size(); i++) {
      const MC::mcVec3f &v = mesh.vertices[i];
      int idx = i * 3;
      pointlist[idx + 0] = v.x;
      pointlist[idx + 1] = v.y;
      pointlist[idx + 2] = v.z;
    }

    in.pointlist = pointlist;
    in.numberoffacets = mesh.indices.size() / 3;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = NULL;

    int nface = in.numberoffacets;
    for (size_t i = 0; i < nface; i++) {
      tetgenio::facet *f = &in.facetlist[i];
      tetgenio::init(f);
      f->numberofpolygons = 1;
      f->polygonlist = new tetgenio::polygon[1];
      tetgenio::polygon *p = &f->polygonlist[0];
      tetgenio::init(p);
      p->numberofvertices = 3;
      p->vertexlist = new int[3];

      int idx = i * 3;
      p->vertexlist[0] = mesh.indices[idx + 0];
      p->vertexlist[1] = mesh.indices[idx + 1];
      p->vertexlist[2] = mesh.indices[idx + 2];
    }

    tetrahedralize(const_cast<char *>("pV"), &in, &out);

    std::cout << "out numberofpoints :" << out.numberofpoints << "\n";
    std::cout << "out tets :" << out.numberoftetrahedra << "\n";

    mfem::Mesh mesh(3, out.numberofpoints, out.numberoftetrahedra);

    for (int i = 0; i < out.numberofpoints; i++) {
      mesh.AddVertex(&out.pointlist[i * 3]);
    }
    for (int i = 0; i < out.numberoftetrahedra; i++) {
      mesh.AddTet(&out.tetrahedronlist[i * 4]);
    }
    // mesh.Save("mcsphere.mesh");

    {
      mfem::ParaViewDataCollection paraview_dc("MCSphere", &mesh);
      paraview_dc.SetPrefixPath("ParaView");
      paraview_dc.SetLevelsOfDetail(2);
      paraview_dc.SetDataFormat(mfem::VTKFormat::BINARY);
      paraview_dc.Save();
    }
  }

  return 0;
}
