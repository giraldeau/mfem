#include <iostream>

#include "mfem.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace mfem;

int main(int argc, char **argv) {

  Eigen::IOFormat dbgfmt(3, Eigen::DontAlignCols);

  // Load triangular surface mesh in 3D
  Mesh surf(argv[1]);

  MFEM_ASSERT(surf.SpaceDimension() != 1, "unsupported for 1D mesh");

  // Compute normal at each vertex
  // Compute normal per element, then average in vnorm
  std::vector<Eigen::Vector3d> vnorm(surf.GetNV());
  for (int i = 0; i < surf.GetNE(); i++) {
    Element *e = surf.GetElement(i);
    int nv = e->GetNVertices();
    int *v = e->GetVertices();

    Eigen::Vector3d p1, p2, p3;

    real_t *t1 = surf.GetVertex(v[0]);
    real_t *t2 = surf.GetVertex(v[1]);
    real_t *t3 = surf.GetVertex(v[2]);

    for (int d = 0; d < surf.SpaceDimension(); d++) {
      p1[d] = t1[d];
      p2[d] = t2[d];
      p3[d] = t3[d];
    }

    Eigen::Vector3d v1 = p2 - p1;
    Eigen::Vector3d v2 = p3 - p1;

    Eigen::Vector3d normal = v1.cross(v2).normalized();
    for (int j = 0; j < nv; j++) {
      int id = v[j];
      vnorm[id] += normal;
    }
  }

  // Establish the extrusion thickness vector
  // Note: at this time, we assume a uniform unit thickness
  double thick = 1.0;
  for (int i = 0; i < surf.GetNV(); i++) {
    vnorm[i] = vnorm[i].normalized() * thick;
    std::cout << "norm " << i << " " << vnorm[i].transpose().format(dbgfmt)
              << "\n";
  }

  Eigen::Vector3d up(0, 0, 1);

  // Create volume mesh by extruding each surface element
  int nz = 10;
  int nvz = nz + 1;
  int nvt = surf.GetNV() * nvz;
  Mesh vol(3, nvt, surf.GetNE() * nz);

  // Create vertices and elements. This code is a modified
  // version of Extrude2D function, but extruding in the
  // direction of the vertex normal.
  real_t vcoords[3];
  for (int i = 0; i < surf.GetNV(); i++) {
    vcoords[0] = surf.GetVertex(i)[0];
    vcoords[1] = surf.GetVertex(i)[1];
    vcoords[2] = 0.0;
    if (surf.SpaceDimension() == 3) {
      vcoords[2] = surf.GetVertex(i)[2];
    }

    if (i == 428) {
      std::cout << "DEBUG NODE\n";
    }

    Eigen::Vector3d p0(vcoords);

    // if we extrude in the z direction, the thickness is incorrect
    // adjust the extrusion length to preserve the thickness normal
    // to the surface
    // Note: does not work, because we project the point too far away
    // Eigen::Vector3d up(0, 0, 1);
    // Eigen::Vector3d normal = vnorm[i];
    // double thick = 1.0;
    // double t = thick / normal.dot(up);
    // Eigen::Vector3d extrude_vec = t * up;
    // Eigen::Vector3d offset = extrude_vec / nz;

    Eigen::Vector3d offset = vnorm[i] / nz;
    for (int j = 0; j < nvz; j++) {
      Eigen::Vector3d p1 = p0 + offset * j;
      vol.AddVertex(p1.data());
    }
  }

  // Create elements
  Array<int> vert;
  for (int i = 0; i < surf.GetNE(); i++) {
    const Element *elem = surf.GetElement(i);
    elem->GetVertices(vert);
    const int attr = elem->GetAttribute();
    Geometry::Type geom = elem->GetGeometryType();
    if (geom == Geometry::TRIANGLE) {
      for (int j = 0; j < nz; j++)
      {
        int pv[6];
        pv[0] = vert[0] * nvz + j;
        pv[1] = vert[1] * nvz + j;
        pv[2] = vert[2] * nvz + j;
        pv[3] = vert[0] * nvz + (j + 1) % nvz;
        pv[4] = vert[1] * nvz + (j + 1) % nvz;
        pv[5] = vert[2] * nvz + (j + 1) % nvz;

        vol.AddWedge(pv, attr);
      }
    } else if (geom == Geometry::SQUARE) {
      for (int j = 0; j < nz; j++)
      {
        int hv[8];
        hv[0] = vert[0] * nvz + j;
        hv[1] = vert[1] * nvz + j;
        hv[2] = vert[2] * nvz + j;
        hv[3] = vert[3] * nvz + j;
        hv[4] = vert[0] * nvz + (j + 1) % nvz;
        hv[5] = vert[1] * nvz + (j + 1) % nvz;
        hv[6] = vert[2] * nvz + (j + 1) % nvz;
        hv[7] = vert[3] * nvz + (j + 1) % nvz;

        vol.AddHex(hv, attr);
      }
    }

  }

  vol.FinalizeMesh(0, false);

  // Save output mesh
  vol.Save("triwedge.mesh");
  std::ofstream ofs("triwedge.vtk");
  vol.PrintVTK(ofs);
  return 0;
}
