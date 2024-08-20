#include "surfmeshcompat.h"

#include <pmp/algorithms/curvature.h>
#include <pmp/surface_mesh.h>
#include <pmp/io/io.h>

SurfMeshCompat::SurfMeshCompat() : p_imp(new pmp::SurfaceMesh()) {}

SurfMeshCompat::~SurfMeshCompat() {}

void SurfMeshCompat::reserve(int nv, int nedge, int ne) {
  p_imp->reserve(nv, nedge, ne);
}

int SurfMeshCompat::add_vertex(double x, double y, double z) {
  pmp::Point point;
  point[0] = x;
  point[1] = y;
  point[2] = z;
  pmp::Vertex v = p_imp->add_vertex(point);
  return v.idx();
}

int SurfMeshCompat::add_face(int *index, int n) {
  std::vector<pmp::Vertex> vertices(n);
  for (int i = 0; i < n; i++) {
    vertices[i] = pmp::Vertex(index[i]);
  }
  pmp::Face f = p_imp->add_face(vertices);
  return f.idx();
}

void SurfMeshCompat::get_curvature(std::vector<double> &vec, Curvature c,
                                   int smoothing_step, bool use_tensor,
                                   bool use_two_ring) {
  pmp::curvature(*p_imp, static_cast<pmp::Curvature>(c), smoothing_step,
                 use_tensor, use_two_ring);
  pmp::VertexProperty<pmp::Scalar> vp =
      p_imp->get_vertex_property<pmp::Scalar>("v:curv");
  vec.resize(p_imp->vertices_size());
  for (auto v : p_imp->vertices()) {
    pmp::Vertex vertex = v;
    pmp::Scalar x = vp[v];
    vec[v.idx()] = x;
  }
}

void SurfMeshCompat::write(const std::string &filename)
{
  pmp::write(*p_imp, filename);
}

