#include <iostream>

#include "mfem.hpp"
#include <Eigen/Core>

using namespace mfem;

void project_to_unit_sphere(Mesh &mesh) {
  for (int i = 0; i < mesh.GetNV(); i++) {
    real_t *v = mesh.GetVertex(i);
    Eigen::Vector3d p(v);
    real_t n = p.norm();
    p = (1.0 / n) * p;
    for (int j = 0; j < 3; j++) {
      v[j] = p[j];
    }
  }
}

Mesh icosahedron() {

  Mesh mesh(2, 12, 20, 0, 3);

  real_t phi = (1.0f + sqrt(5.0f)) * 0.5f; // golden ratio
  real_t a = 1.0f;
  real_t b = 1.0f / phi;
  real_t z = 0.0;

  // add vertices
  auto v1 = mesh.AddVertex(Vector({z, b, -a}));
  auto v2 = mesh.AddVertex(Vector({b, a, z}));
  auto v3 = mesh.AddVertex(Vector({-b, a, z}));
  auto v4 = mesh.AddVertex(Vector({z, b, a}));
  auto v5 = mesh.AddVertex(Vector({z, -b, a}));
  auto v6 = mesh.AddVertex(Vector({-a, z, b}));
  auto v7 = mesh.AddVertex(Vector({z, -b, -a}));
  auto v8 = mesh.AddVertex(Vector({a, z, -b}));
  auto v9 = mesh.AddVertex(Vector({a, z, b}));
  auto v10 = mesh.AddVertex(Vector({-a, z, -b}));
  auto v11 = mesh.AddVertex(Vector({b, -a, z}));
  auto v12 = mesh.AddVertex(Vector({-b, -a, z}));

  project_to_unit_sphere(mesh);

  // add triangles
  mesh.AddTriangle(v3, v2, v1);
  mesh.AddTriangle(v2, v3, v4);
  mesh.AddTriangle(v6, v5, v4);
  mesh.AddTriangle(v5, v9, v4);
  mesh.AddTriangle(v8, v7, v1);
  mesh.AddTriangle(v7, v10, v1);
  mesh.AddTriangle(v12, v11, v5);
  mesh.AddTriangle(v11, v12, v7);
  mesh.AddTriangle(v10, v6, v3);
  mesh.AddTriangle(v6, v10, v12);
  mesh.AddTriangle(v9, v8, v2);
  mesh.AddTriangle(v8, v9, v11);
  mesh.AddTriangle(v3, v6, v4);
  mesh.AddTriangle(v9, v2, v4);
  mesh.AddTriangle(v10, v3, v1);
  mesh.AddTriangle(v2, v8, v1);
  mesh.AddTriangle(v12, v10, v7);
  mesh.AddTriangle(v8, v11, v7);
  mesh.AddTriangle(v6, v12, v5);
  mesh.AddTriangle(v11, v9, v5);

  return mesh;
}

Mesh icosphere(int refinement) {
  Mesh mesh = icosahedron();

  for (int i = 0; i < refinement; i++) {
    mesh.UniformRefinement();
    project_to_unit_sphere(mesh);
  }
  return mesh;
}

Mesh plane(int refinement) {
  Mesh mesh(2, 4, 1, 0, 3);
  auto v1 = mesh.AddVertex(Vector({0, 0, 0}));
  auto v2 = mesh.AddVertex(Vector({0, 1, 0}));
  auto v3 = mesh.AddVertex(Vector({1, 1, 0}));
  auto v4 = mesh.AddVertex(Vector({1, 0, 0}));
  mesh.AddTriangle(v1, v2, v3);
  mesh.AddTriangle(v1, v3, v4);

  for (int i = 0; i < refinement; i++) {
    mesh.UniformRefinement();
  }

  return mesh;
}

Mesh wave(int refinement) {
  Mesh mesh = plane(refinement);

  real_t w = 2 * M_PI;

  for (int i = 0; i < mesh.GetNV(); i++) {
    real_t *v = mesh.GetVertex(i);
    // z = cos(x) + sin(y)
    v[0] *= w;
    v[1] *= w;
    v[2] = std::cos(v[0]);
  }
  return mesh;
}

Mesh freq(int refinement) {


  int ny = 1;
  int nx = 10;
  double end = 10.0 * 2 * M_PI;
  Mesh base(Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL, false, end, 10.0));

  for (int i = 0; i < refinement; i++) {
    base.UniformRefinement();
  }

  Mesh mesh(2, base.GetNV(), base.GetNE(), 0, 3);
  real_t node[3];
  for (int i = 0; i < base.GetNV(); i++) {
    real_t *v = base.GetVertex(i);
    node[0] = v[0];
    node[1] = v[1];
    node[2] = std::sin(v[0] * v[0] * 0.02);
    mesh.AddVertex(node);
  }

  for (int i = 0; i < base.GetNE(); i++) {
    Element *element = base.GetElement(i);
    int *vertices = element->GetVertices();
    mesh.AddElement(new Quadrilateral(vertices));
  }

  return mesh;
}

int main(int argc, char **argv) {

  int refinement = 2;
  int geom = 0;
  std::string outname = "surf";
  int simplex = 0;

  OptionsParser args(argc, argv);
  args.AddOption(&refinement, "-r", "--refinement", "Refinement steps");
  args.AddOption(&geom, "-g", "--geometry",
                 "Geometry: icosphere = 0, plate = 1, wave = 2, freq = 3");
  args.AddOption(&simplex, "-s", "--simplex", "Transform mesh to simplical");
  args.AddOption(&outname, "-o", "--output", "Output file basename");

  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return 1;
  }
  args.PrintOptions(std::cout);

  Mesh mesh;
  if (geom == 0) {
    // generate icosphere
    mesh = icosphere(refinement);
  } else if (geom == 1) {
    // generate plane
    mesh = plane(refinement);
  } else if (geom == 2) {
    // generate wave
    mesh = wave(refinement);
  } else if (geom == 3) {
    // generate plane with varying frequency
    mesh = freq(refinement);
  } else {
    std::cout << "unkown geometry\n";
  }

  if (simplex && mesh.HasGeometry(Geometry::SQUARE)) {
    mesh = Mesh::MakeSimplicial(mesh);
  }

  mesh.PrintInfo();
  mesh.Save(outname + ".mesh");

  std::ofstream ofs(outname + ".vtk");
  mesh.PrintVTK(ofs);

  return 0;
}
