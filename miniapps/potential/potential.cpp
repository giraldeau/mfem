#include <mfem.hpp>

using namespace mfem;

const double k = 1.0;

void pointcloud_save(const std::string& filename, std::vector<Vector> &points) {
  std::ofstream out(filename);
  if (!out.is_open()) {
    std::cerr << "Error: Could not open the file for writing.\n";
    return;
  }

  // Writing the VTK header
  out << "# vtk DataFile Version 3.0\n";
  out << "Point cloud data\n";
  out << "ASCII\n";
  out << "DATASET POLYDATA\n";

  // Writing the points
  out << "POINTS " << points.size() << " float\n";
  for (const auto& point : points) {
    out << point[0] << " " << point[1] << " " << point[2] << "\n";
  }

  // Writing the vertices
  int nv = points.size() * 2;
  out << "VERTICES " << points.size() << " " << nv << "\n";
  for (int i = 0; i < points.size(); i++) {
    out << 1 << " " << i << "\n";
  }

  out.close();
}

struct Particle {

  // Constructeur: position, vitesse et charge
  Particle(const Vector &p_, double q_)
      : m_x(p_), m_q(q_) {}

  Vector potential_at(const Vector &loc) const {
    double r = loc.DistanceTo(m_x) + 1e-6;
    Vector u = loc;
    u -= m_x;
    u *= k * m_q / (r*r);
    return u;
  };

  Vector m_x; // position courante
  double m_q;   // charge
};

static int calls = 0;
std::vector<Vector> gauss_points;

std::function<void(const Vector &, Vector &)> f_vec(std::vector<Particle> &lst)
{
  return [&](const Vector &xvec, Vector &f)
  {
    calls++;
    gauss_points.push_back(xvec);
    f.SetSize(3);
    f = 0.0;
    for (int i = 0; i < lst.size(); i++) {
      const Particle &p = lst.at(i);
      f += p.potential_at(xvec);
    }
    std::cout << calls << std::endl;
    xvec.Print();
  };
}

int main()
{
  Particle p1(Vector({0.4, 0.5, 0.5}), 1.0);
  Particle p2(Vector({0.6, 0.5, 0.5}), -1.0);
  std::vector<Particle> lst = {p1, p2};

  Mesh mesh(Mesh::MakeCartesian3D(1, 1, 1, Element::QUADRILATERAL));
  H1_FECollection fec(3, mesh.Dimension());
  FiniteElementSpace fespace(&mesh, &fec, mesh.Dimension());
  VectorFunctionCoefficient field(mesh.Dimension(), f_vec(lst));
  GridFunction gf(&fespace);
  gf.ProjectCoefficient(field);

  {
    ParaViewDataCollection pvd("Particle", &mesh);
    pvd.SetLevelsOfDetail(3);
    pvd.RegisterField("e", &gf);
    pvd.SetDataFormat(VTKFormat::ASCII);
    pvd.Save();
  }

  mesh.Save("particle.mesh");
  gf.Save("particle.gf");

  std::cout << "calls: " << calls << std::endl;

  pointcloud_save("gauss.vtk", gauss_points);

  return 0;
}
