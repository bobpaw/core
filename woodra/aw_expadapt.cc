#include <iostream>

// Parallel Control Utility
#include <PCU.h>
#include <pcu_util.h>

// Attached parallel fields
#include <apf.h>
#include <apfMDS.h>

// Geometry model
#include <gmi_mesh.h>

// Mesh adapt
#include <ma.h>

class ExpLengthFunc : public ma::IsotropicFunction {
  enum class Axis { X, Y, Z } axis;
  double midpoint, avg_edge_len;
  ma::Mesh *m;

  double ldist (const ma::Vector &v) const {
    switch (axis) {
      case Axis::X:
        return std::fabs(v.x() - midpoint);
      case Axis::Y:
        return std::fabs(v.y() - midpoint);
      case Axis::Z:
        return std::fabs(v.z() - midpoint);
    }
  }

  public:
  ExpLengthFunc(ma::Mesh *mesh) : m(mesh) {
    // Find mesh bounds.
    ma::Vector mins(HUGE_VAL, HUGE_VAL, HUGE_VAL), maxs(-HUGE_VAL, -HUGE_VAL, -HUGE_VAL);
    ma::Iterator *it = m->begin(0);
    while (ma::Entity *ent = m->iterate(it)) {
      ma::Vector v = ma::getPosition(m, ent);

      if (v.x() < mins.x()) mins.x() = v.x();
      if (v.x() > maxs.x()) maxs.x() = v.x();
      if (v.y() < mins.y()) mins.y() = v.y();
      if (v.y() > maxs.y()) maxs.y() = v.y();
      if (v.z() < mins.z()) mins.z() = v.z();
      if (v.z() > maxs.z()) maxs.z() = v.z();
    }
    m->end(it);

    ma::Vector length = maxs - mins;
    if (length.x() >= length.y() && length.x() >= length.z()) {
      // X is the longest axis.
      axis = Axis::X;
      midpoint = (mins.x() + maxs.x()) / 2.0;
    } else if (length.y() >= length.x() && length.y() >= length.z()) {
      // Y is the longest axis.
      axis = Axis::Y;
      midpoint = (mins.y() + maxs.x()) / 2.0;
    } else {
      // Z is the longest axis.
      PCU_DEBUG_ASSERT(length.z() >= length.x() && length.z() >= length.y());
      axis = Axis::Z;
      midpoint = (mins.z() + maxs.z()) / 2.0;
    }

    avg_edge_len = ma::getAverageEdgeLength(m);
  }

  ~ExpLengthFunc () {}

  double getValue(ma::Entity *vert) {
    ma::Vector v = ma::getPosition(m, vert);
    double dist = ldist(v);
    double h_min = avg_edge_len / 4.0, h_max = avg_edge_len * 2.0;
    double a = 0.1542;
    return h_max + (h_max - h_min) * std::exp(-a*ldist(v)/avg_edge_len);
  }
};

int main (int argc, char *argv[]) {
  // Check arguments
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " mesh.smb model.dmg ";
    std::cerr << "outmesh.smb outmodel.dmg" << std::endl;
    return -1;
  }

  // Initialize parallelism
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  const char *imesh = argv[1], *imodel = argv[2], *omesh = argv[3], *omodel = argv[4];

  // Loading mesh
  std::cout << "Loading mesh and model." << std::endl;
  gmi_register_mesh();
  ma::Mesh *m = apf::loadMdsMesh(imodel, imesh);
  m->verify();

  // Adapt mesh
  std::cout << "Creating size function." << std::endl;
  ExpLengthFunc elf(m);
  std::cout << "Adapting mesh." << std::endl;
  ma::adapt(m, &elf);
  m->verify();
  
  // Write mesh
  std::cout << "Writing mesh." << std::endl;
  m->writeNative(omesh);
  std::cout << "Writing model." << std::endl;
  gmi_write_dmg(m->getModel(), omodel);

  // Destroy objects
  m->destroyNative();
  apf::destroyMesh(m);

  // Exit calls
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

