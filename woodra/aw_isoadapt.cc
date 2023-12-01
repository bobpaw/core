#include <iostream>

// Parallel Control Utility
#include <PCU.h>
#include <pcu_util.h>

// Parallel fields
#include <apf.h>
#include <apfMDS.h>

// Geometry model
#include <gmi_mesh.h>

// Mesh adapt
#include <ma.h>

int main (int argc, char* argv[]) {
  // Check arguments
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " inmesh.smb inmodel.dmg ";
    std::cerr << "outmesh.smb outmodel.dmg" << std::endl;
    return -1;
  }

  // Initialize parallelism
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  const char *inmesh = argv[1], *inmodel = argv[2],
             *outmesh = argv[3], *outmodel = argv[4];

  std::cout << "Loading mesh." << std::endl;
  gmi_register_mesh();
  ma::Mesh *m = apf::loadMdsMesh(inmodel, inmesh);
  m->verify();

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

  double L;
  ma::Vector length = maxs - mins;
  if (length.x() >= length.y() && length.x() >= length.z()) {
    // X is the longest axis.
    L = length.x();
  } else if (length.y() >= length.x() && length.y() >= length.z()) {
    // Y is the longest axis.
    L = length.y();
  } else {
    // Z is the longest axis.
    PCU_DEBUG_ASSERT(length.z() >= length.x() && length.z() >= length.y());
    L = length.z();
  }
  std::cout << "Adapting with L = " << L << std::endl;

  class IF : public ma::IsotropicFunction {
    double val_;
    public:
      IF (double v) : val_(v) {
        std::cout << "Creating IsoFunc with v = " << v << std::endl;
      }
      virtual double getValue(ma::Entity*) {
        return val_;
      }
      virtual ~IF() {}
  } func(L/20);

  std::cout << "Adapting mesh." << std::endl;
  ma::adapt(m, &func);
  m->verify();

  std::cout << "Writing mesh." << std::endl;
  m->writeNative(outmesh);

  std::cout << "Writing model." << std::endl;
  gmi_write_dmg(m->getModel(), outmodel);

  // Destroy objects
  m->destroyNative();
  apf::destroyMesh(m);

  // Exit calls
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

