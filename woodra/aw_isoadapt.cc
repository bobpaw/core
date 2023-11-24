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

  double L = ma::getMaximumEdgeLength(m);

  class IF : public ma::IsotropicFunction {
    double val_;
    public:
      IF (double v) : val_(v) {}
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

