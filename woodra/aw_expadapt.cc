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

#include "aw_explf.h"

class IsoExpLengthFunc : public ma::IsotropicFunction, public ExpLengthFunc {

  public:
  IsoExpLengthFunc (ma::Mesh *mesh) : ExpLengthFunc(mesh) {}

  ~IsoExpLengthFunc () {}

  double getValue(ma::Entity *vert) {
    ma::Vector v = ma::getPosition(m, vert);
    return eval(v);
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
  IsoExpLengthFunc ielf(m);
  std::cout << "Adapting mesh." << std::endl;
  ma::adapt(m, &ielf);
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

