#include <iostream>

// Parallel Control
#include <PCU.h>
#include <pcu_util.h> // pcu assert

// Attached parallel fields
#include <apf.h>
#include <apfMDS.h>

// Geometry model
#include <gmi_mesh.h>

// Mesh adapt
#include <ma.h>

#include "aw_explf.h"

class HessianFunc {
  public:
  virtual void getHessian (const ma::Vector &v, ma::Matrix &H) = 0;
};

class DummyHessian : public HessianFunc {
  void getHessian (const ma::Vector &v, ma::Matrix &H) {
    H =  ma::Matrix(2, 0, 0,
                    0, 2, 0,
                    0, 0, 2);
  }
};

class SpeedSolution : public ma::AnisotropicFunction, public ExpLengthFunc {
  HessianFunc *hf;
  public:
  SpeedSolution(ma::Mesh *m, HessianFunc *h) : ExpLengthFunc(m), hf(h) {}
  ~SpeedSolution () {}

  void getValue(ma::Entity *vert, ma::Matrix &r, ma::Vector &h) {
    ma::Vector v = ma::getPosition(m, vert);
    ma::Matrix H;
    hf->getHessian(v, H);
    ma::Vector evecs[3];
    double evals[3];
    eigen(H, evecs, evals);
    // Copy eigen-values and eigen-vectors into adapt frame and scale
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) r[i][j] = evecs[i][j];
      h[i] = evals[i];
    }
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
  DummyHessian d;
  SpeedSolution sol(m, &d);
  std::cout << "Adapting mesh." << std::endl;
  ma::adapt(m, &sol);
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

