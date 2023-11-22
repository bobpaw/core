#include <iostream>

// Parallel Control Utility
#include <PCU.h>
#include <pcu_util.h>

// Parallel Fields/Meshes
#include <apf.h>
#include <apfMDS.h>

// Geometric model interface
#include <gmi_mesh.h>

// Mesh adapt
#include <ma.h>

int main (int argc, char* argv[]) {
  // Check arguments
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " mesh.smb model.dmg vtkPrefix" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Initalize parallel communication
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  const char *meshFile = argv[1], *modelFile = argv[2], *vtkPrefix = argv[3];

  // Register GMI for .dmg mesh file.
  gmi_register_mesh();

  ma::Mesh *m = apf::loadMdsMesh(modelFile, meshFile);
  m->verify();
  const ma::Input *in = ma::configureUniformRefine(m);
  ma::adapt(in);
  m->verify();

  // Write output
  apf::writeVtkFiles(vtkPrefix, m);

  // Destroy objects
  m->destroyNative();
  apf::destroyMesh(m);

  // Exit calls
  PCU_Comm_Free();
  MPI_Finalize();
}

