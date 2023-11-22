#include <iostream>

// Parallel Control Utility
#include <PCU.h>
#include <pcu_util.h>

// Parallel fields
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>

// Geometry model
#include <gmi_mesh.h>

int main (int argc, char* argv[]) {
  // Check arguments
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " mesh.smb model.dmg vtkPrefix" << std::endl;
    return -1;
  }

  // Initialize parallelism
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  const char *mesh = argv[1], *model = argv[2], *vtk = argv[3];

  std::cout << "Loading mesh." << std::endl;
  gmi_register_mesh();
  apf::Mesh2 *m = apf::loadMdsMesh(model, mesh);
  m->verify();

  std::cout << "Writing VTKs." << std::endl;
  apf::writeVtkFiles(vtk, m);

  // Destroy objects
  m->destroyNative();
  apf::destroyMesh(m);

  // Exit calls
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

