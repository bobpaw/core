#include <iostream>

// Parallel Control Utility
#include <PCU.h>
#include <pcu_util.h> // asserts

// Parallel Fields/Meshes
#include <apf.h>
#include <apfCAP.h>
#include <apfConvert.h>
#include <apfMDS.h>

// Geometric model interface
#include <gmi_cap.h>
#include <gmi_null.h>

// Capstone
#include "CapstoneModule.h"

using namespace CreateMG;

int main (int argc, char* argv[]) {
  // Check arguments
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " mesh.smb model.dmg mesh.cre..." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Initalize parallel communication
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  const char* meshFile = argv[1], *modelFile = argv[2];

  // Start gmi_cap
  gmi_cap_start();
  gmi_register_cap();

  CapstoneModule cs("cre2vtk",
                    "Geometry Database : SMLIB",
                    "Mesh Database : Create",
                    "Attribution Database : Create");

  
  Geometry::GeometryDatabaseInterface *gdi = cs.get_geometry();
  Mesh::MeshDatabaseInterface *mdi = cs.get_mesh();

  PCU_ALWAYS_ASSERT(gdi);
  PCU_ALWAYS_ASSERT(mdi);

  v_string filenames;
  for (int i = 3; i < argc; ++i) filenames.push_back(argv[i]);

  M_GModel gmodel = cs.load_files(filenames);

  // Pick the volume mesh model from associated mesh models to this geom model
  std::vector<M_MModel> mmodels;
  MG_API_CALL(mdi, get_associated_mesh_models(gmodel, mmodels));
  std::cout << "Mesh model count: " << mmodels.size() << std::endl;
  if (mmodels.size() != 0) {
    PCU_ALWAYS_ASSERT(mmodels.size() == 1);
    MG_API_CALL(mdi, set_current_model(mmodels[0]));
  } else {
    // We must reconstruct the mesh model from the geometry.
    M_MModel mmodel = cs.generate_mesh();
    if (mmodel.is_invalid()) {
      std::cerr << "Failed to mesh the model." << std::endl;
      return -1;
    }
    MG_API_CALL(mdi, set_current_model(mmodel));
  }
  /* SET THE ADJACENCIES */
  MG_API_CALL(mdi, set_adjacency_state(REGION2FACE|
				     REGION2EDGE|
				     REGION2VERTEX|
				     FACE2EDGE|
				     FACE2VERTEX));
  MG_API_CALL(mdi, set_reverse_states());
  MG_API_CALL(mdi, compute_adjacency());

  size_t num_v, num_e, num_f, num_r; 
  MG_API_CALL(mdi, get_num_topos(Mesh::TOPO_VERTEX, num_v));
  MG_API_CALL(mdi, get_num_topos(Mesh::TOPO_EDGE, num_e));
  MG_API_CALL(mdi, get_num_topos(Mesh::TOPO_FACE, num_f));
  MG_API_CALL(mdi, get_num_topos(Mesh::TOPO_REGION, num_r));
  std::cout << "Mesh model vertices: " << num_v << std::endl;
  std::cout << "Mesh model edges: " << num_e << std::endl;
  std::cout << "Mesh model faces: " << num_f << std::endl;
  std::cout << "Mesh model regions: " << num_r << std::endl;

  /* CONVERT THE MESH TO APF::MESH2 */
  apf::Mesh2* m = apf::createMesh(mdi, gdi);
  std::cout << "Verifying." << std::endl;
  m->verify();

  // Convert to MDS.
  int dim;
  MG_API_CALL(mdi, get_dimension(dim));

  gmi_register_null();
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m2 = apf::makeEmptyMdsMesh(model, dim, true);
  apf::convert(m, m2);
  m2->verify();
  std::cout << "Writing mesh." << std::endl;
  m2->writeNative(meshFile);
  gmi_write_dmg(model, modelFile);

  // Destroy objects
  apf::destroyMesh(m);

  // Exit calls
  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}

