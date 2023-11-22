#include <iostream>
#include <map>

// Parallel Control Utility
#include <PCU.h>
#include <pcu_util.h> // asserts
 
// Parallel Fields/Meshes
#include <apfCAP.h>
#include <apfMDS.h>
#include <apfConvert.h>

// Geometric model interface
#include <gmi_cap.h>
#include <gmi_mesh.h>

// Capstone
#include "CapstoneModule.h"

using namespace CreateMG;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;

int main (int argc, char* argv[]) {
  // Check arguments
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " mesh.dmg model.smb output.cre" << std::endl;
    return -1;
  }

  // Initialize parallel communication
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  const char* model = argv[1], *mesh = argv[2], *creOut = argv[3];

  // Start gmi
  gmi_register_mesh();

  // Start gmi_cap
  gmi_cap_start();
  gmi_register_cap();

  // Load MDS mesh from file
  apf::Mesh2* m = apf::loadMdsMesh(model, mesh);
  m->verify();

  // Load CapstoneModule with SMLIB for geometry and CREATE for Mesh and Attribution
  CapstoneModule cs("mds2cre",
                    "Geometry Database : SMLIB",
                    "Mesh Database : Create",
                    "Attribution Database : Create");

  GeometryDatabaseInterface *gdi = cs.get_geometry();
  MeshDatabaseInterface *mdi = cs.get_mesh();

  PCU_ALWAYS_ASSERT(gdi);
  PCU_ALWAYS_ASSERT(mdi);

  // Create empty geometry model
  M_GModel gmodel;
  MG_API_CALL(gdi, create_model(gmodel));
  M_GBRep brep;
  MG_API_CALL(gdi, create_brep(brep));

  // Copy vertices
  std::map<apf::MeshEntity*,M_GVertex> v_map;
  apf::MeshIterator* it = m->begin(0);
  while (apf::MeshEntity* e = m->iterate(it)) {
    apf::Vector3 v;
    m->getPoint(e, 0, v);
    CreateMG::vec3d v2;
    v.toArray(v2.ptr());
    M_GVertex c_v;
    MG_API_CALL(gdi, create_vertex(brep, v2, c_v));
    auto result = v_map.insert({e, c_v});
    if (!result.second) {
      std::cerr << "Failed to update map." << std::endl;
      MPI_Finalize();
      return -1;
    } 
  }
  m->end(it);
  std::cout << "Copied " << v_map.size() << " vertices." << std::endl;

  // Copy edges
  std::map<apf::MeshEntity*, M_GEdge> e_map;
  it = m->begin(1);
  while (apf::MeshEntity* e = m->iterate(it)) {
    apf::Downward d;
    m->getDownward(e, 0, d);
    M_GEdge c_e;
    MG_API_CALL(gdi, create_edge(brep, v_map.at(d[0]), v_map.at(d[1]),  c_e));
    auto result = e_map.insert({e, c_e});
    if (!result.second) {
      std::cerr << "Failed to update edge map." << std::endl;
      MPI_Finalize();
      return -1;
    }
  }
  m->end(it);
  std::cout << "Copied " << e_map.size() << " edges." << std::endl;

  // Create mesh model
  M_MModel mmodel = cs.generate_mesh();
  MG_API_CALL(mdi, set_current_model(mmodel));

  // Set up mesh model
  MG_API_CALL(mdi, set_adjacency_state(REGION2FACE|REGION2EDGE|REGION2VERTEX|
                                       FACE2EDGE|FACE2VERTEX));
  MG_API_CALL(mdi, set_reverse_states());
  MG_API_CALL(mdi, compute_adjacency());

  // Write mesh to CRE.
  std::string fname(creOut);
  cs.save_file(fname, gmodel);

  // Destroy meshes.
  m->destroyNative();
  apf::destroyMesh(m);

  // Exit calls
  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}

