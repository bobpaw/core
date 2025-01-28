#include <PCU.h>
#include <apfCAP.h>
#include <apfMDS.h>
#include <ma.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <apf.h>
#include <apfMesh2.h>
#include <ma.h>
#include <pcu_util.h>
#include <iostream>
#include <vector>


#include "CapstoneModule.h"

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;


int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  if (argc != 3) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0]
        << " <input.cre> <output.mds>\n";
    return EXIT_FAILURE;
  }

  const char* creFileName = argv[1];
  const char* mdsFileName = argv[2];

  const std::string gdbName("Geometry Database : SMLIB");// Switch Create with SMLIB for CAD
  const std::string mdbName("Mesh Database : Create");
  const std::string adbName("Attribution Database : Create");

  CapstoneModule  cs("cap2mds", gdbName.c_str(), mdbName.c_str(), adbName.c_str());

  GeometryDatabaseInterface     *g = cs.get_geometry();
  MeshDatabaseInterface         *m = cs.get_mesh();

  v_string filenames;
  filenames.push_back(creFileName);

  M_GModel gmodel = cs.load_files(filenames);

  int numbreps = 0;
  MG_CALL(g->get_num_breps(numbreps));
  std::cout << "number of b reps is " << numbreps << std::endl;
  if(numbreps == 0)
      error(HERE, ERR_INVALID_INPUT, "Model is empty");

  M_MModel mmodel;
  // Pick the volume mesh model from associated mesh models to this geom model
  std::vector<M_MModel> mmodels;
  MG_API_CALL(m, get_associated_mesh_models(gmodel, mmodels));
  for(std::size_t i = 0; i < mmodels.size(); ++i) {
    M_MModel ammodel = mmodels[i];
    std::size_t numregs = 0;
    std::size_t numfaces = 0;
    std::size_t numedges = 0;
    std::size_t numverts = 0;
    MG_API_CALL(m, set_current_model(ammodel));
    MG_API_CALL(m, get_num_topos(TOPO_REGION, numregs));
    MG_API_CALL(m, get_num_topos(TOPO_FACE, numfaces));
    MG_API_CALL(m, get_num_topos(TOPO_EDGE, numedges));
    MG_API_CALL(m, get_num_topos(TOPO_VERTEX, numverts));
    std::cout << "num regions is " << numregs << std::endl;
    std::cout << "num faces   is " << numfaces << std::endl;
    std::cout << "num edges   is " << numedges << std::endl;
    std::cout << "num verts   is " << numverts << std::endl;
    std::cout << "-----------" << std::endl;
    if(numregs > 0) {
      mmodel = ammodel;
      break;
    }
  }

  /* SET THE ADJACENCIES */
  MG_API_CALL(m, set_adjacency_state(REGION2FACE|
                                     REGION2EDGE|
                                     REGION2VERTEX|
                                     FACE2EDGE|
                                     FACE2VERTEX));
  MG_API_CALL(m, set_reverse_states());
  MG_API_CALL(m, set_adjacency_scope(TOPO_EDGE, SCOPE_FULL));
  MG_API_CALL(m, set_adjacency_scope(TOPO_FACE, SCOPE_FULL));
  MG_API_CALL(m, compute_adjacency());

  gmi_cap_start();
  gmi_register_cap();

  // convert the mesh to apf/mds mesh
  apf::Mesh2* mesh = apf::createMesh(m,g);

  apf::Mesh2* mdsMesh = apf::createMdsMesh(mesh->getModel(), mesh, true);
  apf::disownMdsModel(mdsMesh);

  mdsMesh->writeNative(mdsFileName);

  apf::destroyMesh(mdsMesh);
  apf::destroyMesh(mesh);

  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}
