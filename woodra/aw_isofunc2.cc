#include <fstream>

#include <apf.h>
#include <apfCAP.h>
#include <ma.h>
#include <gmi_mesh.h>
#include <gmi_cap.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>

#include "CapstoneModule.h"
using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;

class AWFunc : public ma::IsotropicFunction {
public:
	AWFunc(ma::Mesh* mesh): m(mesh) {
		rescan();
	}
	virtual double getValue(ma::Entity* v) {
		ma::Vector p = ma::getPosition(m, v);
		return eval(p);
	}

	void rescan() {
		scan_h0();
		scan_xmid();
	}

	void scan_h0() {
		h_0 = ma::getAverageEdgeLength(m);
	}

	void scan_xmid() {
		double x_min = INFINITY, x_max = -INFINITY;
		ma::Iterator* it = m->begin(0);
		ma::Entity* e;
		while (e = m->iterate(it)) {
			ma::Vector v = ma::getPosition(m, e);
			if (v.x() < x_min) x_min = v.x();
			if (v.x() > x_max) x_max = v.x();
		}
		m->end(it);
		x_mid = (x_max + x_min) / 2;
	}

	double eval(const ma::Vector& p) const {
		double x = p.x();
		double a = 0.1542;
		double h_min = h_0/4.0, h_max = h_0 * 2;
		return h_max + (h_min - h_max)*std::exp(-a*std::fabs(x - x_mid)/h_0);
	}

	ma::Mesh* m;
	double h_0;
	double x_mid;
};

int main (int argc, char ** argv) {
	MPI_Init(&argc, &argv);
	PCU_Comm_Init();
	lion_set_verbosity(1);

	ma::Mesh* m;
	CapstoneModule *cs = nullptr;

	if (argc == 2) {
		// CRE
		gmi_cap_start();
		gmi_register_cap();
		cs = new CapstoneModule("test", "Geometry Database : SMLIB", "Mesh Database : Create", "Attribution Database : Create");
		std::vector<std::string> files = {argv[1]};
		M_GModel gmodel = cs->load_files(files);
		MeshDatabaseInterface *mdi = cs->get_mesh();
		std::vector<M_MModel> mmodels;
		MG_API_CALL(mdi, get_associated_mesh_models(gmodel, mmodels));
  PCU_ALWAYS_ASSERT(mmodels.size() == 1);
  MG_API_CALL(mdi, set_current_model(mmodels[0]));

	MG_API_CALL(mdi, set_adjacency_state(REGION2FACE|
				     REGION2EDGE|
				     REGION2VERTEX|
				     FACE2EDGE|
				     FACE2VERTEX));
  MG_API_CALL(mdi, set_reverse_states());
  MG_API_CALL(mdi, set_adjacency_scope(TOPO_EDGE, SCOPE_FULL));
  MG_API_CALL(mdi, set_adjacency_scope(TOPO_FACE, SCOPE_FULL));
  MG_API_CALL(mdi, compute_adjacency());
	m = apf::createMesh(mdi, cs->get_geometry());
	} else if (argc == 3) {
		// dmg mesh
		gmi_register_mesh();
		m = apf::loadMdsMesh(argv[1], argv[2]);

  m->verify();
	} else {
		// Unknown mesh.
		if (PCU_Comm_Self() == 0) {
			fprintf(stderr, "Usage: %s <create_file.cre>\nOR\n", argv[0]);
			fprintf(stderr, "Usage: %s <model.dmg> <mesh.smb>\n", argv[0]);
		}
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
	AWFunc aw(m);
	apf::Field *f = apf::createFieldOn(m, "isofield", apf::SCALAR);
	ma::Iterator *it = m->begin(0);
	for (ma::Entity *e = m->iterate(it); e != 0; e = m->iterate(it)) {
		apf::setScalar(f, e, 0, aw.getValue(e));
	}
	m->end(it);
	apf::writeVtkFiles("before", m);
	apf::destroyField(f);
	const ma::Input* in = ma::configure(m, &aw);
	ma::adapt(in);
	m->verify();
	f = apf::createFieldOn(m, "isofield", apf::SCALAR);
	it = m->begin(0);
	for (ma::Entity *e = m->iterate(it); e != 0; e = m->iterate(it)) {
		apf::setScalar(f, e, 0, aw.getValue(e));
	}
	m->end(it);
	apf::writeVtkFiles("after", m);
	aw.scan_h0();
	apf::destroyField(f);
	m->destroyNative();
	apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
	if (cs != nullptr) {
		delete cs;
		gmi_cap_stop();
	}
	PCU_Comm_Free();
  MPI_Finalize();
}
