#include <apf.h>
#include <ma.h>
#include <gmi_mesh.h>
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
	PCU_ALWAYS_ASSERT(argc == 3);
	const char* modelFile = argv[1];
	const char* meshFile = argv[2];
	MPI_Init(&argc, &argv);
	PCU_Comm_Init();
	lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
	MS_init();
	SimModel_start();
	Sim_readLicenseFile(0);
	gmi_sim_start();
	gmi_register_sim();
#endif
	gmi_register_mesh();
	ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
	apf::writeVtkFiles("before", m);
	const ma::Input* in = ma::configureUniformRefine(m, 1);
	ma::adapt(in);
	AWFunc aw(m);
	apf::Field *f = apf::createFieldOn(m, "isofield", apf::SCALAR);
	ma::Iterator *it = m->begin(0);
	for (ma::Entity *e = m->iterate(it); e != 0; e = m->iterate(it)) {
		apf::setScalar(f, e, 0, aw.getValue(e));
	}
	m->end(it);
	apf::writeVtkFiles("uniform", m);
	apf::destroyField(f);
	in = ma::configure(m, &aw);
	ma::adapt(in);
	m->verify();
	f = apf::createFieldOn(m, "isofield", apf::SCALAR);
	it = m->begin(0);
	for (ma::Entity *e = m->iterate(it); e != 0; e = m->iterate(it)) {
		apf::setScalar(f, e, 0, aw.getValue(e));
	}
	m->end(it);
	apf::writeVtkFiles("after", m);
	apf::destroyField(f);
	m->destroyNative();
	apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
	PCU_Comm_Free();
  MPI_Finalize();
}
