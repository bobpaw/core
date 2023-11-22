#include <CapstoneModule.h>

#include <CreateMG_Variable.h>
#include <CreateMG_Framework_Geometry.h>

int main () {
  CapstoneModule cs("list_readers", "Geometry Database : SMLIB",
                    "Mesh Database : Create", "Attribution Database : Create");
  CreateMG::v_string readers;
  cs.list_readers(readers);
  for (auto r : readers) {
    std::cout << r << std::endl;
  }
  CreateMG::v_string functions;
  cs.list_functions(functions);
  for (auto f : functions) {
    std::string info;
    std::cout << f << std::endl;
    cs.get_function_info(f, info);
    std::cout << "- " << info << std::endl;
  }

  CreateMG::v_string meshers;
  
  CreateMG::Geometry::GeometryDatabaseInterface *gdi = cs.get_geometry();
  CreateMG::M_GModel gmodel;
  MG_API_CALL(gdi, create_model(gmodel));
  CreateMG::Variables vars, out;
  vars.add("Model", gmodel);
  cs.execute(std::string("GetValidMeshers"), vars, &out);
  out.get("Meshers")->get_value(meshers);
  for (auto m : meshers) {
    std::cout << m << std::endl;
  }
}
 
