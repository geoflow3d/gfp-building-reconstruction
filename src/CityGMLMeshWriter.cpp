#include "stepedge_nodes.hpp"

namespace geoflow::nodes::stepedge {

  void CityGMLMeshWriterNode::process() {
    auto& meshes = vector_input("mesh");

    //write citygml header

    for (size_t i=0; i< meshes.size(); ++i) {
      //write citygml piece for one building

      auto& mesh = meshes.get<Mesh>(i);

      // both of these have the same length:
      auto& faces = mesh.get_polygons();
      auto& labels = mesh.get_labels();

      for (size_t j=0 ; j<faces.size(); ++j) {
        //write citygml piece for one semantic surface

        const LinearRing& face = faces[j];
        const int& label = labels[j];
        
        // label
        std::cout << "Label=";
        if (label==FLOOR) std::cout << "FLOOR\n";
        else if (label==ROOF) std::cout << "ROOF\n";
        else if (label==INNERWALL) std::cout << "INNERWALL\n";
        else if (label==OUTERWALL) std::cout << "OUTERWALL\n";

        // exterior coordinates (NB no repetition of first point)
        for (auto& p : face) {
          std::cout << p[0] << " " << p[1] << " " << p[2] << "\n";
        }
        // interior rings
        for (auto& ring : face.interior_rings()) {
          std::cout << "hole coming\n";
          for (auto& p : ring) {
            std::cout << p[0] << " " << p[1] << " " << p[2] << "\n";
          }
        }
      }
    }

    //write citygml footer
  }

}