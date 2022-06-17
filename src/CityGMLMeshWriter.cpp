// This file is part of gfp-building-reconstruction
// Copyright (C) 2018-2022 Ravi Peters

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include "stepedge_nodes.hpp"

namespace geoflow::nodes::stepedge {

  void CityGMLMeshWriterNode::process() {
    auto& meshes = vector_input("mesh");

    // get output file path
    auto filepath = manager.substitute_globals(filepath_);

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