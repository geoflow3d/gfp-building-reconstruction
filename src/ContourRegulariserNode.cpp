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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_contours.h>

#include <fstream>
#include <iostream>

namespace geoflow::nodes::stepedge {

  void ContourRegulariserNode::process() {
    // Typedefs.
    using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
    using FT      = typename Kernel::FT;
    using Point_2 = typename Kernel::Point_2;
    using Contour = std::vector<Point_2>;
    using Contour_directions =
      CGAL::Shape_regularization::Contours::Multiple_directions_2<Kernel, Contour>;


    auto& polygons = vector_input("polygons");
    auto& opolygons = vector_output("regularised_polygons");

    // Set parameters.
    const FT min_length_2 = FT(min_length);
    const FT max_angle_2 = FT(max_angle);
    const FT max_offset_2 = FT(max_offset);
    
    for (size_t i = 0; i < polygons.size(); ++i) {
      std::vector<Point_2> contour;
      auto& polygon = polygons.get<LinearRing&>(i);
      for (auto& p : polygon) {
        contour.push_back(Point_2(p[0], p[1]));
      }

      // Regularize.
      const bool is_closed = true;
      Contour_directions directions(
        contour, is_closed, CGAL::parameters::
        minimum_length(min_length_2).maximum_angle(max_angle_2));
      std::vector<Point_2> regularized;
      CGAL::Shape_regularization::Contours::regularize_closed_contour(
        contour, directions, std::back_inserter(regularized),
        CGAL::parameters::maximum_offset(max_offset_2));
      std::cout << "* number of directions = " <<
        directions.number_of_directions() << std::endl;
      
      LinearRing opoly;
      for (auto& point : regularized) {
        opoly.push_back(
          {
            float(point.x()),
            float(point.y()),
            0
          }
        );
      }
      opolygons.push_back(opoly);
    }

  }
}