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

// Let 'vertices' be an array of N pairs (x,y), indexed from 0
// Let 'area' = 0.0
// for i = 0 to N-1, do
//   Let j = (i+1) mod N
//   Let area = area + vertices[i].x * vertices[j].y
//   Let area = area - vertices[i].y * vertices[j].x
// end for
// Return 'area'

#include "stepedge_nodes.hpp"
#include <cmath>

namespace geoflow::nodes::stepedge {

template<typename T> float compute_ring_area(const T& ring) {
  
  size_t n = ring.size();
  float area = 0;
  for(size_t i=0; i<n; ++i) {
    size_t j = (i+1) % n;
    area += ring[i][0] * ring[j][1];
    area -= ring[i][1] * ring[j][0];
  }
  return std::fabs(area) / 2;

}

float compute_polygon_area(const LinearRing& polygon) {
  
  // total exterior ring area
  float area = compute_ring_area(polygon);

  // subtract hole areas
  for (auto& iring : polygon.interior_rings()) {
    area -= compute_ring_area(iring);
  }

  return area;

}

void DataCoverageCalcNode::process() {
  auto& footprint_polygon = input("footprint_polygon").get<LinearRing>();
  auto data_area = input("data_area").get<float>();
  auto& groundparts = vector_input("ground_parts");

  // non_ground_area is sum of the area of the roofparts, thus excluding all groundparts
  auto non_ground_area = compute_polygon_area(footprint_polygon);
  
  for(size_t i=0; i< groundparts.size(); ++i) {
    auto& groundpart = groundparts.get<LinearRing>(i);
    non_ground_area -= compute_polygon_area(groundpart);
  }

  output("data_coverage").set( data_area / non_ground_area );
}

}