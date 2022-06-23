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
#include "point_edge.h"
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <unordered_map>
#include <boost/tuple/tuple.hpp>

#include <CGAL/Handle_hash_function.h>

Polygon_2 ring_to_cgal_polygon(geoflow::LinearRing& ring) {
  typedef AK::Point_2 Point_2;
  std::vector<Point_2> footprint_pts;
  for (auto p : ring) {
    footprint_pts.push_back(Point_2(p[0], p[1]));
  }
  return Polygon_2(footprint_pts.begin(), footprint_pts.end());
}