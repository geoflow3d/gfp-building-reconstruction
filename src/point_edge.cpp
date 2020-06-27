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