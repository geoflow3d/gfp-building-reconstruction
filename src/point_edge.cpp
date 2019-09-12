#include "point_edge.h"
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <unordered_map>
#include <boost/tuple/tuple.hpp>

Polygon_2 ring_to_cgal_polygon(geoflow::LinearRing& ring) {
  std::vector<Point_2> footprint_pts;
  for (auto p : ring) {
    footprint_pts.push_back(Point_2(p[0], p[1]));
  }
  return Polygon_2(footprint_pts.begin(), footprint_pts.end());
}


void arrangementface_to_polygon(Face_handle face, vec2f& polygons){
  // if(extract_face){ // ie it is a face on the interior of the footprint
  auto he = face->outer_ccb();
  auto first = he;

  while(true){
    // if (!he->source()- at_infinity())
      polygons.push_back({
        float(CGAL::to_double(he->source()->point().x())),
        float(CGAL::to_double(he->source()->point().y()))
      });

    he = he->next();
    if (he==first) break;
  // }
  }
}

// helper functions
void arr_dissolve_seg_edges(Arrangement_2& arr)
{
  std::vector<Halfedge_handle> to_remove;
  for (auto he : arr.edge_handles()) {
    auto d1 = he->face()->data();
    auto d2 = he->twin()->face()->data();
    if ((d1.segid == d2.segid ) && (d1.in_footprint && d2.in_footprint) && d1.segid != 0)
      to_remove.push_back(he);
  }
  for (auto he : to_remove) {
    arr.remove_edge(he);
  }
}

void arr_dissolve_step_edges(Arrangement_2& arr, float step_height_threshold)
{
  std::vector<Arrangement_2::Halfedge_handle> to_remove;
  for (auto& edge : arr.edge_handles()) {
    auto f1 = edge->face();
    auto f2 = edge->twin()->face();

    if((f1->data().in_footprint && f2->data().in_footprint) && (f1->data().segid!=0 && f2->data().segid!=0)) {
      if(std::abs(f1->data().elevation_avg - f2->data().elevation_avg) < step_height_threshold){
        // Face_merge_observer takes care of data merge
        // if (f2->data().elevation_avg < f1->data().elevation_avg) {
        //   f2->data()= f1->data();
        // } else {
        //   f1->data() = f2->data();
        // }
        to_remove.push_back(edge);
      }
    }
  }
  for (auto edge : to_remove) {
    arr.remove_edge(edge);
  }
}