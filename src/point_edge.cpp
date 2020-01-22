#include "point_edge.h"
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <unordered_map>
#include <boost/tuple/tuple.hpp>

#include <CGAL/Handle_hash_function.h>

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

void arr_dissolve_step_edges_naive(Arrangement_2& arr, float step_height_threshold, bool compute_on_edge)
{
  std::vector<Arrangement_2::Halfedge_handle> to_remove;
  for (auto& edge : arr.edge_handles()) {
    auto f1 = edge->face();
    auto f2 = edge->twin()->face();

    if((f1->data().in_footprint && f2->data().in_footprint) && (f1->data().segid!=0 && f2->data().segid!=0)) {
      double d;
      if (compute_on_edge) {
        auto& s = edge->source()->point();
        auto& t = edge->target()->point();
        auto& pl1 = f1->data().plane;
        auto& pl2 = f2->data().plane;
        double h1_pl1 = CGAL::to_double((pl1.a()*s.x() + pl1.b()*s.y() + pl1.d()) / (-pl1.c()));
        double h2_pl1 = CGAL::to_double((pl1.a()*t.x() + pl1.b()*t.y() + pl1.d()) / (-pl1.c()));
        double h1_pl2 = CGAL::to_double((pl2.a()*s.x() + pl2.b()*s.y() + pl2.d()) / (-pl2.c()));
        double h2_pl2 = CGAL::to_double((pl2.a()*t.x() + pl2.b()*t.y() + pl2.d()) / (-pl2.c()));
        d = std::max(std::abs(h1_pl1-h1_pl2), std::abs(h2_pl1-h2_pl2));
      } else {
        d = std::abs(f1->data().elevation_avg - f2->data().elevation_avg);
      }
      if(d < step_height_threshold){
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

auto HandleHash = CGAL::Handle_hash_function{};
void arr_dissolve_step_edges(Arrangement_2& arr, float step_height_threshold)
{
  struct FacePair {
      Arrangement_2::Face_handle f_lo;
      Arrangement_2::Face_handle f_hi;

      FacePair(){};
      FacePair(Arrangement_2::Face_handle f1, Arrangement_2::Face_handle f2) {
        if (HandleHash(f1) < HandleHash(f1)) {
          f_lo = f1;
          f_hi = f2;
        } else {
          f_lo = f2;
          f_hi = f1;
        }
      };
  };
    
  struct KeyEqual {
    bool operator()(const FacePair& lhs, const FacePair& rhs) const
    {
      return lhs.f_hi == rhs.f_hi && lhs.f_lo == rhs.f_lo;
    }
  };
  struct KeyHash
  {
    std::size_t operator()(FacePair const& p) const
    {
      std::size_t h1 = HandleHash(p.f_lo);
      std::size_t h2 = HandleHash(p.f_hi);
      return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion)
    }
  };
  
  std::unordered_map<
    FacePair, 
    std::vector<Arrangement_2::Halfedge_handle>,
    KeyHash, 
    KeyEqual
  > step_boundaries;

  while (true) {
    double d_min = step_height_threshold;
    step_boundaries.clear();
    for (auto& edge : arr.edge_handles()) {
      auto f1 = edge->face();
      auto f2 = edge->twin()->face();
      if((f1->data().in_footprint && f2->data().in_footprint) && (f1->data().segid!=0 && f2->data().segid!=0)) {
        step_boundaries[FacePair(f1,f2)].push_back(edge);
      }
    }
    FacePair facepair_min;
    for (auto& [faces, edges] : step_boundaries) {
      double d = std::abs(faces.f_hi->data().elevation_avg - faces.f_lo->data().elevation_avg);
      if (d < d_min) {
        d_min = d;
        facepair_min = faces;
      }
    }
    if (d_min == step_height_threshold) break;
    std::vector<Halfedge_handle> to_remove;
    for (auto& edge : step_boundaries[facepair_min]) {
      arr.remove_edge(edge);
    }
  }


}