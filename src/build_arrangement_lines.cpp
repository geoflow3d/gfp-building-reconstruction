#include "stepedge_nodes.hpp"

namespace geoflow::nodes::stepedge {

void BuildArrFromLinesNode::process() {
  
  auto fp_term = input("footprint");
  linereg::Polygon_with_holes_2 footprint;
  if (fp_term.is_connected_type(typeid(linereg::Polygon_with_holes_2)))
    footprint = fp_term.get<linereg::Polygon_with_holes_2>();
  else {
    auto& lr = fp_term.get<LinearRing&>();
    linereg::Polygon_2 poly2;
    for (auto& p : lr) {
      poly2.push_back(linereg::Point_2(p[0], p[1]));
    }
    footprint = linereg::Polygon_with_holes_2(poly2);
  }

  Arrangement_2 arr_base;
  Face_split_observer obs (arr_base);
  {
    insert(arr_base, footprint.outer_boundary().edges_begin(), footprint.outer_boundary().edges_end());
    // arr_insert_polygon(arr_base, footprint);
    // insert_non_intersecting_curves(arr_base, footprint.edges_begin(), footprint.edges_end());
    if (!footprint.outer_boundary().is_simple()) {
      arr_filter_biggest_face(arr_base, rel_area_thres);
    }
    obs.set_hole_mode(true);
    for(auto hole = footprint.holes_begin(); hole != footprint.holes_end(); ++hole) {
      insert(arr_base, hole->edges_begin(), hole->edges_end());
    }
    obs.set_hole_mode(false);
  }

  // output only empty footprint if there are too many lines/faces (makes the graph-cut optimisation too slow)
  auto& lines_term = vector_input("lines");
  
  typedef std::pair<Point_2, Point_2> PointPair;


  int arr_complexity = lines_term.size();
  output("arr_complexity").set(arr_complexity);

  if (lines_term.is_connected_type(typeid(linereg::Segment_2))) {
    for(size_t i=0; i<lines_term.size(); ++i) {
      auto& s = lines_term.get<linereg::Segment_2>(i);
      insert(arr_base, s);
    }

  } else {

    std::vector<PointPair> segments;
    for(size_t i=0; i<lines_term.size(); ++i) {
      auto& s = lines_term.get<Segment>(i);
      Point_2 a(s[0][0],s[0][1]);
      Point_2 b(s[1][0],s[1][1]);
      segments.push_back(std::make_pair(a,b));
    }

    if (arr_complexity > max_arr_complexity) {
      std::sort(segments.begin(), segments.end(), [](PointPair& a, PointPair& b) {
        return CGAL::squared_distance(a.first,a.second) > CGAL::squared_distance(b.first,b.second);
      });
    }

    {
      std::vector<X_monotone_curve_2> lines;
      size_t i=0;
      for(auto& s : segments) {
        // Check for 0 segment length (quick fix for linux crash)
        if (i++ == max_arr_complexity) break;
        if (CGAL::squared_distance(s.first,s.second) > 0.001) lines.push_back(Segment_2(s.first,s.second));
      }
      insert(arr_base, lines.begin(), lines.end());
    }
  }
  
  //remove dangling edges
  {
    std::vector<Arrangement_2::Halfedge_handle> to_remove;
    for (auto he : arr_base.edge_handles()) {
      if (he->face()==he->twin()->face())
        to_remove.push_back(he);
    }
    for (auto he : to_remove) {
      arr_base.remove_edge(he);
    }
  }

  // if (snap_clean) arr_snapclean(arr_base, snap_dist, snap_detect_only);

  output("arrangement").set(arr_base);
}

} //namespace geoflow::nodes::stepedge