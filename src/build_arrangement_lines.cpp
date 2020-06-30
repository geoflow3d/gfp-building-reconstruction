#include "stepedge_nodes.hpp"

namespace geoflow::nodes::stepedge {

Arrangement_2::Vertex_handle arr_cross_edge(Arrangement_2& arr, Segment_2& segment, Arrangement_2::Halfedge_handle& edge, double& dist_threshold) {
  auto& e_source = edge->source()->point();
  auto& e_target = edge->target()->point();
  std::cout << "arr_cross_edge:" << std::endl;
  std::cout << "\te source " << e_source << std::endl;
  std::cout << "\te target " << e_target << std::endl;
  auto result = CGAL::intersection(Segment_2(e_source, e_target), segment);
  if (result) {
    if (auto p = boost::get<Point_2>(&*result)) {
      std::cout << "\tp " << *p << std::endl;
      if (CGAL::squared_distance(*p, e_source) < dist_threshold) {
        return edge->source();
      } else if (CGAL::squared_distance(*p, e_target) < dist_threshold) {
        return edge->target();
      } else {
        auto e = arr.split_edge(edge, Segment_2(e_source, *p), Segment_2(*p, e_target));
        return e->target(); // e->target is the split vertex
      }
    }
  }
  return Arrangement_2::Vertex_handle();
}

void arr_insert(Arrangement_2& arr, Segment_2 segment, int& dist_threshold_exp) {
  std::list<CGAL::Object> zone_elems;
  
  Arrangement_2::Halfedge_handle edge;
  Arrangement_2::Vertex_handle vertex;
  Arrangement_2::Face_handle face;
  CGAL::zone(arr, segment, std::back_inserter(zone_elems));
  std::cout << "Zone has " << zone_elems.size() << " elems\n";

  if (zone_elems.size()==0) {
    return;
  } else {

    // collect the vertices in the crossing of segment with the arrangement. Some are existing, some are created in this function.
    std::list<Arrangement_2::Vertex_handle> vertices;
    // std::pair<Arrangement_2::Vertex_handle, Arrangement_2::Vertex_handle> vertices;
    auto p_source = segment.min();
    auto p_target = segment.max();
    double dist_threshold = std::pow(10, -dist_threshold_exp);
    dist_threshold *= dist_threshold; // we compare only squared distances later

    // take care of a first face
    if(assign(face, zone_elems.front())) {
      if (zone_elems.size()==1) {
        arr.insert_in_face_interior(segment, face);
        zone_elems.pop_front();
      } else {
        // remove first face and see what comes next
        zone_elems.pop_front();
        if(assign(edge, zone_elems.front())) {
          vertex = arr_cross_edge(arr, segment, edge, dist_threshold);
        } else {
          assign(vertex, zone_elems.front());
        }
        arr.insert_from_right_vertex(Segment_2(segment.min(), vertex->point()), vertex);
        vertices.push_back(vertex);
      }
    }

    // collect and create further vertices
    while (zone_elems.size()) {

      // for a face we only are left with handling the case were face is the last element in the zone, otherwise skip faces
      if ( assign(face, zone_elems.front()) && (zone_elems.size()==1) ) {
        auto& v = vertices.back();
        arr.insert_from_left_vertex(Segment_2(v->point(), segment.max()), v);
        vertices.pop_back(); //vertices should be empty after this

      // compute intersection and split edge
      } else if ( assign(edge, zone_elems.front()) ) {
        vertices.push_back(arr_cross_edge(arr, segment, edge, dist_threshold));

      // vertex is easiest
      } else if ( assign(vertex, zone_elems.front()) ) {
        vertices.push_back(vertex);
      }
      zone_elems.pop_front();
      if (vertices.size()==2) {
        auto& s = vertices.front();
        auto& t = vertices.back();
        std::cout << "insert_at_vertices...\n";
        std::cout << s->point() << std::endl;
        std::cout << t->point() << std::endl;
        if(s==t) {
          std::cout << "same vertex" << std::endl;
        } else {
          arr.insert_at_vertices(Segment_2(s->point(), t->point()), s, t);
        }
        vertices.pop_front();
      }
    }
  }

  // size_t v_cnt=0, f_cnt=0, e_cnt=0;
  // for ( auto& object : zone_elems ) {
  //   if ( assign(face, object) ) {
  //     // std::cout << "Zone face, segid= " << face->data().segid << "\n";
  //     ++f_cnt;
  //   } else if ( assign(edge, object) ) {
  //     // std::cout << "Zone edge\n";
  //     ++e_cnt;
  //   } else if ( assign(vertex, object) ) {
  //     // std::cout << "Zone vertex\n";
  //     ++v_cnt;
  //   }
  // }
}

void BuildArrFromLinesNode::process() {
  
  // prepare footprint segments
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

  // insert footprint segments
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

  auto& lines_term = vector_input("lines");
  
  typedef std::pair<Point_2, Point_2> PointPair;

  int arr_complexity = lines_term.size();
  output("arr_complexity").set(arr_complexity);

  if (lines_term.is_connected_type(typeid(linereg::Segment_2))) {
    for(size_t i=0; i<lines_term.size(); ++i) {
      auto& s = lines_term.get<linereg::Segment_2>(i);
      arr_insert(arr_base, s, dist_threshold_exp);
      // insert(arr_base, s);
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