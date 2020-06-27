#include "stepedge_nodes.hpp"

namespace geoflow::nodes::stepedge {


std::tuple<size_t, size_t, size_t> arr_checkzone(Arrangement_2& arr, Segment_2 segment, Arrangement_2::Face_handle& face) {
  std::vector<CGAL::Object> zone_elems;
  
  Arrangement_2::Halfedge_handle edge;
  Arrangement_2::Vertex_handle vertex;
  CGAL::zone(arr, segment, std::back_inserter(zone_elems));
  // std::cout << "Zone has " << zone_elems.size() << " elems\n";
  size_t v_cnt=0, f_cnt=0, e_cnt=0;
  for ( auto& object : zone_elems ) {
    if ( assign(face, object) ) {
      // std::cout << "Zone face, segid= " << face->data().segid << "\n";
      ++f_cnt;
    } else if ( assign(edge, object) ) {
      // std::cout << "Zone edge\n";
      ++e_cnt;
    } else if ( assign(vertex, object) ) {
      // std::cout << "Zone vertex\n";
      ++v_cnt;
    }
  }
  return std::make_tuple(v_cnt, e_cnt, f_cnt);
}

void BuildArrFromRingsExactNode::arr_snapclean_from_fp(Arrangement_2& arr) {
  typedef Arrangement_2::Traits_2 AT;
  typedef std::variant<Arrangement_2::Vertex_handle, Arrangement_2::Halfedge_handle> Candidate;
  typedef std::unordered_map<Arrangement_2::Vertex_handle, std::vector<Candidate>> CandidateMap;
  float snap_dist_sq = snap_dist*snap_dist;

  PointCollection snap_to_v, snap_v;
  SegmentCollection snap_to_e;

  CandidateMap candidates;
  std::vector<Arrangement_2::Vertex_handle> vertices_to_snap;
  // collect non fp vertices
  for (auto& v :  arr.vertex_handles()){
    auto vhe = v->incident_halfedges();
    auto vdone = vhe;
    // check if v is not on the fp boundary
    bool v_is_on_fp = false;
    do {
      v_is_on_fp |= (!vhe->face()->data().in_footprint) || (!vhe->twin()->face()->data().in_footprint);
    } while (++vhe!=vdone);
    if (v_is_on_fp)
      vertices_to_snap.push_back(v);
  }
  // find candidate vertices/edges to snap to for each vertex
  auto obs = Snap_observer(arr);

  for (auto& v : vertices_to_snap) {
    auto vhe = v->incident_halfedges();
    auto vdone = vhe;
    do { //for all incident faces f of v
      auto f = vhe->face();
      if(f->data().in_footprint) {
        auto fhe = f->outer_ccb();
        auto fdone = fhe;
        do { //for all edges in outer_ccb of f
          // only care if one side of the edge is outside the fp, ie the edge is part of the fp
          if (fhe->source() == v || fhe->target() == v) continue;

          bool check_this_edge;
          // both sides of the edge must be in the fp
          check_this_edge = (fhe->face()->data().in_footprint) && (fhe->twin()->face()->data().in_footprint);

          if( check_this_edge ) {
            // compute distance and compare to threshold
            auto s = AT::Segment_2(fhe->source()->point(), fhe->target()->point());
            if (snap_dist_sq > CGAL::squared_distance(v->point(), s)) {
              candidates[v].push_back(fhe);
              snap_to_e.push_back({
                arr3f{float(CGAL::to_double(s.source().x())), float(CGAL::to_double(s.source().y())), 0},
                arr3f{float(CGAL::to_double(s.target().x())), float(CGAL::to_double(s.target().y())), 0}
              });
              snap_v.push_back({float(CGAL::to_double(v->point().x())), float(CGAL::to_double(v->point().y())), 0});
            }
            if (snap_dist_sq > CGAL::squared_distance(v->point(), fhe->source()->point())) {
              candidates[v].push_back(fhe->source());
              snap_to_v.push_back({float(CGAL::to_double(s.source().x())), float(CGAL::to_double(s.source().y())), 0});
              snap_v.push_back({float(CGAL::to_double(v->point().x())), float(CGAL::to_double(v->point().y())), 0});
            }
          }
        } while (++fhe!=fdone);
      }
    } while (++vhe!=vdone);    

    // perform snapping for this vertex
    if (candidates.count(v) && !snap_detect_only) {
      auto& cvec = candidates[v];
      // std::cout << cvec.size() << "\n";
      // merge v with an edge
      if (cvec.size() == 1) {
        if(auto he_ptr = std::get_if<Arrangement_2::Halfedge_handle>(&cvec[0])) {

          auto& source = (*he_ptr)->source()->point();
          auto& target = (*he_ptr)->target()->point();
          auto line = AT::Line_2(source, target);
          auto split_point = line.projection(v->point());
          
          // check zone to target edge
          auto blocking_segment = AT::Segment_2(split_point, v->point());
          Arrangement_2::Face_handle face;
          auto [v_cnt, e_cnt, f_cnt] = arr_checkzone(arr, blocking_segment, face);
          bool empty_zone = (v_cnt==1 && e_cnt==1) && f_cnt==1;

          if (empty_zone && face->data().in_footprint) {
            // split
            if (CGAL::do_intersect(AT::Segment_2(source, target), split_point)) {
              AT::Segment_2 s1(source, split_point);
              AT::Segment_2 s2(split_point, target);
              auto e_split = arr.split_edge((*he_ptr), s1, s2);
              auto v_split = e_split->target();

              // create new edge to split_vertex
              auto he_fix = arr.insert_at_vertices(blocking_segment, v_split, v);
            }
          }
        }
      } else if (cvec.size() > 1) {
        // pick the 1st vertex
        bool found_vertex=false;
        Arrangement_2::Vertex_handle v_target;
        for (auto& obj : cvec) {
          if(auto v_ptr = std::get_if<Arrangement_2::Vertex_handle>(&obj)) {
            v_target = *v_ptr;
            found_vertex = true;
            break;
          }
        }
        if (found_vertex) {
          auto blocking_segment = AT::Segment_2(v_target->point(), v->point());
          Arrangement_2::Face_handle face;
          auto [v_cnt, e_cnt, f_cnt] = arr_checkzone(arr, blocking_segment, face);
          bool empty_zone = (v_cnt==2 && e_cnt==0) && f_cnt==1;
          if (empty_zone) {
            // create new edge to target vertex
            if (face->data().segid==0)
              auto he_fix = arr.insert_at_vertices(blocking_segment, v_target, v);
          }
        }
      }
    }
  }

  output("snap_fp_v").set(snap_v);
  output("snap_fp_to_v").set(snap_to_v);
  output("snap_fp_to_e").set(snap_to_e);
}

void arr2segments(Face_handle& face, LineStringCollection& segments) {
  auto he = face->outer_ccb();
  auto first = he;

  while(true){
    segments.push_back({
      {
        float(CGAL::to_double(he->source()->point().x())),
        float(CGAL::to_double(he->source()->point().y())),
        0
      },{
        float(CGAL::to_double(he->target()->point().x())),
        float(CGAL::to_double(he->target()->point().y())),
        0
      }
    });

    he = he->next();
    if (he==first) break;
  }
  // segments.push_back({
  //   {
  //     float(CGAL::to_double(he->source()->point().x())),
  //     float(CGAL::to_double(he->source()->point().y())),
  //     0
  //   },{
  //     float(CGAL::to_double(he->target()->point().x())),
  //     float(CGAL::to_double(he->target()->point().y())),
  //     0
  //   }
  // });
}

void arr_snapclean(Arrangement_2& arr, float& snap_dist, bool& snap_detect_only) {
  typedef Arrangement_2::Traits_2 AT;
  typedef std::variant<Arrangement_2::Vertex_handle, Arrangement_2::Halfedge_handle> Candidate;
  typedef std::unordered_map<Arrangement_2::Vertex_handle, std::map<double,Candidate>> CandidateMap;
  float snap_dist_sq = snap_dist*snap_dist;

  PointCollection snap_to_v, snap_v;
  SegmentCollection snap_to_e;

  CandidateMap candidates;
  std::vector<Arrangement_2::Vertex_handle> vertices_to_snap;
  // collect non fp vertices
  for (auto& v :  arr.vertex_handles()){
    auto vhe = v->incident_halfedges();
    auto vdone = vhe;
    // check if v is not on the fp boundary
    bool on_fp = false;
    do {
      on_fp |= (!vhe->face()->data().in_footprint) || (!vhe->twin()->face()->data().in_footprint);
    } while (++vhe!=vdone);
    if (!on_fp)
      vertices_to_snap.push_back(v);
  }
  // find candidate vertices/edges to snap to for each vertex
  auto obs = Snap_observer(arr);

  for (auto& v : vertices_to_snap) {
    auto vhe = v->incident_halfedges();
    auto vdone = vhe;
    do { //for all incident faces f of v
      auto f = vhe->face();
      if(f->data().in_footprint) {
        auto fhe = f->outer_ccb();
        auto fdone = fhe;
        do { //for all edges in outer_ccb of f
          // only care if one side of the edge is outside the fp, ie the edge is part of the fp
          bool check_this_edge;
          // one side of the edge must be outside the fp
          check_this_edge = (!fhe->face()->data().in_footprint) || (!fhe->twin()->face()->data().in_footprint);

          if( check_this_edge ) {
            // compute distance and compare to threshold
            auto s = AT::Segment_2(fhe->source()->point(), fhe->target()->point());
            double d = CGAL::to_double(CGAL::squared_distance(v->point(), s));
            if (snap_dist_sq > d) {
              candidates[v][d] = fhe;
              snap_to_e.push_back({
                arr3f{float(CGAL::to_double(s.source().x())), float(CGAL::to_double(s.source().y())), 0},
                arr3f{float(CGAL::to_double(s.target().x())), float(CGAL::to_double(s.target().y())), 0}
              });
              snap_v.push_back({float(CGAL::to_double(v->point().x())), float(CGAL::to_double(v->point().y())), 0});
            }
            d = CGAL::to_double(CGAL::squared_distance(v->point(), fhe->source()->point()));
            if (snap_dist_sq > d) {
              candidates[v][d] = fhe->source();
              snap_to_v.push_back({float(CGAL::to_double(s.source().x())), float(CGAL::to_double(s.source().y())), 0});
              snap_v.push_back({float(CGAL::to_double(v->point().x())), float(CGAL::to_double(v->point().y())), 0});
            }
          }
        } while (++fhe!=fdone);
      }
    } while (++vhe!=vdone);    

    // perform snapping for this vertex
    if (candidates.count(v) && !snap_detect_only) {
      auto& cvec = candidates[v];
      // std::cout << cvec.size() << "\n";
      // merge v with an edge
      if (cvec.size() == 1) {
        if(auto he_ptr = std::get_if<Arrangement_2::Halfedge_handle>(&(cvec.begin()->second))) {

          auto& source = (*he_ptr)->source()->point();
          auto& target = (*he_ptr)->target()->point();
          auto line = AT::Line_2(source, target);
          auto split_point = line.projection(v->point());
          
          // check zone to target edge
          auto blocking_segment = AT::Segment_2(split_point, v->point());
          Arrangement_2::Face_handle face;
          auto [v_cnt, e_cnt, f_cnt] = arr_checkzone(arr, blocking_segment, face);
          bool empty_zone = (v_cnt==1 && e_cnt==1) && f_cnt==1;

          if (empty_zone && face->data().in_footprint && face->data().segid==0) {
            // split
            AT::Segment_2 s1(source, split_point);
            AT::Segment_2 s2(split_point, target);
            if(CGAL::do_intersect(AT::Segment_2(source, target), split_point)) {
              auto e_split = arr.split_edge((*he_ptr), s1, s2);
              auto v_split = e_split->target();

              // create new edge to split_vertex
              auto he_fix = arr.insert_at_vertices(blocking_segment, v_split, v);
            }
          }
        }
      } else if (cvec.size() > 1) {
        // pick the 1st (closest) vertex
        bool found_vertex=false;
        Arrangement_2::Vertex_handle v_target;
        for (auto& [d, obj] : cvec) {
          if(auto v_ptr = std::get_if<Arrangement_2::Vertex_handle>(&obj)) {
            v_target = *v_ptr;
            found_vertex = true;
            break;
          }
        }
        if (found_vertex) {
          auto blocking_segment = AT::Segment_2(v_target->point(), v->point());
          Arrangement_2::Face_handle face;
          auto [v_cnt, e_cnt, f_cnt] = arr_checkzone(arr, blocking_segment, face);
          bool empty_zone = (v_cnt==2 && e_cnt==0) && f_cnt==1;
          if (empty_zone) {
            // create new edge to target vertex
            if (face->data().segid==0)
              auto he_fix = arr.insert_at_vertices(blocking_segment, v_target, v);
          }
        }
      }
    }
  }

  // output("snap_v").set(snap_v);
  // output("snap_to_v").set(snap_to_v);
  // output("snap_to_e").set(snap_to_e);
}

void BuildArrFromRingsExactNode::arr_process(Arrangement_2& arr) {

  if (flood_to_unsegmented) {
    std::map<float, Face_handle> face_map;
    for (auto& face : arr.face_handles()) {
      if (face->data().segid!=0 && face->data().in_footprint)
        face_map[face->data().elevation_avg] = face;
    }
    for (auto& kv : face_map) {
      std::stack<Face_handle> candidate_stack;
      // std::cout << "Growing face with elevation=" << kv.first << "\n";
      auto& cur_data = kv.second->data();
      candidate_stack.push(kv.second);
      while (!candidate_stack.empty()) {
        auto fh = candidate_stack.top(); candidate_stack.pop();
        auto he = fh->outer_ccb();
        auto done = he;
        do {
          // std::cout << &(*curr) << "\n";
          // ignore weird nullptrs (should not be possible...)
          if (he==nullptr) {
            // std::cout << "nullptr detected!\n";
            break;
          }
          // skip blocking edges (see snapping)
          if (!he->data().blocks) {
            auto candidate = he->twin()->face();
            if (candidate->data().segid == 0 && candidate->data().in_footprint) {
              candidate->data() = cur_data;
              candidate_stack.push(candidate);
            }
          }
          he = he->next();
        } while (he != done);
      }
    }
  }
  
  {
    Face_merge_observer obs(arr);
    if (dissolve_stepedges) {
      arr_dissolve_step_edges_naive(arr, step_height_threshold, false);
    }
    //remove edges that have the same segid on both sides
    if (dissolve_edges) {
      arr_dissolve_seg_edges(arr);
    }
  }
  // if (remove_low_holes) {
  //   std::vector<Arrangement_2::Halfedge_handle> to_remove;
  //   for (auto& face : arr.face_handles()) {
  //     if (face->data().segid!=0 && face->data().in_footprint) {
  //       for (auto iccb=face->holes_begin(); iccb != face->holes_end(); ++iccb) {
  //         auto inner_face = (*iccb)->twin()->face();
  //         if (inner_face->data().elevation_avg < face->data().elevation_avg)
  //           to_remove.push_back(inner_face);
  //       }
  //     }
  //   }
  //   for (auto& f : to_remove) {
  //     arr.rem
  //   }
  // }
}


std::pair<double, double> arr_measure_nosegid(Arrangement_2& arr) {
  // check number of faces
  double total_area=0;
  double no_segid_area=0;
  for (auto& fh : arr.face_handles()) {
    if (fh->data().in_footprint) {
      auto poly = arr_cell2polygon(fh);
      double area = CGAL::to_double(CGAL::abs(poly.area()));
      if (fh->data().segid == 0) {
        no_segid_area += area;
      }
      total_area += area;
    }
  }
  return std::make_pair(no_segid_area, no_segid_area/total_area);
}

void BuildArrFromRingsExactNode::arr_assign_pts_to_unsegmented(Arrangement_2& arr, std::vector<Point>& points) {
  typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Point_location;

  std::unordered_map<Face_handle, std::vector<Point>> points_per_face;

  // collect for each face the points it contains
  Point_location pl(arr);
  for (auto& p : points){
    auto obj = pl.locate( Point_2(p.x(), p.y()) );
    if (auto f = boost::get<Face_const_handle>(&obj)) {
      auto fh = arr.non_const_handle(*f);
      if (fh->data().segid==0 && fh->data().in_footprint) {
        points_per_face[fh].push_back(p);
      }
    }
  }
  // find elevation percentile for each face
  for(auto& ppf : points_per_face) {
    double area = CGAL::to_double(CGAL::abs(arr_cell2polygon(ppf.first).area()));
    if (ppf.second.size()/area > extrude_mindensity) {
      std::sort(ppf.second.begin(), ppf.second.end(), [](linedect::Point& p1, linedect::Point& p2) {
        return p1.z() < p2.z();
      });
      auto pid = int(z_percentile*float(ppf.second.size()-1));
      // auto pid = get_percentile(ppf.second, percentile);
      ppf.first->data().segid = -1;
      ppf.first->data().elevation_avg = ppf.second[pid].z();
    }
  }
}

void arr_insert_polygon(Arrangement_2& arr, const linereg::Polygon_2& polygon) {
  auto tr = Traits_2();
  auto compare = tr.compare_xy_2_object();

  auto e = polygon.edges_begin();
  bool source_was_left = (compare(e->source(),e->target()) == CGAL::SMALLER);
  auto he = arr.insert_in_face_interior(*e, arr.unbounded_face());
  auto v_prev = source_was_left ? he->target() : he->source();
  auto first_he = he;
  auto last_e = --polygon.edges_end();

  ++e;
  for (auto eit = e; eit != last_e; ++eit) {
    // Traits_2::Compare_xy_2()(eit->source(),eit->target());
    auto source_is_left = compare(eit->source(),eit->target()) == CGAL::SMALLER;
    if (source_is_left)
      he = arr.insert_from_left_vertex(*eit, v_prev);
    else
      he = arr.insert_from_right_vertex(*eit, v_prev);
    v_prev = he->target();
  }
  auto v1 = first_he->next()==first_he->twin() ? first_he->target() : first_he->source();
  auto v2 = he->next()==he->twin() ? he->target() : he->source();
  arr.insert_at_vertices(*last_e, v1, v2);
}

void BuildArrFromRingsExactNode::process() {
  // Set up vertex data (and buffer(s)) and attribute pointers
  auto rings = input("rings").get<std::unordered_map<size_t, linereg::Polygon_2>>();
  // auto plane_idx = input("plane_idx").get<vec1i>();
  auto points_per_plane = input("pts_per_roofplane").get<IndexedPlanesWithPoints>();

  auto fp_in = input("footprint");
  linereg::Polygon_2 footprint;
  if (fp_in.is_connected_type(typeid(linereg::Polygon_2)))
    footprint = fp_in.get<linereg::Polygon_2>();
  else {
    auto& lr = fp_in.get<LinearRing&>();
    for (auto& p : lr) {
      footprint.push_back(linereg::Point_2(p[0], p[1]));
    }
  }

  Arrangement_2 arr_base;
  {
    Face_index_observer obs (arr_base, true, 0, 0, Plane());
    insert(arr_base, footprint.edges_begin(), footprint.edges_end());
    // arr_insert_polygon(arr_base, footprint);
    // insert_non_intersecting_curves(arr_base, footprint.edges_begin(), footprint.edges_end());
    if (!footprint.is_simple()) {
      arr_filter_biggest_face(arr_base, rel_area_thres);
    }
  }
  // insert step-edge lines
  {
    Arrangement_2 arr_overlay;
    size_t i=0;
    // NOTE: rings and points_per_plane must be aligned!! (matching length and order)
    for (auto& kv : rings) {
      auto plane_id = kv.first;
      if (plane_id<1) continue;
      auto& polygon = kv.second;
      if (polygon.size()>2) {
        
        auto& points = points_per_plane[plane_id].second;
        auto& plane = points_per_plane[plane_id].first;
        if (points.size()==0) continue;
        std::sort(points.begin(), points.end(), [](linedect::Point& p1, linedect::Point& p2) {
          return p1.z() < p2.z();
        });
        int elevation_id = std::floor(z_percentile*float(points.size()-1));

        // wall_planes.push_back(std::make_pair(Plane(s.first, s.second, s.first+Vector(0,0,1)),0));
        Arrangement_2 arr;
        Face_index_observer obs (arr, false, plane_id, points[elevation_id].z(), plane);
        insert(arr, polygon.edges_begin(), polygon.edges_end());
        // arr_insert_polygon(arr, polygon);

        if (!polygon.is_simple()) {
          arr_filter_biggest_face(arr, rel_area_thres);
        }

        Overlay_traits overlay_traits;
        arr_overlay.clear();
        overlay(arr_base, arr, arr_overlay, overlay_traits);
        arr_base = arr_overlay;
      }
    }
  }
  if(snap_clean) arr_snapclean(arr_base, snap_dist, snap_detect_only);
  if(snap_clean_fp) arr_snapclean_from_fp(arr_base);

  if(extrude_unsegmented && points_per_plane.count(-1)) {
    arr_assign_pts_to_unsegmented(arr_base, points_per_plane[-1].second);
  }
  auto nosegid_area = arr_measure_nosegid(arr_base);
  
  arr_process(arr_base);

  // compute vertical errors roof points wrt arrangement
  typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Point_location;
  Point_location pl(arr_base);
  for (auto& [plane_id, plane_pts] : points_per_plane) {
    if (plane_id<1) continue;
    for (auto& p : plane_pts.second) {
      auto obj = pl.locate( Point_2(p.x(), p.y()) );
      if (auto f = boost::get<Face_const_handle>(&obj)) {
        auto fh = arr_base.non_const_handle(*f);
        if (fh->data().segid!=0 && fh->data().in_footprint) {
          double d = CGAL::squared_distance(fh->data().plane, p);
          // std::cerr << d << "\n";
          fh->data().rms_error_to_avg += d;
          ++fh->data().inlier_count;
        }
      }
    }
  }
  for (auto& fh : arr_base.face_handles()) {
    if (fh->data().segid!=0 && fh->data().in_footprint) {
      if (fh->data().inlier_count)
        fh->data().rms_error_to_avg = CGAL::sqrt(fh->data().rms_error_to_avg/fh->data().inlier_count);
    }
  }
  
  arr_is_valid = arr_base.is_valid();
  vcount = arr_base.number_of_vertices();
  ecount = arr_base.number_of_edges();

  LineStringCollection segments;
  for (auto& face: arr_base.face_handles()){
    if (face->data().in_footprint)
      arr2segments(face, segments);
  }
  output("noseg_area_a").set(float(CGAL::sqrt(nosegid_area.first)));
  output("noseg_area_r").set(float(nosegid_area.second));
  output("arr_segments").set(segments);
  output("arrangement").set(arr_base);
}

}