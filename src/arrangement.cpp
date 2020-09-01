#include "arrangement.hpp"

template<typename E, typename P> void ccb_to_polygon_3(E he, P& polygon, double h=0) {
  auto first = he;

  while(true){
    // if (!he->source()- at_infinity())
      polygon.push_back({
        float(CGAL::to_double(he->source()->point().x())),
        float(CGAL::to_double(he->source()->point().y())),
        float(h)
      });

    he = he->next();
    if (he==first) break;
  // }
  }
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
void arrangementface_to_polygon(Face_handle face, geoflow::LinearRing& polygon, double h){
  // if(extract_face){ // ie it is a face on the interior of the footprint
  auto he = face->outer_ccb();
  ccb_to_polygon_3(he, polygon);

  for (auto ccb = face->inner_ccbs_begin(); ccb != face->inner_ccbs_end(); ++ccb) {
    geoflow::vec3f ring;
    ccb_to_polygon_3(*ccb, ring, h);
    polygon.interior_rings().push_back(ring);
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
        d = std::abs(f1->data().elevation_70p - f2->data().elevation_70p);
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
      double d = std::abs(faces.f_hi->data().elevation_70p - faces.f_lo->data().elevation_70p);
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

void arr_snap_duplicates(Arrangement_2& arr, double dupe_threshold) {
  std::vector<Arrangement_2::Halfedge_handle> to_remove;
  double dupe_threshold_sq = dupe_threshold*dupe_threshold;
  for (auto he : arr.edge_handles()) {
    auto source = he->source();
    auto target = he->target();
    if (CGAL::squared_distance(source->point(), target->point()) < dupe_threshold_sq) {
      if ((source->degree()==2 && target->degree()>2) || (target->degree()==2 && source->degree()>2))
        to_remove.push_back(he);
      else
        std::cout << "skipping and edge in duplicate snapping. Degrees are " << target->degree() << " and " << source->degree() << "\n";
    }
  }
  for (auto he : to_remove) {
    Vertex_handle vy, v_other;
    Halfedge_handle he_other;
    auto source = he->source();
    auto target = he->target();
    if (source->degree()==2 && target->degree()>2) {
      vy = target;
      he_other = he->prev();
      v_other = he_other->source();
    } else if (target->degree()==2 && source->degree()>2) {
      vy = source;
      he_other = he->next();
      v_other = he_other->target();
    }
    arr.merge_edge(he, he_other, Segment_2(vy->point(), v_other->point()));
  }
}


// {
//   for (auto& v :  arr.vertex_handles()){
//     auto vhe = v->incident_halfedges();
//     auto vdone = vhe;
//     // check if v is not on the fp boundary
//     bool on_fp = false;
//     do {
//       on_fp |= (!vhe->face()->data().in_footprint) || (!vhe->twin()->face()->data().in_footprint);
//     } while (++vhe!=vdone);
//     if (!on_fp)
//       vertices_to_snap.push_back(v);
//   }
// }

void arr_dissolve_fp(Arrangement_2& arr, bool inside, bool outside) {
  {
    std::vector<Arrangement_2::Halfedge_handle> to_remove;
    for (auto he : arr.edge_handles()) {
      auto d1 = he->face()->data();
      auto d2 = he->twin()->face()->data();
      if(outside)
        if (!d1.in_footprint && !d2.in_footprint)
          to_remove.push_back(he);
      if(inside)
        if (d1.in_footprint && d2.in_footprint)
          to_remove.push_back(he);
    }
    for (auto he : to_remove) {
      arr.remove_edge(he);
    }
  }
  // cleanup vertices with degree==2
  {
    std::vector<Arrangement_2::Vertex_handle> to_remove;
    for (auto v : arr.vertex_handles()) {
      if(v->degree()==2)
        to_remove.push_back(v);
    }
    for (auto v : to_remove)
      CGAL::remove_vertex(arr, v);
  }
}

void arr_filter_biggest_face(Arrangement_2& arr, const float& rel_area_thres) {
  // check number of faces
  typedef std::pair<Polygon_2, double> polyar;
  std::vector<polyar> polygons;
  double total_area=0;
  for (auto& fh : arr.face_handles()) {
    if (fh->data().segid != 0 || fh->data().in_footprint == true) {
      auto poly = arr_cell2polygon(fh);
      double area = CGAL::to_double(CGAL::abs(poly.area()));
      total_area += area;
      polygons.push_back(std::make_pair(poly, area));
    }
  }
  std::sort(polygons.begin(), polygons.end(), [](const polyar& a, const polyar& b) {
    return a.second < b.second;   
  });
  arr.clear();
  for (auto& poly_a : polygons) {
    if (poly_a.second > rel_area_thres * total_area)
      insert_non_intersecting_curves(arr, poly_a.first.edges_begin(), poly_a.first.edges_end());
  }
}

Polygon_2 arr_cell2polygon(const Face_handle& fh) {
  Polygon_2 poly;
  auto he = fh->outer_ccb();
  auto first = he;
  do {
    poly.push_back(he->target()->point());
    he = he->next();
  } while (he!=first);
  return poly;
}