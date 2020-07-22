
#include "stepedge_nodes.hpp"
// graph cut
#include "alpha_expansion_graphcut.h"
#include <CGAL/graph_traits_dual_arrangement_2.h>
#include <CGAL/Arr_face_index_map.h>
#include "Alpha_expansion_MaxFlow_tag.h"

#include <algorithm>

class FootprintGraph {
  public:
  typedef typename Arrangement_2::Face_handle vertex_descriptor;
  typedef typename Arrangement_2::Halfedge_handle edge_descriptor;
  typedef std::vector<vertex_descriptor> vertex_container;
  typedef std::vector<edge_descriptor> edge_container;
  typedef vertex_container::const_iterator vertex_iterator;
  typedef edge_container::const_iterator edge_iterator;

  private:
  vertex_container vertices_;
  edge_container edges_;

  public:
  FootprintGraph(vertex_container& vc, edge_container& ec)
    : vertices_(vc), edges_(ec) {}

  friend vertex_descriptor source(edge_descriptor& e, const FootprintGraph& g) {
    return e->face();
  };
  friend vertex_descriptor target(edge_descriptor& e, const FootprintGraph& g) {
    return e->twin()->face();
  };
  friend size_t num_vertices(const FootprintGraph& g) {
    return g.vertices_.size();
  };
  friend std::pair<vertex_iterator, vertex_iterator> vertices(const FootprintGraph& g) {
    return std::make_pair(g.vertices_.begin(), g.vertices_.end());
  };
  friend std::pair<edge_iterator, edge_iterator> edges(const FootprintGraph& g) {
    return std::make_pair(g.edges_.begin(), g.edges_.end());
  };
};

namespace boost {
  template<> struct graph_traits<FootprintGraph> {
    typedef typename Arrangement_2::Face_handle vertex_descriptor;
    typedef typename Arrangement_2::Halfedge_handle edge_descriptor;
    typedef boost::disallow_parallel_edge_tag edge_parallel_category;
    typedef boost::edge_list_graph_tag traversal_category;
    typedef boost::directed_tag directed_category;
  };
}


namespace geoflow::nodes::stepedge {

// A property map that reads/writes the information to/from the extended 
// face.
class Vertex_label_cost_property_map {
public:
  typedef typename Arrangement_2::Face_handle     Face_handle;
  // Boost property type definitions.
  typedef boost::read_write_property_map_tag      category;
  typedef std::vector<double>                     value_type;
  typedef value_type&                             reference;
  typedef Face_handle                             key_type;
  // The get function is required by the property map concept.
  friend reference get(const Vertex_label_cost_property_map&, key_type key)
  { return key->data().vertex_label_cost; }
  // The put function is required by the property map concept.
  friend void put(const Vertex_label_cost_property_map&,
                  key_type key, value_type val)
  { key->data().vertex_label_cost=val; }
};
// A property map that reads/writes the information to/from the extended 
// face.
class Vertex_label_property_map {
public:
  typedef typename Arrangement_2::Face_handle     Face_handle;
  // Boost property type definitions.
  typedef boost::read_write_property_map_tag      category;
  typedef size_t                                  value_type;
  typedef value_type&                             reference;
  typedef Face_handle                             key_type;
  // The get function is required by the property map concept.
  friend reference get(const Vertex_label_property_map&, key_type key)
  { return key->data().label; }
  // The put function is required by the property map concept.
  friend void put(const Vertex_label_property_map&,
                  key_type key, value_type val)
  { key->data().label=val; }
};
// A property map that reads/writes the information to/from the extended 
// face.
class Vertex_index_map {
public:
  typedef typename Arrangement_2::Face_handle     Face_handle;
  // Boost property type definitions.
  typedef boost::readable_property_map_tag        category;
  typedef size_t                                  value_type;
  typedef value_type&                             reference;
  typedef Face_handle                             key_type;
  // The get function is required by the property map concept.
  friend reference get(const Vertex_index_map&, key_type key)
  { return key->data().v_index; }
};
// A property map that reads/writes the information to/from the extended 
// edge.
class Edge_weight_property_map {
public:
  typedef typename Arrangement_2::Halfedge_handle   Halfedge_handle;
  // Boost property type definitions.
  typedef boost::read_write_property_map_tag      category;
  typedef double                                  value_type;
  typedef value_type&                             reference;
  typedef Halfedge_handle                         key_type;
  // The get function is required by the property map concept.
  friend reference get(const Edge_weight_property_map&, key_type key)
  { return key->data().edge_weight; }
  // The put function is required by the property map concept.
  friend void put(const Edge_weight_property_map&,
                  key_type key, value_type val)
  { key->data().edge_weight = val; }
};

inline double edge_length(const Arrangement_2::Halfedge_handle& e)
{
  const double     x1 = CGAL::to_double (e->source()->point().x());
  const double     y1 = CGAL::to_double (e->source()->point().y());
  const double     x2 = CGAL::to_double (e->target()->point().x());
  const double     y2 = CGAL::to_double (e->target()->point().y());
  const double     diff_x = x2 - x1;
  const double     diff_y = y2 - y1;
  return std::sqrt(diff_x*diff_x + diff_y*diff_y);
}

inline double rmse_plane_points(Plane& plane, std::vector<Point>& points) {
  double dist_sum = 0;
  for (auto& p : points) {
    dist_sum += CGAL::squared_distance(plane, p);
  }
  return CGAL::sqrt(dist_sum/points.size());
}
inline double volume_to_plane(Plane& plane, vec3f& points) {
  double volume = 0;
  for (auto& p : points) {
    volume += 
    std::abs( p[2] - (-plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c()) );
  }
  return volume;
}

struct InFootprintFaceFilter {
  InFootprintFaceFilter() {}

  template <typename Face_handle>
  bool operator()(const Face_handle& f) const {
    return f->data().in_footprint;
  }
};
struct InFootprintEdgeFilter {
  InFootprintEdgeFilter() {}

  template <typename Edge_handle>
  bool operator()(const Edge_handle& e) const {
    bool infp_l = e->twin()->face()->data().in_footprint;
    bool infp_r = e->face()->data().in_footprint;
    return infp_l && infp_r;
  }
};

void OptimiseArrangmentGridNode::process() {

  auto arr = input("arrangement").get<Arrangement_2>();
  auto& planes = input("pts_per_roofplane").get<IndexedPlanesWithPoints>();
  auto& heightfield = input("heightfield").get<RasterTools::Raster>();

  std::vector<std::tuple<Plane, size_t>> points_per_plane; // plane, points, seg_id
  for (auto& [plane_id, plane_pts] : planes) {
    if (plane_id<1) continue; // ignore unclassified points
    points_per_plane.push_back(std::make_tuple(plane_pts.first, plane_id));
  }
  size_t roofplane_cnt = points_per_plane.size();
  if (use_ground) {
    auto& ground_planes = input("ground_pts_per_roofplane").get<IndexedPlanesWithPoints>();
    for (auto& [plane_id, plane_pts] : ground_planes) {
      if (plane_id<1) continue; // ignore unclassified points
      size_t plane_id_ = (roofplane_cnt-1) + plane_id;
      points_per_plane.push_back(std::make_tuple(plane_pts.first, plane_id_));
    }
  }

  // compute vertex_label_cost (data term)
  // 1 compute for each face the error to each plane
  double max_cost = 0;
  size_t face_i=0;
  size_t label = 0;
  std::vector<Face_handle> faces;
  for (auto face: arr.face_handles()) {
    if(face->data().in_footprint) {
      vec2f polygon;
      arrangementface_to_polygon(face, polygon);
      auto height_points = heightfield.rasterise_polygon(polygon, false);
      
      // calculate percentile elevation
      if (height_points.size()==0) {
        face->data().elevation_avg = 0;
      } else {
        std::sort(height_points.begin(), height_points.end(), [](auto& p1, auto& p2) {
          return p1[2] < p2[2];
        });
        int elevation_id = std::floor(z_percentile*float(height_points.size()-1));
        face->data().elevation_avg = height_points[elevation_id][2];
      }
      face->data().inlier_count = height_points.size();

      for (auto& [plane, plane_id] : points_per_plane) {
        double volume = data_multiplier * volume_to_plane(plane, height_points);
        face->data().vertex_label_cost.push_back(volume);
        max_cost = std::max(max_cost, volume);
      }
      face->data().v_index = face_i++;
      faces.push_back(face);
      if(preset_labels)
        face->data().label = label; //assign an initial label
      // also compute average elevation for all inliers (needed for LoD1.3 later)
    }
    ++label;
  }
  // normalise
  if(do_normalise) {
    for (auto face : faces) {
      for (auto& c : face->data().vertex_label_cost) {
        c = (c/max_cost);
      }
    }
  }

  // std::cerr << "  Graph-cut with " << faces.size() << " faces and " << points_per_plane.size() << " plane labels\n";
  // don't bother if there are no planes or faces
  if (roofplane_cnt==0 || faces.size()==0) {
    arr_dissolve_fp(arr, true, true);
    output("arrangement").set(arr);
    return;
  }

  // compute edge_weights (smoothness term)
  double max_weight = 0;
  std::vector<Halfedge_handle> edges;
  for (auto edge : arr.edge_handles()) {
    bool fp_u = edge->twin()->face()->data().in_footprint;
    bool fp_l = edge->face()->data().in_footprint;
    if (fp_u && fp_l) { // only edges with both neighbour faces inside the footprint
      double l = smoothness_multiplier * edge_length(edge);
      edge->data().edge_weight = l;
      edge->twin()->data().edge_weight = l;
      max_weight = std::max(max_weight, l);
      edges.push_back(edge);
      edges.push_back(edge->twin());
    }
  }
  // normalise
  if(do_normalise) {
    for (auto edge : edges) {
      bool fp_l = edge->face()->data().in_footprint;
      double n_w =  edge->data().edge_weight/max_weight;
      edge->data().edge_weight = n_w;
    }
  }
  
  FootprintGraph graph(faces, edges);
  
  // assign initial labels?
  // auto graph = Dual_arrangement(arr);
  // InFootprintFaceFilter face_filter;
  // InFootprintEdgeFilter edge_filter;
  // boost::filtered_graph<Dual_arrangement, InFootprintEdgeFilter, InFootprintFaceFilter >
  //   filtered_graph(graph, edge_filter, face_filter);

  // boost::adjacency_list<> g(N);

  double result;

  if (graph_cut_impl==0) {
    result = CGAL::alpha_expansion_graphcut(
      graph, 
      Edge_weight_property_map(),
      Vertex_index_map(),
      Vertex_label_cost_property_map(),
      Vertex_label_property_map(),
      CGAL::Alpha_expansion_boost_adjacency_list_tag(),
      n_iterations
    );
  } else if (graph_cut_impl==1) {
    result = CGAL::alpha_expansion_graphcut(
      graph, 
      Edge_weight_property_map(),
      Vertex_index_map(),
      Vertex_label_cost_property_map(),
      Vertex_label_property_map(),
      CGAL::Alpha_expansion_boost_compressed_sparse_row_tag(),
      n_iterations
    );
  } else if (graph_cut_impl==2) {
    result = CGAL::alpha_expansion_graphcut(
      graph, 
      Edge_weight_property_map(),
      Vertex_index_map(),
      Vertex_label_cost_property_map(),
      Vertex_label_property_map(),
      CGAL::Alpha_expansion_MaxFlow_tag(),
      n_iterations
    );
  }

  //  assign planes from label map to arrangement faces
  for (auto& face : faces) {
    size_t i = face->data().label;
    face->data().plane = std::get<0>(points_per_plane[i]);
    face->data().segid = std::get<1>(points_per_plane[i]);
    face->data().rms_error_to_avg = face->data().vertex_label_cost[i];
    if(i>=roofplane_cnt && label_ground_outside_fp) {
      face->data().in_footprint = false;
    }
  }

  output("arrangement").set(arr);
}
void OptimiseArrangmentNode::process() {

  auto arr = input("arrangement").get<Arrangement_2>();
  auto& planes = input("pts_per_roofplane").get<IndexedPlanesWithPoints>();

  std::vector<std::tuple<Plane, std::vector<Point>, size_t, float>> points_per_plane; // plane, points, seg_id, average elevation
  for (auto& [plane_id, plane_pts] : planes) {
    if (plane_id<1) continue; // ignore unclassified points
    // also calculate percentile elevation for all inliers of this plane
    auto points = plane_pts.second;
    std::sort(points.begin(), points.end(), [](auto& p1, auto& p2) {
      return p1.z() < p2.z();
    });
    int elevation_id = std::floor(z_percentile*float(points.size()-1));

    points_per_plane.push_back(std::make_tuple(plane_pts.first, points, plane_id, points[elevation_id].z()));
  }

  // note we can try to use Face_filtered_graph to exclude the inf face

  // compute vertex_label_cost (data term)
  // 1 assign points to faces
  typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Point_location;
  Point_location pl(arr);
  size_t label = 0;
  for (auto& [plane, pts, plane_id, elevation_avg] : points_per_plane) {
    for (auto& p : pts) {
      auto obj = pl.locate( Point_2(p.x(), p.y()) );
      if (auto f = boost::get<Face_const_handle>(&obj)) {
        auto fh = arr.non_const_handle(*f);
        if (fh->data().in_footprint) {
          fh->data().points.push_back(p);
          if(preset_labels)
            fh->data().label = label; //assign an initial label
        }
      }
    }
    ++label;
  }

  // 2 compute for each face the error to each plane
  double max_cost = 0;
  size_t face_i=0;
  std::vector<Face_handle> faces;
  for (auto face: arr.face_handles()) {
    if(face->data().in_footprint) {
      for (auto& [plane, pts, plane_id, elevation_avg] : points_per_plane) {
        double d = rmse_plane_points(plane, face->data().points);
        face->data().vertex_label_cost.push_back(d);
        max_cost = std::max(max_cost, d);
      }
      face->data().v_index = face_i++;
      faces.push_back(face);
      // also compute average elevation for all inliers (needed for LoD1.3 later)
    }
  }
  // std::cerr << "  Graph-cut with " << faces.size() << " faces and " << points_per_plane.size() << " plane labels\n";
  // if (faces.size() > max_face_complexity) {

  //   return
  // }
  // normalise
  for (auto face : faces) {
    for (auto& c : face->data().vertex_label_cost) {
      c = data_multiplier * (c/max_cost);
    }
  }

  // compute edge_weights (smoothness term)
  double max_weight = 0;
  std::vector<Halfedge_handle> edges;
  for (auto edge : arr.edge_handles()) {
    bool fp_u = edge->twin()->face()->data().in_footprint;
    bool fp_l = edge->face()->data().in_footprint;
    if (fp_u && fp_l) { // only edges with both neighbour faces inside the footprint
      double l = smoothness_multiplier * edge_length(edge);
      edge->data().edge_weight = l;
      edge->twin()->data().edge_weight = l;
      max_weight = std::max(max_weight, l);
      edges.push_back(edge);
      edges.push_back(edge->twin());
    }
  }
  // normalise
  for (auto edge : edges) {
    bool fp_l = edge->face()->data().in_footprint;
    double n_w =  edge->data().edge_weight/max_weight;
    edge->data().edge_weight = n_w;
  }

  FootprintGraph graph(faces, edges);
  
  // assign initial labels?
  // auto graph = Dual_arrangement(arr);
  // InFootprintFaceFilter face_filter;
  // InFootprintEdgeFilter edge_filter;
  // boost::filtered_graph<Dual_arrangement, InFootprintEdgeFilter, InFootprintFaceFilter >
  //   filtered_graph(graph, edge_filter, face_filter);

  // boost::adjacency_list<> g(N);

  double result;

  if (graph_cut_impl==0) {
    result = CGAL::alpha_expansion_graphcut(
      graph, 
      Edge_weight_property_map(),
      Vertex_index_map(),
      Vertex_label_cost_property_map(),
      Vertex_label_property_map(),
      CGAL::Alpha_expansion_boost_adjacency_list_tag(),
      n_iterations
    );
  } else if (graph_cut_impl==1) {
    result = CGAL::alpha_expansion_graphcut(
      graph, 
      Edge_weight_property_map(),
      Vertex_index_map(),
      Vertex_label_cost_property_map(),
      Vertex_label_property_map(),
      CGAL::Alpha_expansion_boost_compressed_sparse_row_tag(),
      n_iterations
    );
  } else if (graph_cut_impl==2) {
    result = CGAL::alpha_expansion_graphcut(
      graph, 
      Edge_weight_property_map(),
      Vertex_index_map(),
      Vertex_label_cost_property_map(),
      Vertex_label_property_map(),
      CGAL::Alpha_expansion_MaxFlow_tag(),
      n_iterations
    );
  }

  //  assign planes from label map to arrangement faces
  for (auto& face : faces) {
    size_t i = face->data().label;
    face->data().plane = std::get<0>(points_per_plane[i]);
    face->data().segid = std::get<2>(points_per_plane[i]);
    face->data().elevation_avg = std::get<3>(points_per_plane[i]);
    face->data().inlier_count = std::get<1>(points_per_plane[i]).size();
    face->data().rms_error_to_avg = face->data().vertex_label_cost[i];
  }

  output("arrangement").set(arr);
}

}
