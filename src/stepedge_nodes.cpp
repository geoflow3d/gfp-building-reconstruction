// #include "nlohmann/json.hpp"
// using json = nlohmann::json;

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <cmath>
#include <boost/tuple/tuple.hpp>

#include <numeric>

// line simplification
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

// 2d union
#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/Direction_3.h>

#include "stepedge_nodes.hpp"

typedef std::array<float,3> vertex;
vertex get_normal(vertex v0, vertex v1, vertex v2) {
    // assuming ccw winding order
    auto a = glm::make_vec3(v0.data());
    auto b = glm::make_vec3(v1.data());
    auto c = glm::make_vec3(v2.data());
    auto n = glm::cross(b-a, c-b);
    return {n.x,n.y,n.z};
}

// interval list
#include "interval.hpp"

#include <unordered_set>
#include <stack>
#include <utility>

namespace geoflow::nodes::stepedge {

void Ring2SegmentsNode::process() {
  auto rings = input("rings").get<LinearRingCollection>();
  SegmentCollection segments;
  std::vector<std::vector<size_t>> ring_idx(rings.size());
  size_t ring_i=0;
  size_t seg_i=0;
  for (auto& ring : rings) {
    for (size_t i=0; i< ring.size()-1; ++i) {
      segments.push_back({ring[i], ring[i+1]});
      ring_idx[ring_i].push_back(seg_i++);
    }
    segments.push_back({ring[ring.size()-1], ring[0]});
    ring_idx[ring_i++].push_back(seg_i++);
  }
  output("edge_segments").set(segments);
  output("ring_idx").set(ring_idx);
}

void PointCloudMeanZNode::process(){
  auto point_clouds = input("point_clouds").get<std::vector<PointCollection>>();

  vec1f heights;
  for (auto& point_cloud : point_clouds) {
    double sum_elevation = 0;
    for (auto& p : point_cloud) {
      sum_elevation += p[2];
    }
    heights.push_back(sum_elevation/point_cloud.size());
  }
  output("heights").set(heights);
}

void PolygonExtruderNode::process() {
  auto rings = input("polygons").get<LinearRingCollection>();
  auto heights = input("heights").get<vec1f>();

  LinearRingCollection rings_3d;
  vec1i surf_type;
  for (size_t i=0; i<rings.size(); ++i) {
    auto h = heights[i];
    //floor
    rings_3d.push_back(rings[i]);
    surf_type.push_back(0);
    //roof
    vec3f r = rings[i];
    for (auto& p : r) p[2] = h;
    rings_3d.push_back(r);
    surf_type.push_back(2);
    //walls
    size_t j_prev = r.size()-1;
    for (size_t j=0; j<r.size(); ++j) {
      LinearRing wall;
      wall.push_back(rings[i][j]);
      wall.push_back(rings[i][j_prev]);
      auto ha = r[j_prev];
      auto hb = r[j];
      wall.push_back(ha);
      wall.push_back(hb);
      surf_type.push_back(1);
      rings_3d.push_back(wall);
      j_prev=j;
    }
  }

  output("rings_3d").set(rings_3d);
  output("ring_types").set(surf_type);
}

void LOD1ExtruderNode::process() {
  // assume ring is CCW oriented
  auto ring = input("polygon").get<LinearRing>();
  auto h_floor = input("floor_elevation").get<float>();

  auto& rings_3d = vector_output("3d_polygons");
  vec1i surf_type;

  Mesh mesh;
  // mesh.create_attribute_field<int>("surface_type");

  //floor
  LinearRing r_floor = ring;
  for (auto& p : r_floor) p[2] = h_floor;
  for (auto& hole_floor : r_floor.interior_rings()) {
    for (auto& p : hole_floor) p[2] = h_floor;
  }

  //roof
  if(input("roof_elevation").get_data_vec()[0].has_value()) {
    float h_roof = input("roof_elevation").get<float>();
    LinearRing r_roof = ring;
    for (auto& p : r_roof) p[2] = h_roof;
    for (auto& hole_floor : r_roof.interior_rings()) {
      for (auto& p : hole_floor) p[2] = h_roof;
    }
    rings_3d.push_back(r_roof);
    surf_type.push_back(2);
    mesh.push_polygon(r_roof, int(1));
    // mesh.push_attribute("surface_type", int(2));

    //walls
    size_t j_prev = ring.size()-1;
    for (size_t j=0; j<ring.size(); ++j) {
      LinearRing wall;
      wall.push_back(r_floor[j_prev]);
      wall.push_back(r_floor[j]);
      wall.push_back(r_roof[j]);
      wall.push_back(r_roof[j_prev]);

      rings_3d.push_back(wall);
      surf_type.push_back(1);
      mesh.push_polygon(wall, int(2));

      // mesh.push_attribute("surface_type", int(1));
      j_prev=j;
    }
    // walls from holes
    for (auto& hole : r_floor.interior_rings()) {
      size_t j_prev = hole.size()-1;
      for (size_t j=0; j<hole.size(); ++j) {
        LinearRing wall;
        wall.push_back(hole[j_prev]);
        wall.push_back(hole[j]);
        wall.push_back(arr3f{hole[j][0], hole[j][1], h_roof});
        wall.push_back(arr3f{hole[j_prev][0], hole[j_prev][1], h_roof});

        rings_3d.push_back(wall);
        surf_type.push_back(1);
        mesh.push_polygon(wall, int(2));
        j_prev=j;
      }
    }
  }
  
  //floor
  std::reverse(r_floor.begin(), r_floor.end());
  for (auto& hole : r_floor.interior_rings()) {
    std::reverse(hole.begin(), hole.end());
  }
  rings_3d.push_back(r_floor);
  surf_type.push_back(0);
  mesh.push_polygon(r_floor, int(0));
  // mesh.push_attribute("surface_type", int(0));
  std::unordered_map<int, Mesh> meshmap;
  meshmap[0] = mesh;

  output("surface_types").set(surf_type);
  output("mesh").set(meshmap);
}

inline arr3f grow(const arr3f& p_, const arr3f& q_, const arr3f& r_, const float& extension) {
  auto p = SCK::Point_2(p_[0], p_[1]);
  auto q = SCK::Point_2(q_[0], q_[1]);
  auto r = SCK::Point_2(r_[0], r_[1]);

  auto pq = q-p;
  pq = pq/pq.squared_length();
  auto pr = r-p;
  pr = pr/pr.squared_length();
  SCK::Vector_2 v;
  if (CGAL::collinear(p,q,r)) {
    v = SCK::Vector_2(pq.y(), -pq.x());
  } else {
    v = pq+pr;
    v = v/v.squared_length();
  }
  auto np = p + v*extension;
  return {float(np.x()), float(np.y()), 0};
}

void PolygonGrowerNode::process(){
  auto rings = input("rings").get<LinearRingCollection>();

  LinearRingCollection new_rings;
 
  for (auto& ring : rings) {
    size_t rsize = ring.size();
    LinearRing new_ring;
    for (size_t i=0; i<rsize; ++i) {
      auto& p = ring[i];
      auto& p_prev = ring[(i-1) % rsize];
      auto& p_next = ring[(i+1) % rsize];
      new_ring.push_back(
        grow(p, p_prev, p_next, extension)
      );
    }
    new_rings.push_back(new_ring);
  }

  output("rings").set(new_rings);
}

void Arr2LinearRingsNode::process() {
  

  // auto& floor_elevation = input("floor_elevation").get<float&>();
  // auto& mesh_error = input("mesh_error").get<float&>();
  // auto& roof_type = input("roof_type").get<int&>();
  // auto& arr_complexity = input("arr_complexity").get<int&>();
  auto& attributes_in = poly_input("attributes");
  auto& linear_rings = vector_output("linear_rings");

  //create all output fields
  std::unordered_map<std::string, gfSingleFeatureOutputTerminal*> input_attr_map;
  for (auto &iterm : attributes_in.sub_terminals()) {
    auto& oterm = poly_output("attributes").add_vector(iterm->get_name(), iterm->get_type());
    input_attr_map[oterm.get_name()] = &oterm;
  }
  // attributes specific to each roofpart
  auto &attr_elevation_50p_term = poly_output("attributes").add_vector("roof_elevation_50p", typeid(float));
  input_attr_map["roof_elevation_50p"] = &attr_elevation_50p_term;
  auto &attr_elevation_70p_term = poly_output("attributes").add_vector("roof_elevation_70p", typeid(float));
  input_attr_map["roof_elevation_70p"] = &attr_elevation_70p_term;
  auto &attr_elevation_min_term = poly_output("attributes").add_vector("roof_elevation_min", typeid(float));
  input_attr_map["roof_elevation_min"] = &attr_elevation_min_term;
  auto &attr_elevation_max_term = poly_output("attributes").add_vector("roof_elevation_max", typeid(float));
  input_attr_map["roof_elevation_max"] = &attr_elevation_max_term;
  auto &attr_datacov_term = poly_output("attributes").add_vector("data_coverage", typeid(float));
  input_attr_map["data_coverage"] = &attr_datacov_term;
  auto &attr_ground_term = poly_output("attributes").add_vector("is_ground", typeid(bool));
  input_attr_map["is_ground"] = &attr_ground_term;
  auto &part_id_term = poly_output("attributes").add_vector("building_part_id", typeid(int));
  input_attr_map["building_part_id"] = &part_id_term;
  


  // in case we get an invalid roof_type (-1 or -2) push the available attribute and an empty ring
  // if(invalid_rooftype) {
  //   for (auto &iterm : attributes_in.sub_terminals()) {
  //     // if (iterm->get_name()=="roof_type")
  //     if (iterm->has_data()) {
  //       input_attr_map[iterm->get_name()]->push_back_any(iterm->get_data());
  //     } else {
  //       input_attr_map[iterm->get_name()]->push_back_any(std::any());
  //     }
  //   }
  //   input_attr_map["roof_elevation_50p"]->push_back_any(std::any());
  //   input_attr_map["roof_elevation_70p"]->push_back_any(std::any());
  //   input_attr_map["roof_elevation_min"]->push_back_any(std::any());
  //   input_attr_map["roof_elevation_max"]->push_back_any(std::any());
  //   input_attr_map["data_coverage"]->push_back_any(std::any());
  //   input_attr_map["is_ground"]->push_back_any(std::any());
  //   input_attr_map["building_part_id"]->push_back_any(std::any());
  //   linear_rings.push_back_any(std::any());
  //   return;
  // } else {

    auto arr = input("arrangement").get<Arrangement_2>();

    // int j=0;
    // auto& plane_a = vector_output("plane_a");
    // auto& plane_b = vector_output("plane_b");
    // auto& plane_c = vector_output("plane_c");
    // auto& plane_d = vector_output("plane_d");
    for (auto face: arr.face_handles()) {
      bool skip_face = false;
      // if (only_in_footprint && !face->data().in_footprint) continue;
      if (output_groundparts && face->data().is_ground)
        skip_face = false;
      else 
        skip_face = !face->data().in_footprint;
      skip_face = skip_face || (face->is_fictitious() || face->is_unbounded());

      if (!skip_face) {
        LinearRing polygon;
        arrangementface_to_polygon(face, polygon);
        linear_rings.push_back(polygon);
        
        // plane paramters
        // plane_a.push_back(float(CGAL::to_double(face->data().plane.a())));
        // plane_b.push_back(float(CGAL::to_double(face->data().plane.b())));
        // plane_c.push_back(float(CGAL::to_double(face->data().plane.c())));
        // plane_d.push_back(float(CGAL::to_double(face->data().plane.d())));

        // attributes specific to each roofpart
        input_attr_map["roof_elevation_50p"]->push_back((float)face->data().elevation_50p);
        input_attr_map["roof_elevation_70p"]->push_back((float)face->data().elevation_70p);
        input_attr_map["roof_elevation_min"]->push_back((float)face->data().elevation_min);
        input_attr_map["roof_elevation_max"]->push_back((float)face->data().elevation_max);
        input_attr_map["data_coverage"]->push_back((float)face->data().data_coverage);
        input_attr_map["building_part_id"]->push_back((int)face->data().part_id);
        input_attr_map["is_ground"]->push_back((bool)face->data().is_ground);
        
        // input_attr_map["dak_id"]->push_back((int)++j);
        // input_attr_map["building_id"]->push_back(int(i+1));

        //attributes not specific to roofpart
        // input_attr_map["rmse"]->push_back((float)mesh_error);
        // input_attr_map["maaiveld_h"]->push_back((float)floor_elevation);
        // input_attr_map["dak_type"]->push_back((int)roof_type);
        // input_attr_map["arr_complexity"]->push_back(arr_complexity);
        for (auto &iterm : poly_input("attributes").sub_terminals()) {
          if (iterm->has_data()) {
            input_attr_map[iterm->get_name()]->push_back_any(iterm->get_data());
          } else {
            input_attr_map[iterm->get_name()]->push_back_any(std::any());
          }
        }
      }
    }
    auto& groundparts = vector_input("groundparts");
    for(size_t i=0; i<groundparts.size(); ++i) {
      auto lr = groundparts.get<LinearRing>(i);
      linear_rings.push_back(lr);

      input_attr_map["roof_elevation_50p"]->push_back_any(std::any());
      input_attr_map["roof_elevation_70p"]->push_back_any(std::any());
      input_attr_map["roof_elevation_min"]->push_back_any(std::any());
      input_attr_map["roof_elevation_max"]->push_back_any(std::any());
      input_attr_map["data_coverage"]->push_back_any(std::any());
      input_attr_map["building_part_id"]->push_back_any(std::any());
      input_attr_map["is_ground"]->push_back(true);
      
      for (auto &iterm : poly_input("attributes").sub_terminals()) {
        if (iterm->has_data()) {
          input_attr_map[iterm->get_name()]->push_back_any(iterm->get_data());
        } else {
          input_attr_map[iterm->get_name()]->push_back_any(std::any());
        }
      }
    }
  // }
}

void Arr2LinearRingsDebugNode::process() {
  

  // auto& floor_elevation = input("floor_elevation").get<float&>();
  // auto& mesh_error = input("mesh_error").get<float&>();
  // auto& roof_type = input("roof_type").get<int&>();
  // auto& arr_complexity = input("arr_complexity").get<int&>();
  auto& linear_rings = vector_output("linear_rings");

  //create all output fields
  std::unordered_map<std::string, gfSingleFeatureOutputTerminal*> input_attr_map;
  // attributes specific to each roofpart
  input_attr_map["roof_elevation_50p"] = 
  & poly_output("attributes").add_vector("roof_elevation_50p", typeid(float));

  input_attr_map["roof_elevation_70p"] = 
  & poly_output("attributes").add_vector("roof_elevation_70p", typeid(float));

  input_attr_map["roof_elevation_min"] = 
  & poly_output("attributes").add_vector("roof_elevation_min", typeid(float));

  input_attr_map["roof_elevation_max"] = 
  & poly_output("attributes").add_vector("roof_elevation_max", typeid(float));

  input_attr_map["data_coverage"] = 
  & poly_output("attributes").add_vector("data_coverage", typeid(float));

  input_attr_map["is_ground"] = 
  & poly_output("attributes").add_vector("is_ground", typeid(bool));

  input_attr_map["in_footprint"] = 
  & poly_output("attributes").add_vector("in_footprint", typeid(bool));

  input_attr_map["is_footprint_hole"] = 
  & poly_output("attributes").add_vector("is_footprint_hole", typeid(bool));

  input_attr_map["building_part_id"] = 
  & poly_output("attributes").add_vector("building_part_id", typeid(int));
  
  input_attr_map["label"] = 
  & poly_output("attributes").add_vector("label", typeid(int));

  input_attr_map["pixel_count"] = 
  & poly_output("attributes").add_vector("pixel_count", typeid(int));

  

    auto arr = input("arrangement").get<Arrangement_2>();

    // int j=0;
    // auto& plane_a = vector_output("plane_a");
    // auto& plane_b = vector_output("plane_b");
    // auto& plane_c = vector_output("plane_c");
    // auto& plane_d = vector_output("plane_d");
    for (auto face: arr.face_handles()) {
      if (face->is_fictitious() || face->is_unbounded()) continue;

      LinearRing polygon;
      arrangementface_to_polygon(face, polygon);
      linear_rings.push_back(polygon);
      
      // plane paramters
      // plane_a.push_back(float(CGAL::to_double(face->data().plane.a())));
      // plane_b.push_back(float(CGAL::to_double(face->data().plane.b())));
      // plane_c.push_back(float(CGAL::to_double(face->data().plane.c())));
      // plane_d.push_back(float(CGAL::to_double(face->data().plane.d())));

      // attributes specific to each roofpart
      input_attr_map["roof_elevation_50p"]->push_back((float)face->data().elevation_50p);
      input_attr_map["roof_elevation_70p"]->push_back((float)face->data().elevation_70p);
      input_attr_map["roof_elevation_min"]->push_back((float)face->data().elevation_min);
      input_attr_map["roof_elevation_max"]->push_back((float)face->data().elevation_max);
      input_attr_map["data_coverage"]->push_back((float)face->data().data_coverage);
      input_attr_map["building_part_id"]->push_back((int)face->data().part_id);
      input_attr_map["label"]->push_back((int)face->data().label);
      input_attr_map["pixel_count"]->push_back((int)face->data().pixel_count);
      input_attr_map["is_ground"]->push_back((bool)face->data().is_ground);
      input_attr_map["in_footprint"]->push_back((bool)face->data().in_footprint);
      input_attr_map["is_footprint_hole"]->push_back((bool)face->data().is_footprint_hole);
    }

    PointCollection isolated_vertices;
    for (auto vertex: arr.vertex_handles()) {
      if (vertex->is_isolated()) {
        isolated_vertices.push_back({
          float(CGAL::to_double(vertex->point().x())),
          float(CGAL::to_double(vertex->point().y())),
          0
        });
      }
    }
    output("isolated_vertices").set(isolated_vertices);
  // }
}

inline arr3f v2p(Arrangement_2::Vertex_handle v, float h) {
  return {
          float(CGAL::to_double(v->point().x())),
          float(CGAL::to_double(v->point().y())),
          h
        };
}
template<typename P> arr3f p2p(P p) {
  return {
          float(CGAL::to_double(p->x())),
          float(CGAL::to_double(p->y())),
          float(CGAL::to_double(p->z()))
        };
}

typedef std::pair<float, Face_handle> hf_pair;
vec3f get_heights(std::vector<hf_pair>& vertex_column, Vertex_handle v, Face_handle f_a, Face_handle f_b, float& h_a, float& h_b) {
  vec3f v_other;
  float h_prev = -999999;
  bool found_first = false;
  for(auto& [h,face] : vertex_column) {
    if (face==f_a) {
      h_a = h;
      if (!found_first) 
        found_first = true;
      else {
        if (h==h_prev && v_other.size()) {
          v_other.erase(v_other.end()-1);
        }
        break;
      }
    } else if (face==f_b) {
      h_b = h;
      if (!found_first) 
        found_first = true;
      else {
        if (h==h_prev && v_other.size()) {
          v_other.erase(v_other.end()-1);
        }
        break;
      }
    } else if (found_first) {
      if (h!=h_prev)
        v_other.push_back(v2p(v, h));
    } //else break;
    h_prev = h;
  }
  return v_other;
}

template<typename T> void push_ccb(
  T& ring, 
  Halfedge_handle hedge, 
  std::unordered_map<Vertex_handle, std::vector<hf_pair>>& vertex_columns, 
  std::unordered_map<Halfedge_handle, 
  AK::Point_3>& extra_wall_points, 
  float& snap_tolerance) {

  auto first = hedge;
  do {
    auto v = hedge->source();
    if(CGAL::squared_distance(v->point(), hedge->target()->point()) > snap_tolerance) {
      for(auto& [h,f_h] : vertex_columns[v]) {
        if (f_h==hedge->face()) {
          ring.push_back(v2p(v,h));
          break;
        }
      }
      auto p_xtra = extra_wall_points.find(hedge);
      auto q_xtra = extra_wall_points.find(hedge->twin());
      if (p_xtra != extra_wall_points.end()) {
        ring.push_back(p2p(&p_xtra->second));
      } else if (q_xtra != extra_wall_points.end()) {
        ring.push_back(p2p(&q_xtra->second));
      }
    }
    hedge = hedge->next();
  } while (hedge!=first);
}

void ArrExtruderNode::process(){
  typedef Arrangement_2::Traits_2 AT;
  auto arr = input("arrangement").get<Arrangement_2>();
  auto floor_elevation = input("floor_elevation").get<float>();
  float snap_tolerance = std::pow(10,-snap_tolerance_exp);

  // assume we have only one unbounded faces that just has the building footprint as a hole

  // std::vector<LinearRing> faces;
  auto& faces = vector_output("faces");
  auto& surface_labels = vector_output("labels");

  std::unordered_map<int, Mesh> meshmap;

  auto unbounded_face = arr.unbounded_face();
  unbounded_face->data().elevation_70p=floor_elevation;

  // floor
  if (do_floor) {
    // std::cout << "arrangement has " << arr.number_of_unbounded_faces() << "unbounded faces\n";
    // there should only be one hole in the unbounded face (building consists of one part)
    for(Arrangement_2::Hole_iterator floorpart_ = unbounded_face->holes_begin(); floorpart_ != unbounded_face->holes_end(); ++floorpart_ ) {
      LinearRing floor;
      auto he = *floorpart_;
      auto first = he;
      do {
        if(CGAL::squared_distance(he->source()->point(), he->target()->point()) > snap_tolerance)
          floor.push_back(v2p(he->source(), floor_elevation));
        he = he->next();
      } while(he!=first);

      faces.push_back(floor);
      surface_labels.push_back(int(0));

      // create a mesh
      Mesh mesh;
      mesh.push_polygon(floor, int(0));
      meshmap[he->twin()->face()->data().part_id] = mesh;
      // mesh.push_attribute("surface_type", int(0));
    }
    // get the holes
    for (auto face: arr.face_handles()) {
      if (face->data().is_footprint_hole) {
        vec3f hole;
        auto he = face->outer_ccb();
        auto first = he;
        do {
          if(CGAL::squared_distance(he->source()->point(), he->target()->point()) > snap_tolerance)
            hole.push_back(v2p(he->source(), floor_elevation));
          he = he->next();
        } while(he!=first);

        // attach hole to appropriate mesh and linear ring
        const auto part_id = he->twin()->face()->data().part_id;
        if(hole.size()!=0 && meshmap.find(part_id)!=meshmap.end()) {
          auto& floor = meshmap[part_id].get_polygons()[0];
          floor.interior_rings().push_back(hole);
        } else { std::cout << "skipping a hole for which no polygon exists...\n"; }
      }
    }
  }



  // check for no detected planes (if one face inside fp has none then all have none)
  for (const auto f : arr.face_handles()) {
    if (f->data().in_footprint && f->data().segid==0) {
      vector_output("mesh").push_back(meshmap[0]);
      return;
    }
  }
  

  // compute all heights for each vertex
  std::unordered_map<Vertex_handle, std::vector<hf_pair>> vertex_columns;
  for(auto& v : arr.vertex_handles()) {
    auto& p = v->point();
    auto he = v->incident_halfedges();
    auto first = he;
    std::vector<hf_pair> heights;
    do {
      auto f = he->face();
      float h;
      if (f->data().in_footprint && f->data().segid!=0) {
        if(LoD2) {
          auto& plane = f->data().plane;
          h = (plane.a()*CGAL::to_double(p.x()) + plane.b()*CGAL::to_double(p.y()) + plane.d()) / (-plane.c());
        } else {
          h = f->data().elevation_70p;
        }
      } else if (f->data().in_footprint) {
        h = nodata_elevation;
      } else {
        h = floor_elevation;
      }
      heights.push_back(std::make_pair(h, f));
    } while (++he!=first);
    
    if(heights.size()==0) continue;
    
    // sort heights
    std::sort(heights.begin(), heights.end(), [](hf_pair& a, hf_pair& b) {
      return a.first < b.first;   
    });
    // equalise heights that are practically the same
    float h_ref = heights[0].first;
    for (auto& [h,face] : heights) {
      if(std::fabs(h_ref-h)<snap_tolerance) {
        h = h_ref;
      }
      h_ref = h;
    }

    vertex_columns[v]=heights;
  }

  // walls
  // store points that need to be created to do the walls right. We need them later for the roofs
  std::unordered_map<Halfedge_handle, AK::Point_3> extra_wall_points;
  if (do_walls) {
    for (auto edge : arr.edge_handles()) {
      auto e_a = edge->twin();
      auto e_b = edge;
      auto v1 = e_a->target();
      auto v2 = e_a->source();
      auto& p1 = v1->point();
      auto& p2 = v2->point();
      auto f_a = e_a->face();
      auto f_b = e_b->face();
      bool fp_a = f_a->data().in_footprint;
      bool fp_b = f_b->data().in_footprint;
      auto& pl_a = f_a->data().plane;
      auto& pl_b = f_b->data().plane;

      // edge has no length
      if(CGAL::squared_distance(p1,p2)<snap_tolerance)
        continue;

      // edge is not in footprint nor on its boundary
      if (!fp_a & !fp_b) {
        continue;
      }

      // get precomputed heights from vertex column
      float h1a, h1b, h2a, h2b;
      vec3f v1_other = get_heights(vertex_columns[v1], v1, f_a, f_b, h1a, h1b);
      vec3f v2_other = get_heights(vertex_columns[v2], v2, f_a, f_b, h2a, h2b);
      
      // // set base (ground) elevation to vertices adjacent to a face oustide the building fp
      int part_id;
      if (fp_a && !fp_b) {
        h1b=h2b=floor_elevation;
        part_id = f_a->data().part_id;
      } else if (!fp_a && fp_b) {
        h1a=h2a=floor_elevation;
        part_id = f_b->data().part_id;
      } else{ // both sides have the same part_id
        part_id = f_b->data().part_id; 
      }

      auto& mesh = meshmap[part_id];

      int wall_label = 3; //inner wall
      if (!fp_a || !fp_b) wall_label = 2; //outer wall

      LinearRing wall_face_1;
      if ((h1a<h1b) and (h2a<h2b)) {
        wall_face_1.push_back(v2p(v1,h1b));
        // v1_other desc
        wall_face_1.insert(wall_face_1.end(), v1_other.rbegin(), v1_other.rend());
        wall_face_1.push_back(v2p(v1,h1a));
        wall_face_1.push_back(v2p(v2,h2a));
        // v2_other asc
        wall_face_1.insert(wall_face_1.end(), v2_other.begin(), v2_other.end());
        wall_face_1.push_back(v2p(v2,h2b));
        faces.push_back(wall_face_1);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_1, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
      } else 
      if ((h1a>h1b) and (h2a>h2b)) {
        wall_face_1.push_back(v2p(v2,h2a));
        // v2_other desc
        wall_face_1.insert(wall_face_1.end(), v2_other.rbegin(), v2_other.rend());
        wall_face_1.push_back(v2p(v2,h2b));
        wall_face_1.push_back(v2p(v1,h1b));
        // v1_other asc
        wall_face_1.insert(wall_face_1.end(), v1_other.begin(), v1_other.end());
        wall_face_1.push_back(v2p(v1,h1a));
        faces.push_back(wall_face_1);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_1, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
      } else 
      if ((h1a<h1b) and (h2a>h2b)) {
        // compute vx and hx
        auto l_a = AK::Line_3(AK::Point_3(p1.x(), p1.y(), h1a), AK::Point_3(p2.x(), p2.y(), h2a));
        auto l_b = AK::Line_3(AK::Point_3(p1.x(), p1.y(), h1b), AK::Point_3(p2.x(), p2.y(), h2b));
        auto result = CGAL::intersection(l_a, l_b);
        auto px = boost::get<typename AK::Point_3>(&*result);

        // TODO: check if distance from px to p1 and p2 is longer than snap_tolerance?

        extra_wall_points[edge] = *px;
        // AK::Point_2 px_2d(px->x(),px->y());
        // arr.split_edge(edge, AT::Segment_2(p1,px_2d), AT::Segment_2(px_2d,p2));
        // if (result) {
          // auto px = )
        // }

        wall_face_1.push_back(v2p(v1,h1b));
        // v1_other desc
        wall_face_1.insert(wall_face_1.end(), v1_other.rbegin(), v1_other.rend());
        wall_face_1.push_back(v2p(v1,h1a));
        wall_face_1.push_back(p2p(px));

        LinearRing wall_face_2;
        wall_face_2.push_back(v2p(v2,h2a));
        // v2_other desc
        wall_face_2.insert(wall_face_2.end(), v2_other.rbegin(), v2_other.rend());
        wall_face_2.push_back(v2p(v2,h2b));
        wall_face_2.push_back(p2p(px));
        
        faces.push_back(wall_face_1);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_1, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
        faces.push_back(wall_face_2);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_2, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
      } else 
      if ((h1a>h1b) and (h2a<h2b)) {
        // compute vx and hx
        auto l_a = AK::Line_3(AK::Point_3(p1.x(), p1.y(), h1a), AK::Point_3(p2.x(), p2.y(), h2a));
        auto l_b = AK::Line_3(AK::Point_3(p1.x(), p1.y(), h1b), AK::Point_3(p2.x(), p2.y(), h2b));
        auto result = CGAL::intersection(l_a, l_b);
        auto px = boost::get<typename AK::Point_3>(&*result);

        // TODO: check if distance from px to p1 and p2 is longer than snap_tolerance?

        extra_wall_points[edge] = *px;
        // AK::Point_2 px_2d(px->x(),px->y());
        // arr.split_edge(edge, AT::Segment_2(p1,px_2d), AT::Segment_2(px_2d,p2));
        // if (result) {
          // auto px = )
        // }

        wall_face_1.push_back(v2p(v1,h1b));
        // v1_other asc
        wall_face_1.insert(wall_face_1.end(), v1_other.begin(), v1_other.end());
        wall_face_1.push_back(v2p(v1,h1a));
        wall_face_1.push_back(p2p(px));

        LinearRing wall_face_2;
        wall_face_2.push_back(v2p(v2,h2a));
        // v2_other asc
        wall_face_2.insert(wall_face_2.end(), v2_other.begin(), v2_other.end());
        wall_face_2.push_back(v2p(v2,h2b));
        wall_face_2.push_back(p2p(px));
        
        faces.push_back(wall_face_1);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_1, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
        faces.push_back(wall_face_2);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_2, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
      } else 
      if ((h1a>h1b) and (h2a==h2b)) {
        wall_face_1.push_back(v2p(v1,h1b));
        // v1_other asc
        wall_face_1.insert(wall_face_1.end(), v1_other.begin(), v1_other.end());
        wall_face_1.push_back(v2p(v1,h1a));
        wall_face_1.push_back(v2p(v2,h2a));
        faces.push_back(wall_face_1);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_1, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
      } else 
      if ((h1a<h1b) and (h2a==h2b)) {
        wall_face_1.push_back(v2p(v1,h1b));
        // v1_other desc
        wall_face_1.insert(wall_face_1.end(), v1_other.rbegin(), v1_other.rend());
        wall_face_1.push_back(v2p(v1,h1a));
        wall_face_1.push_back(v2p(v2,h2a));
        faces.push_back(wall_face_1);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_1, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
      } else 
      if ((h2b>h2a) and (h1a==h1b)) {
        wall_face_1.push_back(v2p(v2,h2a));
        // v2_other asc
        wall_face_1.insert(wall_face_1.end(), v2_other.begin(), v2_other.end());
        wall_face_1.push_back(v2p(v2,h2b));
        wall_face_1.push_back(v2p(v1,h1a));
        faces.push_back(wall_face_1);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_1, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
      } else 
      if ((h2b<h2a) and (h1a==h1b)) {
        wall_face_1.push_back(v2p(v2,h2a));
        // v2_other desc
        wall_face_1.insert(wall_face_1.end(), v2_other.rbegin(), v2_other.rend());
        wall_face_1.push_back(v2p(v2,h2b));
        wall_face_1.push_back(v2p(v1,h1a));
        faces.push_back(wall_face_1);
        surface_labels.push_back(wall_label);
        mesh.push_polygon(wall_face_1, wall_label);
        // mesh.push_attribute("surface_type", wall_label);
      }
    }
  }

  // roofs
  if (do_roofs) {
    for (auto face: arr.face_handles()) {
      if (face->data().in_footprint) {
        LinearRing roofpart;
        auto he = face->outer_ccb();
        
        push_ccb(roofpart, he, vertex_columns, extra_wall_points, snap_tolerance);

        for(Arrangement_2::Hole_iterator hole = face->holes_begin(); hole != face->holes_end(); ++hole ) {
          vec3f roofpart_hole;
          auto he = *hole;
          push_ccb(roofpart_hole, he, vertex_columns, extra_wall_points, snap_tolerance);
          if (roofpart_hole.size()>2) {
            roofpart.interior_rings().push_back(roofpart_hole);
          }
        }

        if (roofpart.size()>2) {
          faces.push_back(roofpart);
          surface_labels.push_back(int(1));
          meshmap[face->data().part_id].push_polygon(roofpart, int(1));
          // mesh.push_attribute("surface_type", int(1));
        }

      }
    }

  }

  for (size_t i=0; i< meshmap.size(); ++i) {
    vector_output("mesh").push_back(meshmap[i]);
  }
  output("multisolid").push_back(meshmap);
  
}

// void ExtruderNode::process(){
//   // if (!(do_walls || do_roofs)) return;
//   // Set up vertex data (and buffer(s)) and attribute pointers
//   // auto polygons = std::any_cast<input("polygons").get<vec2f>>();
//   // auto elevations = std::any_cast<input("elevations").get<float>>();
//   auto arr = input("arrangement").get<Arrangement_2>();

//   TriangleCollection triangles;
//   vec3f normals;
//   vec1i cell_id_vec1i, plane_id, face_ids;
//   vec1i labels;
//   vec1f rms_errors, max_errors, segment_coverages, elevations;
//   using N = uint32_t;
  
//   size_t cell_id=0, pid, face_id=0; // face_id==0 is the floor
//   float rms_error, max_error, segment_coverage;
//   for (auto face: arr.face_handles()) {
//     // bool extract = param<bool>("in_footprint") ? face->data().in_footprint : face->data().is_finite;
//     bool extract = face->data().in_footprint;
//     if(extract) {
//       cell_id++;
//       rms_error = face->data().rms_error_to_avg;
//       max_error = face->data().max_error;
//       segment_coverage = face->data().segid_coverage;
//       pid = face->data().segid;
//       if (pid==0) {
//         face->data().plane = Plane(Kernel::Point_3(0,0,nodata_elevation), Kernel::Direction_3(0,0,1));
//         face->data().elevation_avg = nodata_elevation;
//       }

//       vec2f polygon, vertices;
//       arrangementface_to_polygon(face, polygon);
//       std::vector<N> indices = mapbox::earcut<N>(std::vector<vec2f>({polygon}));
//       for(auto i : indices) {
//         vertices.push_back({polygon[i]});
//       }
//       // floor triangles
//       for (size_t i=0; i<indices.size()/3; ++i) {
//         Triangle triangle;
//         for (size_t j=0; j<3; ++j) {
//           triangle[j] = {vertices[i*3+j][0], vertices[i*3+j][1], base_elevation};
//           labels.push_back(0);
//           normals.push_back({0,0,-1});
//           cell_id_vec1i.push_back(cell_id);
//           plane_id.push_back(pid);
//           face_ids.push_back(0);
//           rms_errors.push_back(rms_error);
//           elevations.push_back(face->data().elevation_avg);
//           max_errors.push_back(max_error);
//           segment_coverages.push_back(segment_coverage);
//         }
//         triangles.push_back(triangle);
//       }
//       if (do_roofs) {
//         // roof triangles
//         ++face_id;
//         auto& plane = face->data().plane;
//         arr3f roof_normal = {0,0,1};
//         if (LoD2) {
//           auto n_ = plane.orthogonal_vector();
//           if(n_*Vector(0.0f,0.0f,1.0f) <0) n_*=-1;
//           n_ = n_/CGAL::sqrt(n_.squared_length());
//           roof_normal = {float(n_.x()), float(n_.y()), float(n_.z())};
//         }
//         for (size_t i=0; i<indices.size()/3; ++i) {
//           Triangle triangle;
//           for (size_t j=0; j<3; ++j) {
//             auto& px = vertices[i*3+j][0];
//             auto& py = vertices[i*3+j][1];
//             float h;
//             if (LoD2) {
//               h = (plane.a()*px + plane.b()*py + plane.d()) / (-plane.c());
//             } else {
//               h = face->data().elevation_avg;
//             }
//             triangle[j] = {px, py, h};
//             labels.push_back(1);
//             normals.push_back(roof_normal);
//             cell_id_vec1i.push_back(cell_id);
//             face_ids.push_back(face_id);
//             plane_id.push_back(pid);
//             // rms_errors.push_back(rms_error);
//             rms_errors.push_back(rms_error);
//             elevations.push_back(face->data().elevation_avg);
//             max_errors.push_back(max_error);
//             segment_coverages.push_back(segment_coverage);
//           }
//           triangles.push_back(triangle);
//         }
//       }
//     }
//   }
//   if (do_walls) {
//     vertex n;
//     for (auto edge : arr.edge_handles()) {
//       // skip if faces on both sides of this edge are not finite
//       bool fp_u = edge->twin()->face()->data().in_footprint;
//       bool fp_l = edge->face()->data().in_footprint;
//       auto plane_u = edge->twin()->face()->data().plane;
//       auto plane_l = edge->face()->data().plane;
//       if (fp_u || fp_l) {
//         ++face_id;
//         auto& source = edge->source()->point();
//         auto& target = edge->target()->point();
//         int wall_label = 2;
//         if (fp_u && fp_l)
//           wall_label = 3;

//         float u1z,u2z,l1z,l2z;
//         if(LoD2){
//           u1z = 
//           (plane_u.a()*CGAL::to_double(source.x()) + plane_u.b()*CGAL::to_double(source.y()) + plane_u.d()) / (-plane_u.c());
//           u2z = 
//           (plane_u.a()*CGAL::to_double(target.x()) + plane_u.b()*CGAL::to_double(target.y()) + plane_u.d()) / (-plane_u.c());
//           l1z = 
//           (plane_l.a()*CGAL::to_double(source.x()) + plane_l.b()*CGAL::to_double(source.y()) + plane_l.d()) / (-plane_l.c());
//           l2z = 
//           (plane_l.a()*CGAL::to_double(target.x()) + plane_l.b()*CGAL::to_double(target.y()) + plane_l.d()) / (-plane_l.c());
//         } else {
//           l1z = l2z = edge->face()->data().elevation_avg;
//           u1z = u2z = edge->twin()->face()->data().elevation_avg;
//         }
//         // set base (ground) elevation to vertices adjacent to a face oustide the building fp
//         if (fp_u && !fp_l) l1z=l2z=base_elevation;
//         if (!fp_u && fp_l) u1z=u2z=base_elevation;
//         // push 2 triangles to form the quad between lower and upper edges
//         // notice that this is not always topologically correct, but fine for visualisation
        
//         // define four points of the quad between upper and lower edge
//         std::array<float,3> l1,l2,u1,u2;
//         l1 = {
//           float(CGAL::to_double(source.x())),
//           float(CGAL::to_double(source.y())),
//           l1z
//         };
//         l2 = {
//           float(CGAL::to_double(target.x())),
//           float(CGAL::to_double(target.y())),
//           l2z
//         };
//         u1 = {
//           float(CGAL::to_double(source.x())),
//           float(CGAL::to_double(source.y())),
//           u1z
//         };
//         u2 = {
//           float(CGAL::to_double(target.x())),
//           float(CGAL::to_double(target.y())),
//           u2z
//         };

//         // 1st triangle
//         triangles.push_back({u1,l2,l1});
//         labels.push_back(wall_label);
//         // triangles.push_back(l2);
//         labels.push_back(wall_label);
//         // triangles.push_back(l1);
//         labels.push_back(wall_label);
//         face_ids.push_back(face_id);
//         face_ids.push_back(face_id);
//         face_ids.push_back(face_id);

//         n = get_normal(u1,l2,l1);
//         normals.push_back(n);
//         normals.push_back(n);
//         normals.push_back(n);

//         cell_id_vec1i.push_back(0);
//         plane_id.push_back(0);
//         cell_id_vec1i.push_back(0);
//         plane_id.push_back(0);
//         cell_id_vec1i.push_back(0);
//         plane_id.push_back(0);
//         rms_errors.push_back(-1);
//         elevations.push_back(-1);
//         rms_errors.push_back(-1);
//         elevations.push_back(-1);
//         rms_errors.push_back(-1);
//         elevations.push_back(-1);
//         max_errors.push_back(-1);
//         max_errors.push_back(-1);
//         max_errors.push_back(-1);
//         segment_coverages.push_back(-1);
//         segment_coverages.push_back(-1);
//         segment_coverages.push_back(-1);

//         // 2nd triangle
//         triangles.push_back({u1,u2,l2});
//         labels.push_back(wall_label);
//         // triangles.push_back(u2);
//         labels.push_back(wall_label);
//         // triangles.push_back(l2);
//         labels.push_back(wall_label);
//         face_ids.push_back(face_id);
//         face_ids.push_back(face_id);
//         face_ids.push_back(face_id);

//         n = get_normal(u1,u2,l2);
//         normals.push_back(n);
//         normals.push_back(n);
//         normals.push_back(n);

//         cell_id_vec1i.push_back(0);
//         plane_id.push_back(0);
//         cell_id_vec1i.push_back(0);
//         plane_id.push_back(0);
//         cell_id_vec1i.push_back(0);
//         plane_id.push_back(0);
//         rms_errors.push_back(-1);
//         elevations.push_back(-1);
//         rms_errors.push_back(-1);
//         elevations.push_back(-1);
//         rms_errors.push_back(-1);
//         elevations.push_back(-1);
//         max_errors.push_back(-1);
//         max_errors.push_back(-1);
//         max_errors.push_back(-1);
//         segment_coverages.push_back(-1);
//         segment_coverages.push_back(-1);
//         segment_coverages.push_back(-1);
//       }
//     }
//   }
  
//   output("normals_vec3f").set(normals);
//   output("cell_id_vec1i").set(cell_id_vec1i);
//   output("plane_id").set(plane_id);
//   output("rms_errors").set(rms_errors);
//   output("max_errors").set(max_errors);
//   output("elevations").set(elevations);
//   output("segment_coverages").set(segment_coverages);
//   output("triangles").set(triangles);
//   output("labels_vec1i").set(labels);
//   output("face_ids").set(face_ids);
// }

void LinearRingtoRingsNode::process(){
  auto lr = input("linear_ring").get<LinearRing>();
  LinearRingCollection lrc;
  lrc.push_back(lr);
  output("linear_rings").set(lrc);
}

inline void DetectLinesNode::detect_lines_ring_m1(linedect::LineDetector& LD, SegmentCollection& segments_out) {
  LD.dist_thres = dist_thres * dist_thres;
  LD.N = k;
  std::vector<size_t> detected_regions;
  size_t ringsize = LD.point_segment_idx.size();
  RingSegMap ring_seg_map;
  for (size_t i=min_cnt_range_upper; i>=min_cnt_range_lower; --i){
    LD.min_segment_count = i;
    auto new_regions = LD.detect();

    // extract the new regions from the rring
    std::vector<size_t> reset_ids;
    for (const auto& cur_reg : new_regions) {
      bool stitch = (cur_reg != LD.point_segment_idx.back()) && (cur_reg != LD.point_segment_idx[0]);

      std::vector<std::vector<size_t>> other_regs(1);
      size_t j=0;
      for (auto& regid : LD.point_segment_idx) {
        if (regid != cur_reg) {
          other_regs.back().push_back(j);
        } else {
          other_regs.resize(other_regs.size()+1);
        }
        ++j;
      }
      if (stitch) {
        auto& front = other_regs.front();
        auto& back = other_regs.back();
        front.insert(front.begin(), back.begin(), back.end());
        other_regs.pop_back();
      }
      size_t largest = 0, largest_id;
      j=0;
      for (const auto& other_reg : other_regs) {
        auto s = other_reg.size();
        if (s >largest) {
          largest=s;
          largest_id=j;
        } ++j;
      }
      auto idend = other_regs[largest_id].front()-1 % ringsize;
      auto idstart = other_regs[largest_id].back()+1 % ringsize;
      ring_seg_map[std::make_pair(idstart,idend)] = cur_reg;
      reset_ids.push_back(idend);
      reset_ids.push_back(idstart);
    }
    for (const auto& j : reset_ids) {
      LD.point_segment_idx[j] = 0;
    }
  }
  for (auto& [seg,rid] : ring_seg_map) {
    segments_out.push_back(LD.project(seg.first, seg.second));
  }
}
inline void DetectLinesNode::detect_lines(linedect::LineDetector& LD) {
  LD.dist_thres = dist_thres * dist_thres;
  LD.N = k;
  auto& c_upper = min_cnt_range.second;
  auto& c_lower = min_cnt_range.first;
  for (size_t i=c_upper; i>=c_lower; --i){
    LD.min_segment_count = i;
    LD.detect();
  }
}
inline size_t DetectLinesNode::detect_lines_ring_m2(linedect::LineDetector& LD, Plane& plane, SegmentCollection& segments_out) {
  LD.dist_thres = dist_thres * dist_thres;
  LD.N = k;
  auto& c_upper = min_cnt_range.second;
  auto& c_lower = min_cnt_range.first;
  for (size_t i=c_upper; i>=c_lower; --i){
    LD.min_segment_count = i;
    LD.detect();
  }
  size_t ringsize = LD.point_segment_idx.size();
            // chain the detected lines, to ensure correct order
  if (LD.segment_shapes.size()>1) {
    std::vector<std::pair<size_t,size_t>> new_ring_ids;
    bool start_seg = LD.point_segment_idx[0];
    int prev_i=ringsize-1,
      prev_seg=LD.point_segment_idx[prev_i], 
      cur_seg, 
      i_last_seg = -1;
    bool perfect_aligned=false; // is the first point of the ring also the first point of a segment? If yes, we are perfectly aligned!
    for( int i=0; i<ringsize; ++i ) {
      cur_seg = LD.point_segment_idx[i];
      if(cur_seg==prev_seg && cur_seg!=0) { // we are inside a segment

      } else if (cur_seg!=0 && prev_seg==0) { // from unsegmented to segmented
        new_ring_ids.push_back(std::make_pair(i, cur_seg)); // first of cur
        if(i==0) perfect_aligned=true;
        // new_ring_ids.push_back(i); // end of unsegmented linesegment
      } else if (cur_seg!=0 && prev_seg!=0) { // from one segment to another
        new_ring_ids.push_back(std::make_pair(prev_i, prev_seg)); // last of prev
        new_ring_ids.push_back(std::make_pair(i, cur_seg)); // first of cur
      } else if (cur_seg==0 && prev_seg!=0) { //from segment to unsegmented
        new_ring_ids.push_back(std::make_pair(prev_i, prev_seg)); // last of prev
        // new_ring_ids.push_back(prev_i); // begin of unsegmented linesegment
      } // else: we are inside an unsegmented or segmented zone
      prev_seg = cur_seg;
      prev_i = i;
    }
    if (!perfect_aligned) { // add the segment that runs through the first point in the original ring
      new_ring_ids.insert(new_ring_ids.begin(), new_ring_ids.back());
      new_ring_ids.pop_back();
    }
    //ensure the ring is aligned wrt diff region ids around origin of the ring
    if (new_ring_ids.front().second == new_ring_ids.back().second) {
      size_t region = new_ring_ids.front().second;
      do {
        new_ring_ids.push_back(new_ring_ids.front());
        new_ring_ids.erase(new_ring_ids.begin());
      } while (region == new_ring_ids.front().second);
    }
    //merge multiple segments of the same region
    std::unordered_map<size_t,std::pair<size_t,size_t>> map_by_region;
    for (auto el = new_ring_ids.begin(); el<new_ring_ids.end(); ++el) {
      if(!map_by_region.count(el->second)) {
        map_by_region[el->second] = std::make_pair(el->first, el->first);
      } else {
        map_by_region[el->second].second = el->first;
      }
    }
    //sort the segments acc to order in ring
    typedef std::set<std::pair<size_t,size_t>,Cmp> SegSet;
    SegSet sorted_segments;
    for (auto& el : map_by_region) {
      sorted_segments.insert(el.second);
    }
    // TODO: better check for overlapping segments! Following is not 100% robust...
    if (remove_overlap) {
      auto el_prev = sorted_segments.begin();
      // el_prev.first-=sorted_segments.size();
      // el_prev.second-=sorted_segments.size();
      std::vector<SegSet::iterator> to_remove;
      for (auto el = ++sorted_segments.begin(); el != sorted_segments.end(); ++el ){
        el_prev = el;
        --el_prev;
        if (el_prev->second > el->first)
          to_remove.push_back(el);
      }
      for (auto el : to_remove) {
        sorted_segments.erase(el);
      }
    }
    if (perform_chaining) {
      std::vector<SCK::Segment_3> prechain_segments;
      std::vector<size_t> idx; size_t idcnt=0;
      for (auto& [i0,i1] : sorted_segments) {
        // segments_out.push_back( LD.project(i0, i1) );
        prechain_segments.push_back( LD.project_cgal(i0, i1, line_extend) );
        idx.push_back(idcnt++);
      }
      // TODO: chain the ring? for better regularisation results
      SegmentCollection new_ring;
      auto chained_ring_pts = linereg::chain_ring<SCK>(idx, SCK::Plane_3(plane.a(), plane.b(), plane.c(), plane.d()), prechain_segments, snap_threshold, line_extend);

      if (chained_ring_pts.size() > 2) {
        auto first = chained_ring_pts.begin();
        for (auto pit=std::next(first); pit!=chained_ring_pts.end(); ++pit) {
          auto p2 = *pit;
          auto p1 = *std::prev(pit);
          segments_out.push_back({
            arr3f{
              float(CGAL::to_double(p1.x())),
              float(CGAL::to_double(p1.y())),
              float(CGAL::to_double(p1.z()))},
            arr3f{
              float(CGAL::to_double(p2.x())),
              float(CGAL::to_double(p2.y())),
              float(CGAL::to_double(p2.z()))},
          });
        }
        auto p1 = *chained_ring_pts.rbegin();
        auto p2 = *first;
        segments_out.push_back({
          arr3f{
            float(CGAL::to_double(p1.x())),
            float(CGAL::to_double(p1.y())),
            float(CGAL::to_double(p1.z()))},
          arr3f{
            float(CGAL::to_double(p2.x())),
            float(CGAL::to_double(p2.y())),
            float(CGAL::to_double(p2.z()))},
        });
      }

      return segments_out.size();
    } else {
      for (const auto& e : sorted_segments) {
        segments_out.push_back(
          LD.project(e.first,e.second)
        );
      }
      return sorted_segments.size();
    }
  } else return 0;
  
}

void DetectLinesNode::process(){
  auto input_geom = vector_input("edge_points");
  auto planes = input("pts_per_roofplane").get<IndexedPlanesWithPoints>();
  

  SegmentCollection edge_segments, lines3d;
  vec1i ring_order, ring_id, is_start;
  std::unordered_map<size_t,std::vector<size_t>> ring_idx;
  // fit lines in all input points
  // if (input_geom.is_connected_type(typeid(PointCollection))) {
  //   std::vector<linedect::Point> cgal_pts;
  //   auto points = input_geom.get<PointCollection>();
  //   for( auto& p : points ) {
  //     cgal_pts.push_back(linedect::Point(p[0], p[1], p[2]));
  //   }
  //   linedect::LineDetector LD(cgal_pts);
  //   detect_lines(LD);
  //   LD.get_bounded_edges(edge_segments);

  // fit lines per ring
  // } else if (input_geom.is_connected_type(typeid(LinearRingCollection))) {
  auto roofplane_ids = input("roofplane_ids").get<vec1i>();
  int n = k;
  
  size_t seg_cntr=0, plane_id;
  for (size_t i=0; i<input_geom.size(); ++i) {
    auto& ring = input_geom.get<LinearRing>(i);
    plane_id = roofplane_ids[i];
    std::vector<linedect::Point> cgal_pts;
    for( auto& p : ring ) {
      cgal_pts.push_back(linedect::Point(p[0], p[1], p[2]));
    }

    linedect::LineDetector LD(cgal_pts);
    // SegmentCollection ring_edges;
    auto n_detected = detect_lines_ring_m2(LD, planes[plane_id].first, edge_segments);
    LD.get_bounded_edges(lines3d);

    for (size_t j=0; j<n_detected; ++j) {
      // edge_segments.push_back(ring_edges[j]);
      ring_idx[plane_id].push_back(seg_cntr++);
      ring_order.push_back(j);
      ring_id.push_back(plane_id);
      ring_order.push_back(j);
      ring_id.push_back(plane_id);
      is_start.push_back(1);
      is_start.push_back(0);
    }

    // also check the holes
    for (auto& hole : ring.interior_rings()) {
      std::vector<linedect::Point> cgal_pts;
      for( auto& p : hole ) {
        cgal_pts.push_back(linedect::Point(p[0], p[1], p[2]));
      }

      linedect::LineDetector LD(cgal_pts);
      // SegmentCollection ring_edges;
      auto n_detected = detect_lines_ring_m2(LD, planes[plane_id].first, edge_segments);
      LD.get_bounded_edges(lines3d);

      for (size_t j=0; j<n_detected; ++j) {
        // edge_segments.push_back(ring_edges[j]);
        ring_idx[plane_id].push_back(seg_cntr++);
        ring_order.push_back(j);
        ring_id.push_back(plane_id);
        ring_order.push_back(j);
        ring_id.push_back(plane_id);
        is_start.push_back(1);
        is_start.push_back(0);
      }
    }
    // std::cout << "number of shapes: " << LD.segment_shapes.size() <<"\n";
    // std::cout << "number of segments: " << order_cnt <<"\n";
  }
  // }

  output("edge_segments").set(edge_segments);
  output("lines3d").set(lines3d);
  // output("ring_edges").set(ring_edges);
  output("ring_idx").set(ring_idx);
  output("ring_id").set(ring_id);
  output("ring_order").set(ring_order);
  output("is_start").set(is_start);
}

// void ClassifyEdgePointsNode::process(){
//   // Set up vertex data (and buffer(s)) and attribute pointers
//   auto points = input("points").get<PNL_vector>();

//   config c;
//   c.classify_jump_count_min = param<int>("classify_jump_count_min");
//   c.classify_jump_count_max = param<int>("classify_jump_count_max");
//   c.classify_line_dist = param<float>("classify_line_dist");
//   c.classify_jump_ele = param<float>("classify_jump_ele");

//   std::vector<linedect::Point> edge_points;
//   classify_edgepoints(edge_points, points, c);
//   output("edge_points").set(edge_points);

//   vec3f edge_points_vec3f;
//   for(auto& p : edge_points) {
//       std::array<float,3> a = {{
//         float(p.x()), 
//         float(p.y()), 
//         float(p.z())
//       }};
//       edge_points_vec3f.push_back(a);
//     }
//   output("edge_points_vec3f").set(edge_points_vec3f);
// }

float compute_percentile(std::vector<float>& z_vec, float percentile) {
  assert(percentile<=1.);
  assert(percentile>=0.);
  size_t n = (z_vec.size()-1) * percentile;
  std::nth_element(z_vec.begin(), z_vec.begin()+n, z_vec.end());
  return z_vec[n];
}

void BuildingSelectorNode::process() {
  polygon_count = vector_input("point_clouds").size();
  if ((building_id < polygon_count) && 
  (vector_input("point_clouds").size() == vector_input("polygons").size())) {
    auto& point_cloud = vector_input("point_clouds").get<PointCollection&>(building_id);
    auto& ground_point_cloud = vector_input("ground_point_clouds").get<PointCollection&>(building_id);
    auto& polygon = vector_input("polygons").get<LinearRing&>(building_id);
    auto& ground_elevation = vector_input("ground_elevations").get<float&>(building_id);

    output("point_cloud").set(point_cloud);
    output("ground_point_cloud").set(ground_point_cloud);
    output("polygon").set(polygon);
    output("ground_elevation").set(ground_elevation);
  }
};

void RegulariseLinesNode::process(){
  // Set up vertex data (and buffer(s)) and attribute pointers
  auto edges = input("edge_segments").get<SegmentCollection>();
  auto ints_segments = input("ints_segments").get<SegmentCollection>();
  // auto footprint = input("footprint").get<LinearRing>();
  // auto ring_id = input("ring_id").get<vec1i>();
  // auto ring_order = input("ring_order").get<vec1i>();

  // linereg::Polygon_2 cgal_footprint;
  // std::vector<linereg::Polygon_2> ek_holes;
  // for (auto& p : footprint) {
  //   cgal_footprint.push_back(linereg::Point_2(p[0], p[1]));
  // }
  // for (auto& ring : footprint.interior_rings()) {
  //   linereg::Polygon_2 ek_hole;
  //   for (auto& p : ring) {
  //     ek_hole.push_back(linereg::Point_2(p[0], p[1]));
  //   }
  //   ek_holes.push_back(ek_hole);
  // }

  // get clusters from line regularisation 
  auto LR = linereg::LineRegulariser();
  LR.add_segments(0,edges);
  // LR.add_segments(1,cgal_footprint, (double) 0);
  // for (auto& hole : ek_holes) {
  //   LR.add_segments(1,hole, (double) 0);
  // }
  LR.add_segments(2,ints_segments);
  LR.dist_threshold = dist_threshold*dist_threshold;
  LR.angle_threshold = angle_threshold;

  std::cout << "\nangle clustering...";
  LR.perform_angle_clustering();
  std::cout << "distance clustering...";
  LR.perform_distance_clustering();
  std::cout << "...clustering complete\n";

  // output("exact_footprint_out").set(linereg::Polygon_with_holes_2(cgal_footprint, ek_holes.begin(), ek_holes.end()));

  auto& regularised = vector_output("regularised_edges");
  regularised.touch();
  auto& regularised_exact = vector_output("exact_regularised_edges");
  regularised_exact.touch();
  SegmentCollection edges_out_;
  vec1i priorities, angle_cluster_ids, dist_cluster_ids;
  // we should iterate of the distance clusters and output one segment per cluster
  for(auto& line : LR.lines) {
    linereg::Segment_2 segment;
    segment = line.segment;
    auto new_seg = Segment();
    new_seg[0] = {float(CGAL::to_double(segment.source().x())), float(CGAL::to_double(segment.source().y())), 0};
    new_seg[1] = {float(CGAL::to_double(segment.target().x())), float(CGAL::to_double(segment.target().y())), 0};
    edges_out_.push_back(new_seg);
    priorities.push_back(line.priority);
    priorities.push_back(line.priority);
    angle_cluster_ids.push_back(line.angle_cluster_id);
    angle_cluster_ids.push_back(line.angle_cluster_id);
    dist_cluster_ids.push_back(line.dist_cluster_id);
    dist_cluster_ids.push_back(line.dist_cluster_id);
  }
  for(auto& dclust : LR.dist_clusters) {
    size_t dclust_size = dclust->lines.size();
    Segment_2 segment;
    //get lines with highest priority
    size_t max_priority=0;
    for (auto line : dclust->lines) {
      if(line->priority > max_priority) max_priority = line->priority;
    }
    std::vector<linereg::linetype*> prio_lines;
    for (auto line : dclust->lines) {
      if(line->priority == max_priority) {
        prio_lines.push_back(line);
      }
    }
    // TODO: split clusters using common planes if intersection segment(s) are present? 
    // TODO: compute distance clusters with 3D dists? 
    // TODO: output exact segments? quick solve of spikes?
    // TODO: performance optimise clustering algo
    // //skip if cluster only contains fp segments
    // if(max_priority==1 && (dclust_size == prio_lines.size())) continue;
    if(!merge_intersection_lines && (max_priority==2 && (prio_lines.size()>1))) 
    {
      std::vector<linereg::linetype*> other_lines;
      for (auto line : dclust->lines) {
        if(line->priority != 2) {
          other_lines.push_back(line);
        }
      }
      for (auto line : prio_lines) {
        double mean_angle = line->angle;
        auto centroid = line->midpoint;
        segment = calc_segment(centroid, mean_angle, other_lines, extension);
        auto new_seg = Segment();
        new_seg[0] = {float(CGAL::to_double(segment.source().x())), float(CGAL::to_double(segment.source().y())), 0};
        new_seg[1] = {float(CGAL::to_double(segment.target().x())), float(CGAL::to_double(segment.target().y())), 0};
        regularised.push_back(new_seg);
        regularised_exact.push_back(segment);
      }
    } 
    else 
    {
      //compute mean line with small extensions on both ends
      double mean_angle = calc_mean_angle(prio_lines);
      auto centroid = calc_centroid(prio_lines);
      segment = calc_segment(centroid, mean_angle, dclust->lines, extension);
      auto new_seg = Segment();
      new_seg[0] = {float(CGAL::to_double(segment.source().x())), float(CGAL::to_double(segment.source().y())), 0};
      new_seg[1] = {float(CGAL::to_double(segment.target().x())), float(CGAL::to_double(segment.target().y())), 0};
      regularised.push_back(new_seg);
      regularised_exact.push_back(segment);
    }
  }
  output("angle_cluster_id").set(angle_cluster_ids);
  output("dist_cluster_id").set(dist_cluster_ids);
  output("priorities").set(priorities);
  output("edges_out_").set(edges_out_);
  // output("edges_out").set(new_segments);
  // output("rings_out").set(new_rings);
  // output("footprint_out").set(new_fp);
}

void chain(Segment& a, Segment& b, LinearRing& ring, const float& snap_threshold) {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  K::Line_2 l_a(K::Point_2(a[0][0], a[0][1]), K::Point_2(a[1][0], a[1][1]));
  K::Line_2 l_b(K::Point_2(b[0][0], b[0][1]), K::Point_2(b[1][0], b[1][1]));
  K::Segment_2 s(K::Point_2(a[1][0], a[1][1]), K::Point_2(b[0][0], b[0][1]));
  auto result = CGAL::intersection(l_a, l_b);
  if (result) {
    if (auto p = boost::get<K::Point_2>(&*result)) {
      if (CGAL::squared_distance(*p, s) < snap_threshold) {
        ring.push_back({float(p->x()), float(p->y()), 0});
      } else {
        ring.push_back(a[1]);
        ring.push_back(b[0]);
      }
    }
  // } else if (auto l = boost::get<K::Line_2>(&*result)) {

  // }
  } else { // there is no intersection
    ring.push_back(a[1]);
    ring.push_back(b[0]);
  }
}
void RegulariseRingsNode::process(){
  // // Set up vertex data (and buffer(s)) and attribute pointers
  // auto edges = input("edge_segments").get<SegmentCollection>();
  // auto ints_segments = input("ints_segments").get<SegmentCollection>();
  // auto ring_idx = input("ring_idx").get<std::unordered_map<size_t,std::vector<size_t>>>();
  // auto footprint = input("footprint").get<LinearRing>();
  // // auto ring_id = input("ring_id").get<vec1i>();
  // // auto ring_order = input("ring_order").get<vec1i>();

  // // SegmentCollection all_edges;

  // // build vector of all input edges
  // // for(auto edge : edges) {
  // //   all_edges.push_back(edge);
  // // }
  // // size_t fpi_begin = all_edges.size();
  // // SegmentCollection fp_edges;
  // // for(size_t i=0; i<footprint.size()-1; ++i) {
  // //   fp_edges.push_back(Segment({footprint[i], footprint[i+1]}));
  // // }
  // // fp_edges.push_back(
  // //   Segment({footprint[footprint.size()-1], footprint[0]})
  // // );
  // linereg::Polygon_2 ek_footprint;
  // std::vector<linereg::Polygon_2> ek_holes;
  // for (auto& p : footprint) {
  //   ek_footprint.push_back(linereg::Point_2(p[0], p[1]));
  // }
  // for (auto& ring : footprint.interior_rings()) {
  //   linereg::Polygon_2 ek_hole;
  //   for (auto& p : ring) {
  //     ek_hole.push_back(linereg::Point_2(p[0], p[1]));
  //   }
  //   ek_holes.push_back(ek_hole);
  // }
  // // size_t fpi_end = all_edges.size()-1;

  // // get clusters from line regularisation 
  // auto LR = linereg::LineRegulariser();
  // LR.add_segments(0,edges);
  // LR.add_segments(1,ek_footprint, (double) fp_offset);
  // for (auto& hole : ek_holes) {
  //   LR.add_segments(1,hole, (double) fp_offset);
  // }
  // LR.add_segments(2,ints_segments);
  // LR.dist_threshold = dist_threshold;
  // LR.angle_threshold = angle_threshold;
  // LR.cluster(weighted_avg, angle_per_distcluster);

  // if (regularise_fp) {
  //   std::vector<size_t> fp_idx;
  //   for (size_t i=0; i < LR.get_segments(1).size(); ++i) {
  //     fp_idx.push_back(i);
  //   }
  //   auto ek_reg_fp = linereg::chain_ring<linereg::EK>(fp_idx, LR.get_segments(1), snap_threshold);
  //   output("exact_footprint_out").set(linereg::Polygon_with_holes_2(ek_reg_fp));
  // } else {
  //   output("exact_footprint_out").set(linereg::Polygon_with_holes_2(ek_footprint, ek_holes.begin(), ek_holes.end()));
  // }

  // if (recompute_rings) {
  //   std::unordered_map<size_t, linereg::Polygon_2> exact_polygons;
  //   for (auto& idx : ring_idx) {
  //     exact_polygons[idx.first] = 
  //       linereg::chain_ring<linereg::EK>(idx.second, LR.get_segments(0), snap_threshold);
  //     // std::cout << "ch ring size : "<< exact_polygons.back().size() << ", " << exact_polygons.back().is_simple() << "\n";
  //   }
  //   // std::cout << "ch fp size : "<< exact_fp.size() << ", " << exact_fp.is_simple() << "\n";
  //   output("exact_rings_out").set(exact_polygons);


  //   LinearRingCollection lrc;
  //   vec1i plane_ids;
  //   for (auto& poly : exact_polygons) {
  //     LinearRing lr;
  //     for (auto p=poly.second.vertices_begin(); p!=poly.second.vertices_end(); ++p) {
  //       lr.push_back({
  //         float(CGAL::to_double(p->x())),
  //         float(CGAL::to_double(p->y())),
  //         0
  //       });
  //       plane_ids.push_back(poly.first);
  //     }
  //     lrc.push_back(lr);
  //   }
  //   output("rings_out").set(lrc);
  //   output("plane_id").set(plane_ids);

  // }
  // auto& new_segments = vector_output("edges_out");
  // SegmentCollection edges_out_;
  // vec1i priorities, angle_cluster_ids, dist_cluster_ids;
  // for(auto& line : LR.lines) {
  //   auto& ekseg = line.reg_segment;
  //   auto new_seg = Segment();
  //   new_seg[0] = {float(CGAL::to_double(ekseg.source().x())), float(CGAL::to_double(ekseg.source().y())), 0};
  //   new_seg[1] = {float(CGAL::to_double(ekseg.target().x())), float(CGAL::to_double(ekseg.target().y())), 0};
  //   new_segments.push_back(new_seg);
  //   edges_out_.push_back(new_seg);
  //   priorities.push_back(line.priority);
  //   priorities.push_back(line.priority);
  //   angle_cluster_ids.push_back(line.angle_cluster_id);
  //   angle_cluster_ids.push_back(line.angle_cluster_id);
  //   dist_cluster_ids.push_back(line.dist_cluster_id);
  //   dist_cluster_ids.push_back(line.dist_cluster_id);
  // }

  // output("angle_cluster_id").set(angle_cluster_ids);
  // output("dist_cluster_id").set(dist_cluster_ids);
  // output("priorities").set(priorities);
  // output("edges_out_").set(edges_out_);
  // // output("edges_out").set(new_segments);
  // // output("rings_out").set(new_rings);
  // // output("footprint_out").set(new_fp);
}

LinearRing simplify_footprint(const LinearRing& polygon, float& threshold_stop_cost) {
  namespace PS = CGAL::Polyline_simplification_2;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef K::Point_2 Point_2;
  typedef CGAL::Polygon_2<K>                   Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K>        Polygon_with_holes_2;
  typedef PS::Stop_below_count_ratio_threshold Stop_count_ratio;
  typedef PS::Stop_above_cost_threshold        Stop_cost;
  typedef PS::Squared_distance_cost            Cost;

  if (polygon.size()>2) {
      Polygon_2 cgal_exterior;
      Cost cost;

      for (auto& p : polygon) {
        cgal_exterior.push_back(Point_2(p[0], p[1]));
      }

      std::vector<Polygon_2> cgal_holes;
      for (auto& iring : polygon.interior_rings()) {
        Polygon_2 hole;
        for (auto& p : iring) {
          hole.push_back(Point_2(p[0], p[1]));
        }
        cgal_holes.push_back(hole);
      }
      // cgal_polygon.erase(cgal_polygon.vertices_end()-1); // remove repeated point from the boost polygon
      
      // polygon = PS::simplify(polygon, cost, Stop_count_ratio(0.5));

      auto cgal_polygon = PS::simplify(Polygon_with_holes_2(cgal_exterior, cgal_holes.begin(), cgal_holes.end()), cost, Stop_cost(threshold_stop_cost));
      
      LinearRing footprint_vec3f;
      for (auto v = cgal_polygon.outer_boundary().vertices_begin(); v!=cgal_polygon.outer_boundary().vertices_end(); v++){
        footprint_vec3f.push_back({float(v->x()),float(v->y()),0});
      }
      for(auto hole = cgal_polygon.holes_begin(); hole!=cgal_polygon.holes_end(); ++hole) {
        vec3f gf_hole;
        for (auto v = hole->vertices_begin(); v!=hole->vertices_end(); v++){
          gf_hole.push_back({float(v->x()),float(v->y()),0});
        }
        footprint_vec3f.interior_rings().push_back(gf_hole);
      }

      // HACK: CGAL does not seem to remove the first point of the input polygon in any case, so we need to check ourselves
      auto p_0 = *(cgal_polygon.outer_boundary().vertices_begin());
      auto p_1 = *(cgal_polygon.outer_boundary().vertices_begin()+1);
      auto p_end = *(cgal_polygon.outer_boundary().vertices_end()-1);
      // check the distance between the first vertex and the line between its 2 neighbours
      if (CGAL::squared_distance(Point_2(p_0), K::Segment_2(p_end, p_1)) < threshold_stop_cost) {
        footprint_vec3f.erase(footprint_vec3f.begin());
      }

      if (footprint_vec3f.size()>2) {
        return footprint_vec3f;
      } else {
        return polygon;
      }
    } else
      return polygon;
}

bool get_line_extend(
   Kernel::Line_3* l,
   std::vector<Point>& points,
   double& dmin,
   double& dmax,
   Point& pmin,
   Point& pmax,
   float& min_dist_to_line_sq
) {
  size_t cnt = 0;
  bool setminmax=false;  
  double sqd_min = 1 + min_dist_to_line_sq;

  auto lp = l->point();
  auto lv = l->to_vector();
  lv = lv / CGAL::sqrt(lv.squared_length());

  for (auto p : points) {
    auto sqd = CGAL::squared_distance(*l, p);
    if (sqd > min_dist_to_line_sq)
      continue;
    auto av = p-lp;
    auto d = av*lv;
    if (!setminmax) {
      setminmax=true;
      dmin=dmax=d;
      pmin=pmax=p;
    }
    if (d < dmin){
      dmin = d;
      pmin = p;
    }
    if (d > dmax) {
      dmax = d;
      pmax = p;
    }
    sqd_min = std::min(sqd_min, sqd);
    ++cnt;
  }
  return cnt > 1;
}

void PlaneIntersectNode::process() {
  auto pts_per_roofplane = input("pts_per_roofplane").get<IndexedPlanesWithPoints>();
  auto plane_adj = input("plane_adj").get<std::map<size_t, std::map<size_t, size_t>>>();
  // // auto& alpha_rings = vector_input("alpha_rings");//.get<LinearRingCollection>();

  float min_dist_to_line_sq = min_dist_to_line*min_dist_to_line;

  SegmentCollection lines;
  size_t ring_cntr=0;
  for(auto& [id_hi, ids_lo] : plane_adj) {
    auto& plane_hi = pts_per_roofplane[id_hi].first;
    auto& plane_pts_hi = pts_per_roofplane[id_hi].second;
    // // auto& alpha_ring = alpha_rings.get<LinearRing>(ring_cntr++);
    for (auto& [id_lo, cnt] : ids_lo) {
      // skip plain pairs with  low number of neighbouring points
      if (cnt < min_neighb_pts) continue;
      auto& plane_lo = pts_per_roofplane[id_lo].first;
      auto& plane_pts_lo = pts_per_roofplane[id_lo].second;
      auto result = CGAL::intersection(plane_hi, plane_lo);
      if (result) {
        // bound the line to extend of one plane's inliers
        if (auto l = boost::get<typename Kernel::Line_3>(&*result)) {
          Point pmin_lo, pmax_lo;
          Point pmin_hi, pmax_hi;
          double dmin_lo, dmax_lo;
          double dmin_hi, dmax_hi;
          
          // skip this line if it is too far away any of the plane_pts
          if( !get_line_extend(l, plane_pts_hi, dmin_hi, dmax_hi, pmin_hi, pmax_hi, min_dist_to_line_sq) || 
              !get_line_extend(l, plane_pts_lo, dmin_lo, dmax_lo, pmin_lo, pmax_lo, min_dist_to_line_sq) )
            continue;

          // take the overlap between the two extends
          Point ppmin, ppmax;
          if (dmin_lo > dmin_hi)
            ppmin = l->projection(pmin_lo);
          else
            ppmin = l->projection(pmin_hi);
          
          if (dmax_lo < dmax_hi)
            ppmax = l->projection(pmax_lo);
          else
            ppmax = l->projection(pmax_hi);
          
          // Check for infinity (quick fix for linux crash)
          auto sx = float(CGAL::to_double(ppmin.x()));
          auto sy = float(CGAL::to_double(ppmin.y()));
          auto tx = float(CGAL::to_double(ppmax.x()));
          auto ty = float(CGAL::to_double(ppmax.y()));
          if (!((std::isinf(sx) || std::isinf(sy)) || (std::isinf(tx) || std::isinf(ty)))) {
            if(CGAL::squared_distance(ppmin, ppmax) > 1E-10) {
              arr3f source = {
                sx,
                sy,
                float(CGAL::to_double(ppmin.z()))
              };
              arr3f target = {
                tx,
                ty,
                float(CGAL::to_double(ppmax.z()))
              };
              lines.push_back({source,target});
            }
          }
        }
      }
    }
  }

  output("lines").set(lines);
}

void SimplifyPolygonNode::process(){
  // Set up vertex data (and buffer(s)) and attribute pointers

  auto& geom_in = vector_input("polygons");
  auto& geom_out = vector_output("polygons_simp");

  LinearRingCollection polygons_out;
  for (size_t i=0; i< geom_in.size(); ++i) {
    auto& polygon = geom_in.get<LinearRing&>(i);
    geom_out.push_back(
      simplify_footprint(polygon, threshold_stop_cost)
    );
  }
  
}

void ArrDissolveNode::process() {
  auto arr = input("arrangement").get<Arrangement_2>();
  auto& heightfield = input("heightfield").get<RasterTools::Raster>();
  
  if(dissolve_seg_edges) {
    Face_merge_observer obs(arr);
    arr_dissolve_seg_edges(arr);
  }
  
  if (dissolve_all_interior) {
    arr_dissolve_fp(arr, true, false);
  }

  if (dissolve_step_edges) {
    Face_merge_observer obs(arr);
    for (auto face: arr.face_handles()) {
      if (face->data().in_footprint) {
        LinearRing polygon;
        arrangementface_to_polygon(face, polygon);

        auto height_points = heightfield.rasterise_polygon(polygon, true);
        size_t datasize = height_points.size(), data_cnt = 0;
        if(datasize==0) { //polygon was too small to yield any pixels/height_points
          auto& pz = polygon[0][2];
          face->data().elevation_50p = pz;
          face->data().elevation_70p = pz;
          face->data().elevation_min = pz;
          face->data().elevation_max = pz;
          face->data().data_coverage = 1;
        } else {
          auto& plane = face->data().plane;
          // compute elevations based on the assigned plane on this face
          for (auto& p : height_points) {
            if(p[2]!=heightfield.noDataVal_) data_cnt++;
            p[2] = -plane.a()/plane.c() * p[0] - plane.b()/plane.c()*p[1] - plane.d()/plane.c();
          }
          std::sort(height_points.begin(), height_points.end(), [](auto& p1, auto& p2) {
            return p1[2] < p2[2];
          });
          int elevation_id = std::floor(0.5*float(datasize-1));
          face->data().elevation_50p = height_points[elevation_id][2];
          elevation_id = std::floor(0.7*float(datasize-1));
          face->data().elevation_70p = height_points[elevation_id][2];
          face->data().elevation_min = height_points[0][2];
          face->data().elevation_max = height_points[datasize-1][2];
          face->data().data_coverage = float(data_cnt) / float(datasize);
        }
      }
    }
    arr_dissolve_step_edges(arr, step_height_threshold);
  }

  // remove any dangling edges
  {
    std::vector<Arrangement_2::Halfedge_handle> to_remove;
    for (auto he : arr.edge_handles()) {
      if (he->face()==he->twin()->face())
        to_remove.push_back(he);
    }
    for (auto he : to_remove) {
      arr.remove_edge(he);
    }
  }

  // remove all lines not inside footprint
  if (dissolve_outside_fp) {
    arr_dissolve_fp(arr, false, true);
  }

  // label different building parts
  arr_label_buildingparts(arr);
  
  // ensure all holes are marked as such
  auto f_unb = arr.unbounded_face();
  for (auto fh : arr.face_handles()) {
    if(fh!=f_unb) 
      if(!fh->data().in_footprint && !fh->data().is_footprint_hole) {
        fh->data().is_footprint_hole = true;
      }
  }

  // if (remove_duplicates) {
  //   arr_snap_duplicates(arr, (double) std::pow(10,-dupe_threshold_exp));
  // }
  // snap edge shorter than snap_tolerance
  // for (auto he : arr.edge_handles()) {
  //   if(CGAL::squared_distance(he->source()->point(), he->target()->point()) < snap_tolerance) {

  //     arr.remove_edge(he->next());
  //     arr.remove_edge(he->twin()->previous());
  //   }
  // }
  // compute data_area and elevation stats for each face, only for lod22 (lod13 needs to have this before dissolving stepedges)
  if (!dissolve_step_edges) {
    for (auto face: arr.face_handles()) {
      if (face->data().in_footprint) {
        LinearRing polygon;
        arrangementface_to_polygon(face, polygon);

        auto height_points = heightfield.rasterise_polygon(polygon, true);
        size_t data_cnt = 0;
        bool skip_face = height_points.size()==0;
        {
          // auto& plane = face->data().plane;
          // compute elevations based on the assigned plane on this face
          std::vector<arr3f*> pts;
          if (!skip_face){
            for (auto& p : height_points) {
              if(p[2]!=heightfield.noDataVal_) {
                data_cnt++;
                pts.push_back(&p);
              };
            }
          }
          if(data_cnt==0 || skip_face) {
            auto& pz = polygon[0][2];
            face->data().elevation_50p = pz;
            face->data().elevation_70p = pz;
            face->data().elevation_min = pz;
            face->data().elevation_max = pz;
            face->data().data_coverage = 1;
          } else {
            std::sort(pts.begin(), pts.end(), [](auto& p1, auto& p2) {
              return (*p1)[2] < (*p2)[2];
            });
            size_t datasize = pts.size();
            int elevation_id = std::floor(0.5*float(datasize-1));
            face->data().elevation_50p = (*(pts[elevation_id])) [2];
            elevation_id = std::floor(0.7*float(datasize-1));
            face->data().elevation_70p = (*(pts[elevation_id])) [2];
            face->data().elevation_min = (*(pts[0])) [2];
            face->data().elevation_max = (*(pts[datasize-1])) [2];
            face->data().data_coverage = float(data_cnt) / float(height_points.size());
          }
        }
      }
    }
  }

  output("arrangement").set(arr);
}


//this one is not working properly yet
void PolygonUnionNode::process() {
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_2 Point;
  typedef CGAL::Polygon_2<Kernel> Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;

  auto& polygons = vector_input("polygons");

  // use create_interior_skeleton_and_offset_polygons_2 to buffer input?

  std::vector<Polygon_2> cgal_polygons;
  for (size_t i=0; i<polygons.size(); ++i) {
    auto& polygon = polygons.get<LinearRing&>(i);
    Polygon_2 cgal_poly;
    for (auto& p : polygon) {
      cgal_poly.push_back(Point(p[0], p[1]));
    }
    cgal_polygons.push_back(cgal_poly);
  }
  std::vector<Polygon_with_holes_2> res;
  CGAL::join(cgal_polygons.begin(), cgal_polygons.end(), std::back_inserter (res));

  auto& polygons_outer = vector_output("polygons");
  auto& holes = vector_output("holes");
  // polygons_outer.resize<PointCollection>(polygons.size());
  for (auto it = res.begin(); it!=res.end(); ++it) {
    LinearRing l;
    auto& outer_bound = it->outer_boundary();
    for(auto vit = outer_bound.vertices_begin(); vit!=outer_bound.vertices_end(); ++vit) {
      l.push_back({float(vit->x()), float(vit->y()), 0});
    }
    polygons_outer.push_back(l);
    for(auto hole = it->holes_begin(); hole!=it->holes_end(); ++hole) {
      LinearRing h;
      for(auto vit = hole->vertices_begin(); vit!=hole->vertices_end(); ++vit) {
        l.push_back({float(vit->x()), float(vit->y()), 0});
      }
      holes.push_back(h);
    }
  }
}

void PointCloudMergerNode::process() {
  
  for (size_t i=0; i<poly_input("pointclouds").sub_terminals()[0]->size(); ++i)
  {
    PointCollection point_cloud;
    for (auto &iterm : poly_input("pointclouds").sub_terminals()) {
      if (iterm->has_data()) {
        auto& pc = iterm->get<PointCollection>(i);
        point_cloud.insert(point_cloud.begin(), pc.begin(), pc.end());
      }
    }
    output("pointcloud").push_back(point_cloud);
  }
}

void SegmentExtendNode::process() {
  typedef Kernel::Segment_3 Segment;
  typedef Kernel::Point_3 Point;

  auto segments = input("segments").get<SegmentCollection>();

  SegmentCollection segments_out;
  for (auto& segment : segments)  {
    Segment s(
      Point(segment[0][0], segment[0][1], segment[0][2]),
      Point(segment[1][0], segment[1][1], segment[1][2])
    );
    auto va = (s.target()-s.source());
    va = va/CGAL::sqrt(va.squared_length());
    auto se = Segment( 
      s.source() - va*extension,
      s.target() + va*extension  
    );
    segments_out.push_back (
      {
        arr3f{
          float(se.source().x()),
          float(se.source().y()),
          float(se.source().z())
        },
        arr3f{
          float(se.target().x()),
          float(se.target().y()),
          float(se.target().z())
        }
      }
    );
  }
  output("segments").set(segments_out);
}

void PCFilterNode::process() {
  auto& pointcloud = input("pointcloud").get<PointCollection&>();
  

  PointCollection pointcloud_filtered;

  size_t i=0;
  for(auto& p : pointcloud) {
    if(input("values").get<float>(i++) > threshold) {
      pointcloud_filtered.push_back(p);
    }
  }

  output("pointcloud").set(pointcloud_filtered);
}

// void PlaneDetectorNode::process() {
//   auto points = input("point_clouds").get<Feature>();

//   planedect::PlaneDetector PD(points_vec, normals_vec);
//   PD.dist_thres = c.metrics_plane_epsilon * c.metrics_plane_epsilon;
//   PD.normal_thres = c.metrics_plane_normal_threshold;
//   PD.min_segment_count = c.metrics_plane_min_points;
//   PD.N = c.metrics_normal_k;
//   PD.detect();
//   std::cout << PD.segment_shapes.size() << " shapes detected." << std::endl;

//   // // Instantiates shape detection engine.
//   // Region_growing shape_detection;

//   // // Sets parameters for shape detection.
//   // Region_growing::Parameters parameters;
//   // // Sets probability to miss the largest primitive at each iteration.
//   // // parameters.probability = 0.05;
//   // // Detect shapes with at least 500 points.
//   // parameters.min_points = c.metrics_plane_min_points;
//   // // Sets maximum Euclidean distance between a point and a shape.
//   // parameters.epsilon = c.metrics_plane_epsilon;
//   // // Sets maximum Euclidean distance between points to be clustered.
//   // // parameters.cluster_epsilon = 0.01;
//   // // Sets maximum normal deviation.
//   // // 0.9 < dot(surface_normal, point_normal); 
//   // parameters.normal_threshold = c.metrics_plane_normal_threshold;

//   // // Provides the input data.
//   // shape_detection.set_input(points);
//   // // Registers planar shapes via template method.
//   // shape_detection.add_shape_factory<SCPlane>();
//   // // Detects registered shapes with parameters.
//   // std::cout << "points.size: " << points.size() << "\n";
//   // shape_detection.detect(parameters);
//   // // Prints number of detected shapes.
//   // std::cout << shape_detection.shapes().end() - shape_detection.shapes().begin() << " shapes detected." << std::endl;

//   i=1;
//   for(auto seg: PD.segment_shapes){
//     auto& plane = seg.second;
//     auto plane_idx = PD.get_point_indices(seg.first);
//   }

//   output("decomposed_footprints").set(all_cells);
//   output("attributes").set(all_attributes);
// }

}
