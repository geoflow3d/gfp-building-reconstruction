#include "stepedge_nodes.hpp"

namespace as {

  class AlphaShapeRegionGrower {
    Alpha_shape_2 &A;
    enum Mode {
      ALPHA, // stop at alpha boundary
      EXTERIOR // stop at faces labels as exterior
    };
    int label_cnt; // label==-1 means exterior, -2 mean never visiter, 0+ means a regular region
    public:
    std::unordered_map<Face_handle, int> face_map;
    std::unordered_map<int, Vertex_handle> region_map; //label: (boundary vertex)
    AlphaShapeRegionGrower(Alpha_shape_2& as) : A(as), label_cnt(0) {};

    void grow() {
      std::stack<Face_handle> seeds;
      for (auto fh = A.all_faces_begin(); fh!=A.all_faces_end(); ++fh) {
        seeds.push(fh);
        face_map[fh] = -2;
      }
      auto inf_face = A.infinite_face();
      face_map[inf_face] = -1;
      grow_region(inf_face, ALPHA); // this sets label of exterior faces to -1
      while (!seeds.empty()) {
        auto fh = seeds.top(); seeds.pop();
        if (face_map[fh] == -2) {
          face_map[fh] = label_cnt;
          grow_region(fh, EXTERIOR);
          ++label_cnt;
        }
      }
    }

    void grow_region (Face_handle face_handle, Mode mode) {
      std::stack<Face_handle> candidates;
      candidates.push(face_handle);

      while(candidates.size()>0) {
        auto fh = candidates.top(); candidates.pop();
        // check the 3 neighbors of this face
        for (int i=0; i<3; ++i) {
          auto e = std::make_pair(fh,i);
          auto neighbor = fh->neighbor(i);

          if (mode == ALPHA) {
            // add neighbor if it is not on the ohter side of alpha boundary
            // check if this neighbor hasn't been visited before
            if (face_map[neighbor] == -2) {
              auto edge_class = A.classify(e);
              if ( ! (edge_class == Alpha_shape_2::REGULAR || edge_class == Alpha_shape_2::SINGULAR) ) {
                face_map[neighbor] = -1;
                candidates.push(neighbor);
              }
            }
          } else if (mode == EXTERIOR) {
            // check if this neighbor hasn't been visited before and is not exterior
            if (face_map[neighbor] == -2) {
              face_map[neighbor] = label_cnt;
              candidates.push(neighbor);
            // if it is exterior, we store this boundary edge
            } else if (face_map[neighbor] == -1) {
              if( region_map.find(label_cnt)==region_map.end() ) {
                region_map[label_cnt] = fh->vertex(A.cw(i));
              }
            }
          }

        }
      }
    }
  };
}

namespace geoflow::nodes::stepedge {

void AlphaShapeNode::process(){
  auto points_per_segment = input("pts_per_roofplane").get<IndexedPlanesWithPoints>();

  PointCollection edge_points, boundary_points;
  LineStringCollection alpha_edges;
  LinearRingCollection alpha_rings;
  TriangleCollection alpha_triangles;
  std::vector<as::Triangulation_2> alpha_dts;
  vec1i segment_ids, plane_idx;
  for (auto& it : points_per_segment ) {
    if (it.first == -1) continue; // skip points if they put at index -1 (eg if we care not about slanted surfaces for ring extraction)
    auto points = it.second.second;
    as::Triangulation_2 T;
    T.insert(points.begin(), points.end());
    as::Alpha_shape_2 A(T,
                as::FT(thres_alpha),
                as::Alpha_shape_2::GENERAL);
    
    double optimal_alpha = thres_alpha;
    if (optimal_alpha && optimal_only_if_needed) {
      optimal_alpha = std::max(float(*A.find_optimal_alpha(1)), thres_alpha);
    } else if (optimal_alpha) {
      optimal_alpha = *A.find_optimal_alpha(1);
    }
    A.set_alpha(as::FT(optimal_alpha));

    for (auto it = A.alpha_shape_vertices_begin(); it!=A.alpha_shape_vertices_end(); it++) {
      auto p = (*it)->point();
      edge_points.push_back({float(p.x()), float(p.y()), float(p.z())});
    }
    for (auto it = A.alpha_shape_edges_begin(); it!=A.alpha_shape_edges_end(); it++) {
      auto p1 = it->first->vertex(A.cw(it->second))->point();
      auto p2 = it->first->vertex(A.ccw(it->second))->point();

      alpha_edges.push_back({
        {float(p1.x()), float(p1.y()), float(p1.z())},
        {float(p2.x()), float(p2.y()), float(p2.z())}
      });
    }
    
    // flood filling 
    auto grower = as::AlphaShapeRegionGrower(A);
    grower.grow();

    for (auto fh = A.finite_faces_begin(); fh != A.finite_faces_end(); ++fh) {
        arr3f p0 = {float (fh->vertex(0)->point().x()), float (fh->vertex(0)->point().y()), float (fh->vertex(0)->point().z())};
        arr3f p1 = {float (fh->vertex(1)->point().x()), float (fh->vertex(1)->point().y()), float (fh->vertex(1)->point().z())};
        arr3f p2 = {float (fh->vertex(2)->point().x()), float (fh->vertex(2)->point().y()), float (fh->vertex(2)->point().z())
        };
      alpha_triangles.push_back({ p0,p1,p2 });
      segment_ids.push_back(grower.face_map[fh]);
      segment_ids.push_back(grower.face_map[fh]);
      segment_ids.push_back(grower.face_map[fh]);
    }

    for (auto& kv : grower.region_map) {
      auto region_label = kv.first;
      auto v_start = kv.second;
      boundary_points.push_back({
        float(v_start->point().x()),
        float(v_start->point().y()),
        float(v_start->point().z())
      });

      // find edges of outer boundary in order
      LinearRing ring;

      ring.push_back( {float(v_start->point().x()), float(v_start->point().y()), float(v_start->point().z())} );
      // secondly, walk along the entire boundary starting from v_start
      as::Vertex_handle v_next, v_prev = v_start, v_cur = v_start;
      size_t v_cntr = 0;
      do {
        as::Edge_circulator ec(A.incident_edges(v_cur)), done(ec);
        do {
          // find the vertex on the other side of the incident edge ec
          auto v = ec->first->vertex(A.cw(ec->second));
          if (v_cur == v) v = ec->first->vertex(A.ccw(ec->second));
          // find labels of two adjacent faces
          auto label1 = grower.face_map[ ec->first ];
          auto label2 = grower.face_map[ ec->first->neighbor(ec->second) ];
          // check if the edge is on the boundary of the region and if we are not going backwards
          bool exterior = label1==-1 || label2==-1;
          bool region = label1==region_label || label2==region_label;
          if(( exterior && region )  && (v != v_prev)) {
            v_next = v;
            ring.push_back( {float(v_next->point().x()), float(v_next->point().y()), float(v_next->point().z())} );
            break;
          }
        } while (++ec != done);
        v_prev = v_cur;
        v_cur = v_next;

      } while (v_next != v_start);
      // finally, store the ring 
      alpha_rings.push_back(ring);
      plane_idx.push_back(it.first);
      alpha_dts.push_back(T);
    }
  }
  
  output("alpha_rings").set(alpha_rings);
  output("alpha_edges").set(alpha_edges);
  output("alpha_triangles").set(alpha_triangles);
  output("alpha_dts").set(alpha_dts);
  output("segment_ids").set(segment_ids);
  output("edge_points").set(edge_points);
  output("boundary_points").set(boundary_points);
  output("roofplane_ids").set(plane_idx);
}

}