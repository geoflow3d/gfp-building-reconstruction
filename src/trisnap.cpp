#include "stepedge_nodes.hpp"
#include "polygon_util.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Arr_walk_along_line_point_location.h>


namespace geoflow::nodes::stepedge {

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Epeck;
  typedef CGAL::Exact_predicates_tag Tag;

  typedef CGAL::Triangulation_vertex_base_2<K> VertexBase;
  typedef CGAL::Constrained_triangulation_face_base_2<K> FaceBase;
  typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo*, K, FaceBase> FaceBaseWithInfo;
  typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TriangulationDataStructure, Tag> CDT;
  typedef CDT::Edge_circulator Edge_circulator;
  typedef CDT::Face_circulator Face_circulator;
  typedef CDT::Vertex_handle Vertex_handle;
  typedef CDT::Face_handle Face_handle;
  typedef std::pair<Face_handle, int> Edge;


  void TriSnapNode::process() {

    typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;

    auto arr = input("arrangement").get<Arrangement_2>();
    
    CDT cdt;
    float sq_dist_thres = dist_thres*dist_thres;

    // map from arr vertices to cdt vertices
    std::unordered_map<Arrangement_2::Vertex_handle, CDT::Vertex_handle> vertex_map;

    // Segment_list_2 seg_list;
    // Polyline_list_2 output_list;
    for (auto arrVertex : arr.vertex_handles()) {
      // auto& p = v->point();
      vertex_map[arrVertex] = cdt.insert(CDT::Point_2(
        CGAL::to_double(arrVertex->point().x()), 
        CGAL::to_double(arrVertex->point().y())
      ));
    }


    Walk_pl walk_pl (arr);
    for (auto& arrEdge : arr.edge_handles()) {
      cdt.insert_constraint(
        vertex_map[ arrEdge->source() ] ,
        vertex_map[ arrEdge->target() ]
      );
    }

    for (auto fit : cdt.finite_face_handles()) {

      auto p = CGAL::centroid(cdt.triangle(fit));
      auto obj = walk_pl.locate( Walk_pl::Arrangement_2::Point_2(p.x(), p.y()) );

      // std::cout << "The point (" << p << ") is located ";
      if (auto f = boost::get<Face_const_handle>(&obj)) { // located inside a face
        // std::cout << "inside "
        //           << (((*f)->is_unbounded()) ? "the unbounded" : "a bounded")
        //           << " face." << std::endl;
        fit->info() = &(arr.non_const_handle(*f))->data(); 
      }
      // else if (auto e = boost::get<Halfedge_const_handle>(&obj)) // located on an edge
      //   std::cout << "on an edge: " << (*e)->curve() << std::endl;
      // else if (auto v = boost::get<Vertex_const_handle>(&obj)) // located on a vertex
      //   std::cout << "on " << (((*v)->is_isolated()) ? "an isolated" : "a")
      //             << " vertex: " << (*v)->point() << std::endl;
      // else CGAL_error_msg("Invalid object.");
    }

    TriangleCollection triangles_og;
    vec1i segment_ids_og;
    for (auto fh = cdt.finite_faces_begin(); fh != cdt.finite_faces_end(); ++fh) {
      // only export triangles in the interior of a shape (thus excluding holes and exterior)

        arr3f p0 = {float (fh->vertex(0)->point().x()), float (fh->vertex(0)->point().y()), 0};
        arr3f p1 = {float (fh->vertex(1)->point().x()), float (fh->vertex(1)->point().y()), 0};
        arr3f p2 = {float (fh->vertex(2)->point().x()), float (fh->vertex(2)->point().y()), 0};
        triangles_og.push_back({ p0,p1,p2 });
        segment_ids_og.push_back(fh->info()->segid);
        segment_ids_og.push_back(fh->info()->segid);
        segment_ids_og.push_back(fh->info()->segid);
    }
    output("triangles_og").set(triangles_og);
    output("segment_ids_og").set(segment_ids_og);

    // Detect triangles with 3 short edges => collapse triangle to point (remove 2 vertices)
    bool checked_all = false;
    while(!checked_all) {
      checked_all = false;
      for (auto fit : cdt.finite_face_handles()) {
        auto v0 = fit->vertex(0);
        auto v1 = fit->vertex(1);
        auto v2 = fit->vertex(2);
        auto& p0 = v0->point();
        auto& p1 = v1->point();
        auto& p2 = v2->point();
        // auto e0 = std::make_pair(fit, 0);
        // auto e1 = std::make_pair(fit, 1);
        // auto e2 = std::make_pair(fit, 2);
        if (
          (
            CGAL::squared_distance(p0, p1) > sq_dist_thres &&
            CGAL::squared_distance(p1, p2) > sq_dist_thres &&
            CGAL::squared_distance(p2, p0) > sq_dist_thres
          ) //&& (
          //   cdt.is_constrained( e0 ) &&
          //   cdt.is_constrained( e1 ) &&
          //   cdt.is_constrained( e2 )
          // )
        ) {
          std::unordered_map<Vertex_handle, FaceInfo*> constraints_to_restore;
          auto pnew = CGAL::centroid(p0,p1,p2);

          // collect incident constraint edges
          // std::vector<Edge> to_remove;
          for (size_t i=0; i<3; ++i) {
            auto vthis = fit->vertex(i);
            Edge_circulator ec = cdt.incident_edges(vthis),
            done(ec);
            if (ec != nullptr) {
              do {
                if ( cdt.is_constrained( *ec ) ) {
                  // figure out which is the other vertex on this edge
                  // and what is the face counter clockwise to the edge vthis -> vother
                  auto va = ec->first->vertex(cdt.cw(ec->second));
                  auto vb = ec->first->vertex(cdt.ccw(ec->second));
                  Vertex_handle vother;
                  FaceInfo* face_ccw;
                  if (va == vthis) {
                    vother = vb;
                    face_ccw = ec->first->neighbor(ec->second)->info();
                  } else {
                    vother = va;
                    face_ccw = ec->first->info();
                  };

                  // check if edge is not an edge of fit
                  int cw = cdt.cw(i);
                  // std::cout << cw << std::endl;
                  int ccw = cdt.ccw(i);
                  // std::cout << ccw << std::endl;
                  auto vcw = fit->vertex(cw);
                  auto vccw = fit->vertex(ccw);
                  if ( vother != fit->vertex(cdt.cw(i)) && vother != fit->vertex(cdt.ccw(i)) ) {
                    constraints_to_restore[vother] = face_ccw;
                  }
                  // to_remove.push_back(*ec);
                }
              } while ( ++ec != done );
            }
          }
          // for (auto& e : to_remove) {
          //   cdt.remove_constrained_edge(e.first, e.second);
          // }

          // collapse the triangle to centroid
          cdt.remove_incident_constraints(v0);
          cdt.remove(v0);
          cdt.remove_incident_constraints(v1);
          cdt.remove(v1);
          cdt.remove_incident_constraints(v2);
          cdt.remove(v2);
          auto vnew = cdt.insert(pnew);

          // restore constraints
          for (auto& [vh, finfo] : constraints_to_restore) {
            cdt.insert_constraint(vnew, vh);

            // find the ccw face to the edge vnew->vh
            Face_circulator fc = cdt.incident_faces(vnew),
              done(fc);
            if (fc != nullptr) {
              do { 
                if (fc->has_vertex(vh)) {
                  // check if fc is ccw to edge vnew->vh
                  if (vh == fc->vertex( cdt.ccw( fc->index(vnew) ) ) ) {
                    fc->info() = finfo;
                  } else {
                    // fc is cw to edge, so find the neighbor on the other side of edge
                    fc->neighbor( cdt.cw( fc->index(vh) ) )->info() = finfo;
                  }
                }
              } while (++fc != done);
            }
          }


          checked_all = false;
          break;
        } else { checked_all = true; }
      }

    }
    // Detect triangles with 1 short edge  => collapse the short edge to point (remove one vertex)
    // Detect triangles with 1 vertex close to opposing edge  => split opposing edge to contain the vertex

    // void 	remove_incident_constraints (Vertex_handle v)
    // void 	remove_constrained_edge (Face_handle f, int i)
    // void 	remove (Vertex_handle v)


    // // Arrangement_2 arr_overlay;
    // // Overlay_traits overlay_traits;
    // // overlay(arr, arrSR, arr_overlay, overlay_traits);

    // output("arrangement").set(arrSR);

    TriangleCollection triangles_snapped;
    vec1i segment_ids_snapped;
    for (auto fh = cdt.finite_faces_begin(); fh != cdt.finite_faces_end(); ++fh) {
      // only export triangles in the interior of a shape (thus excluding holes and exterior)

        arr3f p0 = {float (fh->vertex(0)->point().x()), float (fh->vertex(0)->point().y()), 0};
        arr3f p1 = {float (fh->vertex(1)->point().x()), float (fh->vertex(1)->point().y()), 0};
        arr3f p2 = {float (fh->vertex(2)->point().x()), float (fh->vertex(2)->point().y()), 0};
        triangles_snapped.push_back({ p0,p1,p2 });
        // segment_ids_snapped.push_back(fh->info()->segid);
        // segment_ids_snapped.push_back(fh->info()->segid);
        // segment_ids_snapped.push_back(fh->info()->segid);
    }
    output("triangles_snapped").set(triangles_snapped);
    output("segment_ids_snapped").set(segment_ids_snapped);

  }
}