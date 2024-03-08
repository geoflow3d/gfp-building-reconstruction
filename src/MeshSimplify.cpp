#include "MeshProcessingNodes.hpp"

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Midpoint placement policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
//Placement wrapper
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/boost/graph/iterator.h>
#include <cmath>
#include <geoflow/common.hpp>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>

#include "tinsimp.hpp"

namespace geoflow::nodes::stepedge {

  typedef SurfaceMesh::Vertex_index       VertexIndex;
  typedef SurfaceMesh::Face_index         FaceIndex;

  typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor    vertex_descriptor;
  typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef boost::graph_traits<SurfaceMesh>::edge_descriptor      edge_descriptor;
  namespace SMS = CGAL::Surface_mesh_simplification;
  // BGL property map which indicates whether an edge is marked as non-removable
  struct Border_is_constrained_edge_map
  {
    const SurfaceMesh* sm_ptr;
    typedef edge_descriptor                                       key_type;
    typedef bool                                                  value_type;
    typedef value_type                                            reference;
    typedef boost::readable_property_map_tag                      category;
    Border_is_constrained_edge_map(const SurfaceMesh& sm) : sm_ptr(&sm) {}
    friend value_type get(const Border_is_constrained_edge_map& m, const key_type& edge) {
      return CGAL::is_border(edge, *m.sm_ptr);
    }
  };
  // Placement class
  typedef SMS::Constrained_placement<SMS::Midpoint_placement<SurfaceMesh>,
                                    Border_is_constrained_edge_map > Placement;
  // namespace d = CGAL::Polygon_mesh_processing;

  struct MeshBuilder {
    std::map<K::Point_3, std::size_t> vertex_map;
    std::set<K::Point_3> vertex_set;
    std::vector<K::Point_3> points;

    std::vector<std::vector<std::size_t> > polygons;

    void add_point(const K::Point_3& p) {
      auto [it, did_insert] = vertex_set.insert(p);
      if (did_insert)
      {
        vertex_map[p] = points.size();
        points.push_back(p);
      }
    }
    void add_points(const K::Point_3& p0, const K::Point_3& p1, const K::Point_3& p2) {
      add_point(p0);
      add_point(p1);
      add_point(p2);
    }

    void add_triangle(const K::Point_3& p0, const K::Point_3& p1, const K::Point_3& p2) {
      // First build a polygon soup
      add_points(p0, p1, p2);
      std::vector<std::size_t> rindices;
      rindices.reserve(3);
      rindices.push_back(vertex_map[p0]);
      rindices.push_back(vertex_map[p1]);
      rindices.push_back(vertex_map[p2]);
      polygons.push_back(rindices);
    }

    void get_mesh(geoflow::nodes::stepedge::SurfaceMesh& smesh) {
      smesh.clear();
      CGAL::Polygon_mesh_processing::repair_polygon_soup(points, polygons);

      // duplicate non-manifold edges (but does not re-orient faces)
      CGAL::Polygon_mesh_processing::duplicate_non_manifold_edges_in_polygon_soup(points, polygons);
      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, smesh);
    }

  };

  const K::Vector_3 up(0,0,1);

  bool is_vertical(const K::Point_3& a, const K::Point_3& b, const K::Point_3& c) {
    auto n = CGAL::orthogonal_vector(a, b, c);
    // check if face is vertical by checking angle of it's normal. Only add constraints for non-vertical faces
    // std::cout << CGAL::abs(CGAL::approximate_angle(n, up))<< std::endl;
    return (CGAL::abs(90 - CGAL::abs(CGAL::approximate_angle(n, up))) == 0);
  }

  void MeshSimplify2DNode::process() {
    auto smesh = input("cgal_surface_mesh").get<SurfaceMesh>();

    // !d::does_self_intersect(smesh)
    if (error_ > 0){
      tinsimp::CDT t;

      // collect vertical faces that make up the vertical wall
      std::vector<std::vector<K::Point_3>> wall_triangles;
      for (auto face : smesh.faces()) {
        std::vector<K::Point_3> triangle;
        for(vertex_descriptor vd : vertices_around_face(smesh.halfedge(face), smesh)){
          triangle.push_back(smesh.point(vd));
        }
        if (is_vertical(triangle[0], triangle[1], triangle[2])) {
          wall_triangles.push_back(triangle);
          CGAL::Euler::remove_face(smesh.halfedge(face), smesh);
        }
      }
      smesh.collect_garbage();

      for(halfedge_descriptor hd : halfedges(smesh))
      {
        if(CGAL::is_border(hd, smesh))
        {
          auto a = smesh.point(source(hd, smesh));
          auto b = smesh.point(target(hd, smesh));
            
          t.insert_constraint(
            tinsimp::Point(
              a.x(),
              a.y(),
              a.z()
            ),
            tinsimp::Point(
              b.x(),
              b.y(),
              b.z()
            )
          );
        }
      }

      std::vector<tinsimp::Point> zpts;
      for(vertex_descriptor vd : smesh.vertices())
      {
        if(!CGAL::is_border(vd, smesh))
        {
          auto p = smesh.point(vd);
          zpts.push_back(
            tinsimp::Point( p.x(), p.y(), p.z() )
          );
        }
      }

      tinsimp::mark_domains(t);
      float sq_area = 0;
      for (auto& fh : t.finite_face_handles()) {
        if(fh->info().in_domain()) {
          sq_area += t.triangle(fh).squared_area();
        }
      }
      float total_area = std::sqrt(sq_area);
      
      tinsimp::greedy_insert(t, zpts, error_, minpts_*total_area);

      // reset and recompute nesting levels
      for (auto& fh : t.all_face_handles()) {
        fh->info().nesting_level = -1;
      }
      tinsimp::mark_domains(t);
      // std::cout << "\nFinished!\n" << r << " edges removed.\n"
                // << smesh.number_of_edges() << " final edges.\n";

      smesh.clear();
      MeshBuilder mb;
      for (auto& fh : t.finite_face_handles()) {
        if(fh->info().in_domain()) {
            mb.add_triangle(
              K::Point_3(
                fh->vertex(0)->point().x(),
                fh->vertex(0)->point().y(),
                fh->vertex(0)->point().z()
              ),
              K::Point_3(
                fh->vertex(1)->point().x(),
                fh->vertex(1)->point().y(),
                fh->vertex(1)->point().z()
              ),
              K::Point_3(
                fh->vertex(2)->point().x(),
                fh->vertex(2)->point().y(),
                fh->vertex(2)->point().z()
              )
            );
        }
      }
      for (auto& triangle : wall_triangles) {
        mb.add_triangle(triangle[0], triangle[1], triangle[2]);
      }
      mb.get_mesh(smesh);

      // TriangleCollection tc;
      // for ( auto& t : wall_triangles ) {
      //   Triangle gft;
      //   gft[0] = arr3f{static_cast<float>(t[0].x()), static_cast<float>(t[0].y()), static_cast<float>(t[0].z())};
      //   gft[1] = arr3f{static_cast<float>(t[1].x()), static_cast<float>(t[1].y()), static_cast<float>(t[1].z())};
      //   gft[2] = arr3f{static_cast<float>(t[2].x()), static_cast<float>(t[2].y()), static_cast<float>(t[2].z())};
      //   tc.push_back(gft);
      // }
      // output("wall_triangles").set(tc);
    }
    output("cgal_surface_mesh").set(smesh);
  }

  void MeshSimplifyNode::process() {
    auto smesh = input("cgal_surface_mesh").get<SurfaceMesh>();

    // !d::does_self_intersect(smesh)
    if (stop_ratio_ < 1.){
      if(!CGAL::is_triangle_mesh(smesh)) CGAL::Polygon_mesh_processing::triangulate_faces(smesh);
      
      SurfaceMesh::Property_map<halfedge_descriptor, std::pair<K::Point_3, K::Point_3> > constrained_halfedges;
      constrained_halfedges = smesh.add_property_map<halfedge_descriptor,std::pair<K::Point_3, K::Point_3> >("h:vertices").first;
      size_t n_border=0;
      for(halfedge_descriptor hd : halfedges(smesh))
      {
        if(CGAL::is_border(hd, smesh))
        {
          constrained_halfedges[hd] = std::make_pair(smesh.point(source(hd, smesh)),
                                                    smesh.point(target(hd, smesh)));
          ++n_border;
        }
      }
      float stop_ratio = stop_ratio_;
      if(border_correction_) {
        size_t n_all = smesh.number_of_halfedges();
        stop_ratio = ((n_all-n_border) * stop_ratio_) / n_all;
      }
      // Contract the surface mesh as much as possible. Correct for the border edges that will not be removed
      SMS::Count_ratio_stop_predicate<SurfaceMesh> stop( stop_ratio );
      Border_is_constrained_edge_map bem(smesh);
      
      // filter that checks if a placement would invert the normal of a face 
      SMS::Bounded_normal_change_filter<> filter;
      // This the actual call to the simplification algorithm.
      // The surface mesh and stop conditions are mandatory arguments.

      int r = SMS::edge_collapse(smesh, stop,
                                CGAL::parameters::edge_is_constrained_map(bem)
                                                  .filter(filter)
                                                  .get_placement(Placement(bem)));
      // std::cout << "\nFinished!\n" << r << " edges removed.\n"
                // << smesh.number_of_edges() << " final edges.\n";
    }
    output("cgal_surface_mesh").set(smesh);
  }

  void SurfaceMesh2OFFNode::process() {
    
    auto smesh = input("cgal_surface_mesh").get<SurfaceMesh>();
    auto fname = manager.substitute_globals(filepath_);

    std::ofstream ofs;
    ofs << std::fixed << std::setprecision(5);
    ofs.open(fname);
    CGAL::IO::write_OFF(ofs, smesh);
    ofs.close();
  }

}