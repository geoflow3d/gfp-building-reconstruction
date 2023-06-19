#include "MeshProcessingNodes.hpp"

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
// Midpoint placement policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
//Placement wrapper
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
// Stop-condition policy
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <CGAL/Surface_mesh/IO/OFF.h>

namespace geoflow::nodes::stepedge {

  typedef SurfaceMesh::Vertex_index       VertexIndex;
  typedef SurfaceMesh::Face_index         FaceIndex;

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
  namespace PMP = CGAL::Polygon_mesh_processing;

  void MeshSimplifyNode::process() {
    auto smesh = input("cgal_surface_mesh").get<SurfaceMesh>();

    // !PMP::does_self_intersect(smesh)
    if (stop_ratio_ < 1.){
      if(!CGAL::is_triangle_mesh(smesh)) CGAL::Polygon_mesh_processing::triangulate_faces(smesh);

      // this prevents potentially getting stuck in infinite loop (see https://github.com/CGAL/cgal/issues/7529)
      CGAL::Polygon_mesh_processing::duplicate_non_manifold_vertices( smesh );
      CGAL::Polygon_mesh_processing::remove_isolated_vertices( smesh	);
      
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
      // This the actual call to the simplification algorithm.
      // The surface mesh and stop conditions are mandatory arguments.

      int r = SMS::edge_collapse(smesh, stop,
                                CGAL::parameters::edge_is_constrained_map(bem)
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