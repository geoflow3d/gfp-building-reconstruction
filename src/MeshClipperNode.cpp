#include "MeshProcessingNodes.hpp"

#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace geoflow::nodes::stepedge {

  void MeshClipperNode::process() { 
    
    typedef SurfaceMesh::Vertex_index       VertexIndex;
    typedef SurfaceMesh::Face_index         FaceIndex;

    auto& gfmesh = input("mesh").get<Mesh>();
    
    SurfaceMesh smesh;
    {
      std::map<arr3f, VertexIndex> vertex_map;
      std::set<arr3f> vertex_set;
      for (const auto &ring : gfmesh.get_polygons())
      {
        for (auto &v : ring)
        {
          auto [it, did_insert] = vertex_set.insert(v);
          if (did_insert)
          {
            vertex_map[v] = smesh.add_vertex(K::Point_3(v[0],v[1],v[2]));;
          }
        }
      }
    
      for (auto& ring : gfmesh.get_polygons()) {
        std::vector<VertexIndex> rindices;
        rindices.reserve(ring.size());
        for(auto& p : ring) {
          rindices.push_back(vertex_map[p]);
        }
        smesh.add_face(rindices);
      }
    }

    if(!skip_clip_){
      auto& gfbox = input("bbox").get<Box>();
      const auto& pmin = gfbox.min();
      const auto& pmax = gfbox.max();
      K::Iso_cuboid_3 cuboid(
        K::Point_3(pmin[0], pmin[1], pmin[2]), 
        K::Point_3(pmax[0], pmax[1], pmax[2])
      );

      // clip
      if(!CGAL::is_triangle_mesh(smesh)) CGAL::Polygon_mesh_processing::triangulate_faces(smesh);
      
      if(!CGAL::Polygon_mesh_processing::does_self_intersect(smesh)) {
        CGAL::Polygon_mesh_processing::clip(
          smesh,
          cuboid
        );
      }
    }
    output("cgal_surface_mesh").set(smesh);
  }

}