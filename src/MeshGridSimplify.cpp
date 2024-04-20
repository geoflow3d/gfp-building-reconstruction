#include "MeshProcessingNodes.hpp"
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace geoflow::nodes::stepedge {

  struct Grid3D {
    Box box;
    // bool x = true;
    unsigned dim[3];
    float cell_size_xy, cell_size_z;
    Grid3D(Box& box, float& cell_size_xy, float& cell_size_z) : box(box), cell_size_xy(cell_size_xy), cell_size_z(cell_size_z) {
      dim[0] = (box.max()[0]-box.min()[0])/cell_size_xy + 1;
      dim[1] = (box.max()[1]-box.min()[1])/cell_size_xy + 1;
      dim[2] = (box.max()[2]-box.min()[2])/cell_size_z + 1;
      // std::cout << "dim[0]: " << dim[0] << "[" << box.min()[0] << ", " << box.max()[0] << "]" << std::endl;
      // std::cout << "dim[1]: " << dim[1] << "[" << box.min()[1] << ", " << box.max()[1] << "]" << std::endl;
      // std::cout << "dim[2]: " << dim[2] << "[" << box.min()[2] << ", " << box.max()[2] << "]" << std::endl;
    }

    unsigned getCellCoordinate(const K::Point_3& p) {
      unsigned c = static_cast<unsigned>( floor((p.x()-box.min()[0]) / cell_size_xy) );
      unsigned r = static_cast<unsigned>( floor((p.y()-box.min()[1]) / cell_size_xy) );
      unsigned s = static_cast<unsigned>( floor((p.z()-box.min()[2]) / cell_size_z) );
      
      return s * dim[0]*dim[1] + r * dim[0] + c;
    }

    arr3f getCellCenterPoint(const unsigned& cell_coordinate) {
      unsigned slice = cell_coordinate / (dim[0]*dim[1]);
      unsigned rest = cell_coordinate % (dim[0]*dim[1]);
      unsigned row = rest / dim[0];
      unsigned col = rest % dim[0];
      arr3f p;
      p[0] = box.min()[0] + col*cell_size_xy + cell_size_xy/2;
      p[1] = box.min()[1] + row*cell_size_xy + cell_size_xy/2;
      p[2] = box.min()[2] + slice*cell_size_z + cell_size_z/2;
      // if (x) {
      //   std::cout << "cell_coordinate: " << cell_coordinate << std::endl;
      //   std::cout << "slice: " << slice << std::endl;
      //   std::cout << "row: " << row << std::endl;
      //   std::cout << "col: " << col << std::endl;
      //   x=false;
      // }
      return p;
    }
  };
  struct Grid2D {
    struct SumCount {
      float sum;
      unsigned cnt;
      SumCount() : sum(0), cnt(0) {};
      SumCount(float val) : sum(val), cnt(1) {};
      void add(float val) { sum+=val; ++cnt; };
      float get_avg() { return sum/cnt; };
    };
    Box box;
    // bool x = true;
    unsigned dim[2];
    
    std::unordered_map<unsigned, SumCount> elevations;
    float cell_size_xy, cell_size_z;
    Grid2D(Box& box, float& cell_size_xy, float& cell_size_z) : box(box), cell_size_xy(cell_size_xy), cell_size_z(cell_size_z) {
      auto dx = fmod(box.min()[0], cell_size_xy);
      auto dy = fmod(box.min()[1], cell_size_xy);
      box.add(arr3f{box.min()[0]-dx, box.min()[1]-dy, box.min()[2]});
      dim[0] = (box.max()[0]-box.min()[0])/cell_size_xy + 1;
      dim[1] = (box.max()[1]-box.min()[1])/cell_size_xy + 1;
      // dim[2] = (box.max()[2]-box.min()[2])/cell_size_z + 1;
      // std::cout << "dim[0]: " << dim[0] << "[" << box.min()[0] << ", " << box.max()[0] << "]" << std::endl;
      // std::cout << "dim[1]: " << dim[1] << "[" << box.min()[1] << ", " << box.max()[1] << "]" << std::endl;
      // std::cout << "dim[2]: " << dim[2] << "[" << box.min()[2] << ", " << box.max()[2] << "]" << std::endl;
    }

    unsigned getCellCoordinate(const K::Point_3& p) {
      unsigned c = static_cast<unsigned>( floor((p.x()-box.min()[0]) / cell_size_xy) );
      unsigned r = static_cast<unsigned>( floor((p.y()-box.min()[1]) / cell_size_xy) );
      unsigned coord = r * dim[0] + c;
      if(elevations.count(coord)) {
        elevations[coord].add(p.z());
      } else {
        elevations[coord] = SumCount(p.z());
      }
      return coord;
    }

    arr3f getCellCenterPoint(const unsigned& cell_coordinate) {
      unsigned row = cell_coordinate / dim[0];
      unsigned col = cell_coordinate % dim[0];
      arr3f p;
      p[0] = box.min()[0] + col*cell_size_xy + cell_size_xy/2;
      p[1] = box.min()[1] + row*cell_size_xy + cell_size_xy/2;
      p[2] = elevations[cell_coordinate].get_avg();
      // if (x) {
      //   std::cout << "cell_coordinate: " << cell_coordinate << std::endl;
      //   std::cout << "slice: " << slice << std::endl;
      //   std::cout << "row: " << row << std::endl;
      //   std::cout << "col: " << col << std::endl;
      //   x=false;
      // }
      return p;
    }
  };

  void MeshGridSimplifyNode::process() {
    typedef SurfaceMesh::Vertex_index       VertexIndex;
    namespace PMP = CGAL::Polygon_mesh_processing;
    SurfaceMesh smesh;
    
    // if(input("mesh").is_connected_type(typeid(Mesh)))
    auto gfmesh = input("mesh").get<Mesh>();
    
    Box box;
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
            box.add(v);
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

    if(!CGAL::is_triangle_mesh(smesh)) PMP::triangulate_faces(smesh);

    // separate vertices and faces

    // assign vertices to 3D grid cell
    // merge vertices in each grid cell, and create a mapping (old vertex indices -> new vertex indices, std::map<unsigned, unsigned>)
    Grid2D G(box, cell_size_xy_, cell_size_z_);
    std::unordered_map<unsigned, unsigned> vertex_map;
    for (auto vi : smesh.vertices()) {
      unsigned c = G.getCellCoordinate( smesh.point(vi) );
      vertex_map[vi] = c;
    }

    // map for each face from old to new vertices and check how many are left, if <3 remove the face]

    std::unordered_map<unsigned, VertexIndex> new_mesh_vertices;
    SurfaceMesh smesh_new;
    {
      // struct Cluster {
      //   float x,y;
      //   std::vector<float> zs;
      //   Cluster() : x(0), y(0) {};
      //   Cluster(float x, float y) : x(x), y(y) {};
      // };
      std::unordered_map<unsigned, arr3f> new_vertices;
      for (auto f : smesh.faces()) {
        for(VertexIndex vi : vertices_around_face(smesh.halfedge(f), smesh)) {
          auto p = G.getCellCenterPoint(vertex_map[vi]);
          new_mesh_vertices[vertex_map[vi]] = smesh_new.add_vertex(K::Point_3(p[0], p[1], p[2]));
        }
      }
    }
    {
      for (auto f : smesh.faces()) {
        std::set<unsigned> face_set;
        for(VertexIndex vi : vertices_around_face(smesh.halfedge(f), smesh)) {
          face_set.insert(vertex_map[vi]);
        }
        if (face_set.size() == 3) {
          std::vector<VertexIndex> rindices;
          rindices.reserve(3);          
          for(VertexIndex vi : vertices_around_face(smesh.halfedge(f), smesh)) {
            rindices.push_back(new_mesh_vertices[vertex_map[vi]]);
          }
          smesh_new.add_face(rindices);
        }
      }
    }
    // map for each face from old to new vertices and check how many are left, if <3 remove the face
    // TriangleCollection triangleCollection;
    // vec3f normals;
    // for (auto& f : smesh.faces()) {      
    //   std::set<unsigned> face_set;
    //   for(VertexIndex vi : vertices_around_face(smesh.halfedge(f), smesh)) {
    //     face_set.insert(vertex_map[vi]);
    //   }
    //   if (face_set.size() == 3) {
    //     Triangle t;
    //     unsigned i=0;
    //     for(VertexIndex vi : vertices_around_face(smesh.halfedge(f), smesh)) {
    //       // auto& p = smesh.point(vi);
    //       auto p = G.getCellCenterPoint(vertex_map[vi]);
    //       t[i++] = p;
    //     }
    //     triangleCollection.push_back(t);
    //   }
      
    // }

    // output("tri").set(triangleCollection);
    output("mesh").set(smesh_new);
  }

}