#include "MeshSimplifyNode.hpp"
#include "Simplify.h"

namespace geoflow::nodes::stepedge {

  void MeshSimplifyNode::process() { 
    auto& triangles = input("triangles").get<TriangleCollection>();

    Simplify::vertices.clear();
    Simplify::triangles.clear();

    std::map<arr3f, size_t> vertex_map;
    // std::vector<arr3f> vertex_vec;
    {
      size_t v_cntr = 0;
      std::set<arr3f> vertex_set;
      for (auto &triangle : triangles)
      {
        for (auto &vertex : triangle)
        {
          auto [it, did_insert] = vertex_set.insert(vertex);
          if (did_insert)
          {
            Simplify::Vertex vo;
            vo.p.x = vertex[0];
            vo.p.y = vertex[1];
            vo.p.z = vertex[2];
            vertex_map[vertex] = v_cntr++;
            // vertex_vec.push_back(vertex);
            Simplify::vertices.push_back(vo);
          }
        }
      }
    }
    
    for (auto& triangle : triangles) {
      Simplify::Triangle t;
      t.v[0] = vertex_map[triangle[0]];
      t.v[1] = vertex_map[triangle[1]];
      t.v[2] = vertex_map[triangle[2]];
      t.attr = 0;
      Simplify::triangles.push_back(t);
    }

    int target_count = round((float)Simplify::triangles.size() * reduce_fraction_);
    
    Simplify::simplify_mesh(target_count, agressiveness_, true);

    TriangleCollection triangles_simplified;

    for (auto& triangle : Simplify::triangles) {
      Triangle t;
      t[0] = arr3f{ 
        (float) Simplify::vertices[triangle.v[0]].p.x,
        (float) Simplify::vertices[triangle.v[0]].p.y,
        (float) Simplify::vertices[triangle.v[0]].p.z
      };
      t[1] = arr3f{ 
        (float) Simplify::vertices[triangle.v[1]].p.x,
        (float) Simplify::vertices[triangle.v[1]].p.y,
        (float) Simplify::vertices[triangle.v[1]].p.z
      };
      t[2] = arr3f{ 
        (float) Simplify::vertices[triangle.v[2]].p.x,
        (float) Simplify::vertices[triangle.v[2]].p.y,
        (float) Simplify::vertices[triangle.v[2]].p.z
      };
      triangles_simplified.push_back(t);
    }

    output("triangles").set(triangles_simplified);

  }

}