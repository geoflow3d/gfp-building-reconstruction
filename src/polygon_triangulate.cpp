// This file is part of gfp-building-reconstruction
// Copyright (C) 2018-2022 Ravi Peters

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.

// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include <CGAL/Kernel/global_functions_3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2.h>

#include <stepedge_nodes.hpp>
#include "point_edge.h"
#include "polygon_util.hpp"
#include "cdt_util.hpp"

namespace geoflow::nodes::stepedge {

namespace tri = tri_util;

typedef CGAL::Polygon_2<tri::K> Polygon_2;
typedef CGAL::Plane_3<tri::K> Plane_3;

glm::vec3 calculate_normal(const LinearRing ring)
{
  glm::vec3 normal(0, 0, 0);
  for (size_t i = 0; i < ring.size(); ++i) {
    const auto &curr = ring[i];
    const auto &next = ring[(i + 1) % ring.size()];
    normal[0] += (curr[1] - next[1]) * (curr[2] + next[2]);
    normal[1] += (curr[2] - next[2]) * (curr[0] + next[0]);
    normal[2] += (curr[0] - next[0]) * (curr[1] + next[1]);
  }
  return glm::normalize(normal);
}

double calculate_volume(const TriangleCollection& triangle_collection) {
  double sum = 0;
  for(const auto& t : triangle_collection) {
    auto a = Vector(t[0][0], t[0][1], t[0][2]);
    auto b = Vector(t[1][0], t[1][1], t[1][2]);
    auto c = Vector(t[2][0], t[2][1], t[2][2]);
    sum += CGAL::scalar_product(a, CGAL::cross_product(b, c));
  }
  return sum/6;
}

// void mark_domains(CDT& cdt) {

//   std::list<CDT::Face_handle> explorables;
//   auto fh_inf = cdt.infinite_face();
//   fh_inf->info().set_interior(false);
//   explorables.push_back(fh_inf);

//   while (!explorables.empty()) {
//     auto fh = explorables.front();
//     explorables.pop_front();
//     // check the three neighbours
//     for (int i = 0; i < 3; ++i) {
//       CDT::Edge e(fh, i);
//       CDT::Face_handle fh_n = fh->neighbor(i);
//       if (fh_n->info().visited == false) {
//         if (cdt.is_constrained(e)) {
//           fh_n->info().set_interior(!fh_n->info().interior);
//         } else {
//           fh_n->info().set_interior(fh_n->info().interior);
//         }
//         explorables.push_back(fh_n);
//       }
//     }
//   }
// }

Polygon_2 project(geoflow::vec3f& ring, Plane_3& plane) {
  Polygon_2 poly_2d;
  for (auto& p : ring) {
    poly_2d.push_back(plane.to_2d(tri::K::Point_3(p[0], p[1], p[2])));
  }
  return poly_2d;
}

void project_and_insert(geoflow::vec3f& ring, Plane_3& plane, tri::CDT& cdt) {
  auto pit_last = ring.end()-1;
  tri::CDT::Vertex_handle vh_next, vh_last, vh = cdt.insert(plane.to_2d(tri::K::Point_3((*pit_last)[0], (*pit_last)[1], (*pit_last)[2])));
  vh_last = vh;
  vh->info().set_point(*pit_last);
  for (auto pit = ring.begin(); pit != ring.end(); ++pit) {
    if(pit==pit_last){
      vh_next=vh_last;
    } else {
      vh_next = cdt.insert(plane.to_2d(tri::K::Point_3((*pit)[0], (*pit)[1], (*pit)[2])));
      vh_next->info().set_point(*pit);
    }
    cdt.insert_constraint(vh, vh_next);
    vh = vh_next;
  }
}

void PolygonTriangulatorNode::triangulate_polygon(LinearRing& poly, vec3f& normals, TriangleCollection& triangles, size_t& ring_id, vec1i& ring_ids) {

  float dupe_threshold = (float) std::pow(10,-dupe_threshold_exp);
  if (is_degenerate(poly, dupe_threshold)) {
    LinearRing new_poly = fix_duplicates(poly, dupe_threshold);
    if(is_degenerate(new_poly, dupe_threshold)) {
      std::cout << "skipping ring with duplicates\n";
      // dupe_rings.push_back(poly);
      return;
    }
    std::cout << "fixed ring with duplicates\n";
    poly = new_poly;
  }
  auto normal = calculate_normal(poly);
  if (std::isnan(normal.x) || std::isnan(normal.y) || std::isnan(normal.z)){
    std::cout << "degenerate normal: " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
    return;
  }
  auto& p0 = poly[0];
  Plane_3 plane(tri::K::Point_3(p0[0], p0[1], p0[2]), tri::K::Vector_3(normal.x, normal.y, normal.z));
  
  // project and triangulate
  tri::CDT triangulation;
  // Polygon_2 poly_2d = project(poly, plane);
  // if(CGAL::abs(poly_2d.area())<1E-4) {
  //   return;
  // }
  project_and_insert(poly, plane, triangulation);
  // triangulation.insert_constraint(poly_2d.vertices_begin(), poly_2d.vertices_end(), true);
  for (auto& ring : poly.interior_rings()) {
    project_and_insert(ring, plane, triangulation);
    // poly_2d = project(poly, plane);
    // triangulation.insert_constraint(poly_2d.vertices_begin(), poly_2d.vertices_end(), true);
  }

  if (triangulation.number_of_faces()==0)
    return;

  mark_domains(triangulation);

  // for (auto& e : triangulation.finite_edges()) {
  //   auto source = e.first->vertex(triangulation.cw(e.second))->info().point;
  //   auto target = e.first->vertex(triangulation.ccw(e.second))->info().point;
    // edges.push_back({
    //   arr3f{source},
    //   arr3f{target}
    // });
    // bool constr = triangulation.is_constrained(e);
    // edges_constr.push_back(constr);
    // edges_constr.push_back(constr);
  // }

  for (tri::CDT::Finite_faces_iterator fit = triangulation.finite_faces_begin();
  fit != triangulation.finite_faces_end(); ++fit) {

    if (!output_all_triangles && !fit->info().in_domain()) continue;

    Triangle triangle;
    triangle = {
      fit->vertex(0)->info().point, 
      fit->vertex(1)->info().point, 
      fit->vertex(2)->info().point
    };
    for (size_t j = 0; j < 3; ++j)
    {
      normals.push_back({normal.x, normal.y, normal.z});
      // values_out.push_back(values_in[vi]);
      ring_ids.push_back(ring_id);
      // nesting_levels.push_back(fit->info().nesting_level);
    }
    triangles.push_back(triangle);
  }
}

void PolygonTriangulatorNode::process()
{
  auto &rings = vector_input("polygons");
  // const auto &values_in = input("valuesf").get<vec1f>();
  typedef uint32_t N;

  TriangleCollection triangles;
  MultiTriangleCollection multitrianglecols;
  vec3f normals;
  // vec1f values_out;
  vec1i ring_ids;
  auto& volumes = output("volumes");
  // vec1i nesting_levels;
  // SegmentCollection edges;
  // vec1i edges_constr;
  // size_t vi = 0;
  // auto& dupe_rings = vector_output("dupe_rings");
  if (rings.is_connected_type(typeid(LinearRing))) {
    for (size_t ri = 0; ri < rings.size(); ++ri)
    {
      auto poly_3d = rings.get<LinearRing>(ri);
      TriangleCollection tc;
      triangulate_polygon(poly_3d, normals, tc, ri, ring_ids);
      triangles.insert(triangles.end(), tc.begin(), tc.end());
    }
    volumes.push_back((float)calculate_volume(triangles));
    multitrianglecols.push_back(triangles);
    output("multi_triangle_collections").set(multitrianglecols);
  } else if (rings.is_connected_type(typeid(std::unordered_map<int, Mesh>))) {
    // We are processing a building part here. We get a building part when we
    // cut a footprint into parts because of cutting off the underground part.
    for (size_t mi = 0; mi < rings.size(); ++mi) {
      auto meshmap = rings.get<std::unordered_map<int, Mesh>>(mi);
      double volume_sum = 0;
      for(auto& [sid, mesh] : meshmap) {
        TriangleCollection mesh_triangles;
        AttributeMap mesh_attributes;
        std::vector<attribute_value> tri_labels;
        for (size_t ri = 0; ri<mesh.get_polygons().size(); ++ri) {
          TriangleCollection tc;
          triangulate_polygon(mesh.get_polygons()[ri], normals, tc, ri, ring_ids);
          int poly_label = mesh.get_labels()[ri];
          // Need to get a label for each triangle that was generated
          for (size_t i = 0; i<tc.size(); i++) tri_labels.emplace_back(poly_label);
          triangles.insert(triangles.end(), tc.begin(), tc.end());
          mesh_triangles.insert(mesh_triangles.end(), tc.begin(), tc.end());
        }
        mesh_attributes["labels"] = tri_labels;
        multitrianglecols.push_back(mesh_triangles);
        multitrianglecols.push_back(mesh_attributes);
        multitrianglecols.building_part_ids_.push_back(sid);
        volume_sum += calculate_volume(mesh_triangles);
      }
      volumes.push_back((float)volume_sum);
    }
    output("multi_triangle_collections").set(multitrianglecols);
  } else if (rings.is_connected_type(typeid(Mesh))) {
    if(output_mtc_for_every_input) {
      multitrianglecols.get_tricollections().clear();
      multitrianglecols.get_attributes().clear();
      multitrianglecols.building_part_ids_.clear();
    }
    for (size_t mi = 0; mi < rings.size(); ++mi) {
      auto mesh = rings.get<Mesh>(mi);
      TriangleCollection mesh_triangles;
      AttributeMap mesh_attributes;
      std::vector<attribute_value> tri_labels;
      for (size_t ri = 0; ri<mesh.get_polygons().size(); ++ri) {
        TriangleCollection tc;
        triangulate_polygon(mesh.get_polygons()[ri], normals, tc, ri, ring_ids);
        int poly_label = mesh.get_labels()[ri];
        // Need to get a label for each triangle that was generated
        for (size_t i = 0; i<tc.size(); i++) tri_labels.emplace_back(poly_label);
        triangles.insert(triangles.end(), tc.begin(), tc.end());
        mesh_triangles.insert(mesh_triangles.end(), tc.begin(), tc.end());
      }
      volumes.push_back((float)calculate_volume(mesh_triangles));
      mesh_attributes["labels"] = tri_labels;
      multitrianglecols.push_back(mesh_triangles);
      multitrianglecols.push_back(mesh_attributes);
      multitrianglecols.building_part_ids_.push_back(mi);
      if(output_mtc_for_every_input) {
        output("multi_triangle_collections").push_back(multitrianglecols);
      }
    }
    if(!output_mtc_for_every_input ||
      output("multi_triangle_collections").size()==0) {
      output("multi_triangle_collections").set(multitrianglecols);
    }
  }

  // set outputs
  output("triangles").set(triangles);
  output("normals").set(normals);
  output("ring_ids").set(ring_ids);
  // output("nesting_levels").set(nesting_levels);
  // output("edges").set(edges);
  // output("edges_constr").set(edges_constr);
  // output("valuesf").set(values_out);
}

void MTC2MMNode::process()
{
  auto& mtc = input("multi_triangle_collections").get<MultiTriangleCollection&>();
  std::unordered_map<int, Mesh> meshmap;

  auto& attrmap = mtc.get_attributes();
  size_t i=0;
  for(auto& tc : mtc.get_tricollections()) {
    Mesh mesh;
    size_t j=0;
    for (auto& triangle : tc) {
      LinearRing lr;
      lr.push_back(triangle[0]);
      lr.push_back(triangle[1]);
      lr.push_back(triangle[2]);
      mesh.push_polygon( lr, std::get<int>( attrmap[i]["labels"][j++] ) );
    }
    
    meshmap[mtc.building_part_ids_[i++]] = mesh;
  }
  output("meshmap").set(meshmap);

}

} // namespace geoflow::nodes::stepedge
