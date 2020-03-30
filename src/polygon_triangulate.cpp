#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_2.h>

#include <stepedge_nodes.hpp>

namespace geoflow::nodes::stepedge {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag Tag;
struct VertexInfoPIT {
  bool hasPoint;
  geoflow::arr3f point;
  VertexInfoPIT() : hasPoint(false)
  {}
};
struct FaceInfoPIT {
  // bool processed;
  int nesting_level;
  FaceInfoPIT() {
    // processed = false;
    // interior = false;
    nesting_level = -1;
  }
  bool in_domain() {
    return nesting_level % 2 == 1;
  }
};
typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfoPIT, K> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<K> FaceBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfoPIT, K, FaceBase> FaceBaseWithInfo;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBaseWithInfo> TriangulationDataStructure;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TriangulationDataStructure, Tag> CDT;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Plane_3<K> Plane_3;

glm::vec3 calculate_normal<typename T>(const T ring)
{
  glm::vec3 normal(0, 0, 0);
  for (size_t i = 0; i < ring.size(); ++i)
  {
    const auto &curr = ring[i];
    const auto &next = ring[(i + 1) % ring.size()];
    normal[0] += (curr[1] - next[1]) * (curr[2] + next[2]);
    normal[1] += (curr[2] - next[2]) * (curr[0] + next[0]);
    normal[2] += (curr[0] - next[0]) * (curr[1] + next[1]);
  }
  return glm::normalize(normal);
}

void mark_domains(CDT& ct,
  CDT::Face_handle start,
  int index,
  std::list<CDT::Edge>& border) {
  if (start->info().nesting_level != -1) {
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while (!queue.empty()) {
    CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if (fh->info().nesting_level == -1) {
      fh->info().nesting_level = index;
      for (int i = 0; i < 3; i++) {
        CDT::Edge e(fh, i);
        CDT::Face_handle n = fh->neighbor(i);
        if (n->info().nesting_level == -1) {
          if (ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}

/**
 * mark the triangles that are inside the original 2D polygon.
 * explore set of facets connected with non constrained edges,
 * and attribute to each such set a nesting level.
 * start from facets incident to the infinite vertex, with a nesting
 * level of 0. Then recursively consider the non-explored facets incident
 * to constrained edges bounding the former set and increase the nesting level by 1.
 * facets in the domain are those with an odd nesting level.
 */
void mark_domains(CDT& cdt) {
  // for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it) {
  //   it->info().nesting_level = -1;
  // }
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border);
  while (!border.empty()) {
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if (n->info().nesting_level == -1) {
      mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
    }
  }
}

void project_and_insert(geoflow::vec3f& ring, Plane_3& plane, CDT& cdt) {
  auto pit_last = ring.end()-1;
  CDT::Vertex_handle new_vh, vh_last, vh = cdt.insert(plane.to_2d(K::Point_3((*pit_last)[0], (*pit_last)[1], (*pit_last)[2])));
  vh_last = vh;
  vh->info().point = *pit_last;
  vh->info().hasPoint = true;
  for (auto pit = ring.begin(); pit != ring.end(); ++pit) {
    if(pit==pit_last){
      new_vh=vh_last;
    } else {
      new_vh = cdt.insert(plane.to_2d(K::Point_3((*pit)[0], (*pit)[1], (*pit)[2])));
      new_vh->info().point = *pit;
      new_vh->info().hasPoint = true;
    }
    cdt.insert_constraint(vh, new_vh);
    vh = new_vh;
  }
}

void PolygonTriangulatorNode::process()
{
  auto &rings = vector_input("polygons");
  // const auto &values_in = input("valuesf").get<vec1f>();
  typedef uint32_t N;

  TriangleCollection triangles;
  vec3f normals;
  vec1f values_out;
  vec1i ring_ids;
  size_t vi = 0;
  for (size_t ri = 0; ri < rings.size(); ++ri)
  {
    auto poly_3d = rings.get<LinearRing>(ri);
    if (poly_3d.size() < 3)
      continue;

    auto normal = calculate_normal(poly_3d);
    if (std::isnan(normal.x) || std::isnan(normal.y) || std::isnan(normal.z)){
      std::cout << "degenerate normal: " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
      continue;
    }
    auto& p0 = poly_3d[0];
    Plane_3 plane(K::Point_3(p0[0], p0[1], p0[2]), K::Vector_3(normal.x, normal.y, normal.z));
    
    // project and triangulate
    CDT triangulation;
    project_and_insert(poly_3d, plane, triangulation);
    for (auto& ring : poly_3d.interior_rings()) {
      project_and_insert(poly_3d, plane, triangulation);
    }

    mark_domains(triangulation);

    for (CDT::Finite_faces_iterator fit = triangulation.finite_faces_begin();
    fit != triangulation.finite_faces_end(); ++fit) {
      if (fit->info().in_domain()) {
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
          ring_ids.push_back(ri);
        }
        triangles.push_back(triangle);
      }
    }
    vi++;
  }

  // set outputs
  output("triangles").set(triangles);
  output("normals").set(normals);
  output("ring_ids").set(ring_ids);
  // output("valuesf").set(values_out);
}

}