#include "stepedge_nodes.hpp"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <list>

typedef CGAL::Simple_cartesian<double> K;

// custom triangle type with
// three pointers to points
struct My_triangle {
    geoflow::Triangle m_t;
    size_t m_face_id;
    My_triangle(geoflow::Triangle t, size_t face_id)
        : m_t(t), m_face_id(face_id) {}
};
// the custom triangles are stored into a vector
typedef std::list<My_triangle>::const_iterator Iterator;
// The following primitive provides the conversion facilities between
// the custom triangle and point types and the CGAL ones
struct My_triangle_primitive {
public:
    // this is the type of data that the queries returns. For this example
    // we imagine that, for some reasons, we do not want to store the iterators
    // of the vector, but raw pointers. This is to show that the Id type
    // does not have to be the same as the one of the input parameter of the 
    // constructor.
    typedef const My_triangle* Id;
    // CGAL types returned
    typedef K::Point_3    Point; // CGAL 3D point type
    typedef K::Triangle_3 Datum; // CGAL 3D triangle type
private:
    Id m_mytri; // this is what the AABB tree stores internally
public:
    My_triangle_primitive() {} // default constructor needed
    // the following constructor is the one that receives the iterators from the 
    // iterator range given as input to the AABB_tree
    My_triangle_primitive(Iterator it)
        : m_mytri(&(*it)) {}
    const Id& id() const { return m_mytri; }
    // on the fly conversion from the internal data to the CGAL types
    Point convert(const geoflow::arr3f& p) const {
      return Point(p[0], p[1], p[2]);
    }
    Datum datum() const {
        return Datum(
          convert(m_mytri->m_t[0]),
          convert(m_mytri->m_t[1]),
          convert(m_mytri->m_t[2])
        );
    }
    // returns a reference point which must be on the primitive
    Point reference_point() const { 
      return convert(m_mytri->m_t[0]);
    }
};

namespace geoflow::nodes::stepedge {
    
  void PC2MeshQualityNode::process() {
    typedef CGAL::AABB_traits<K, My_triangle_primitive> My_AABB_traits;
    typedef CGAL::AABB_tree<My_AABB_traits> Tree;

    auto& points = input("ipoints").get<PointCollection>();
    auto& trin = input("triangles").get<TriangleCollection>();
    auto& face_ids = input("face_ids").get<vec1i>();
    std::list<My_triangle> triangles;
    for (size_t i=0; i < trin.size(); ++i) {
      triangles.push_back(My_triangle(trin[i], face_ids[i*3]));
    }

    // constructs AABB tree
    Tree tree(triangles.begin(), triangles.end());
    tree.accelerate_distance_queries();

    std::map<size_t, std::vector<double>> distances;
    vec1f point_errors;
    size_t i = 0;
    for(auto& p_ : points) {
      auto p = K::Point_3(p_[0], p_[1], p_[2]);
      auto pt_and_id = tree.closest_point_and_primitive(p);
      auto sqd = CGAL::squared_distance(pt_and_id.first, p);
      auto fid = pt_and_id.second->m_face_id;
      distances[fid].push_back(sqd);
      point_errors.push_back(sqd);
    }

    double sum_total = 0;
    size_t len = 0;
    std::map<size_t, float> face_error_map;
    for (auto& [fid, errors] : distances) {
      double sum_face = 0;
      len += errors.size();
      for(double& error : errors) {
        sum_face += error;
      }
      face_error_map[fid] = float(CGAL::sqrt(sum_face/errors.size()));
      sum_total += sum_face;
    }
    float rms_error = float(CGAL::sqrt(sum_total/len));

    vec1f face_errors;
    for (size_t i=0; i < trin.size(); ++i) {
      size_t fid = face_ids[i*3];
      // push zero error if this face has no points/no error defined
      if(face_error_map.find(fid) == face_error_map.end()) {
        face_errors.push_back(0);
        face_errors.push_back(0);
        face_errors.push_back(0);
      } else {
        face_errors.push_back(face_error_map[fid]);
        face_errors.push_back(face_error_map[fid]);
        face_errors.push_back(face_error_map[fid]);
      }
    }

    vec1f mesh_error;
    for (auto& t : triangles) {
      mesh_error.push_back(rms_error);
      mesh_error.push_back(rms_error);
      mesh_error.push_back(rms_error);
    }

    output("point_errors").set(point_errors);
    output("face_errors").set(face_errors);
    output("mesh_error_f").set(rms_error);
    output("mesh_error").set(mesh_error);
  }
}