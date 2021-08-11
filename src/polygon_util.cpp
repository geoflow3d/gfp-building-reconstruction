#include "polygon_util.hpp"

namespace geoflow::nodes::stepedge {
  bool has_duplicates_ring(vec3f& poly, float& dupe_threshold) {
    auto pl = *poly.rbegin();
    for (auto& p : poly) {
      if (std::fabs(pl[0]-p[0])< dupe_threshold && std::fabs(pl[1]-p[1])< dupe_threshold && std::fabs(pl[2]-p[2])< dupe_threshold) {
        return true;
      }
      pl=p;
    }
    return false;
  }

  bool is_degenerate(LinearRing& poly, float& dupe_threshold) {
    if (poly.size() < 3) return true;
    if (has_duplicates_ring(poly, dupe_threshold)) return true;

    for (auto& ring : poly.interior_rings()) {
      if (ring.size() < 3) return true;
      if (has_duplicates_ring(ring, dupe_threshold)) return true;
    }
    return false;
  }
}