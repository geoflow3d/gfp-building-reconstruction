#include <geoflow/geoflow.hpp>

namespace geoflow::nodes::stepedge {
  bool has_duplicates_ring(vec3f& poly, float& dupe_threshold);
  bool is_degenerate(LinearRing& poly, float& dupe_threshold);
}