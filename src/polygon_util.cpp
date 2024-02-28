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


  void fix_duplicates_ring(vec3f& poly, vec3f& new_ring, float& dupe_threshold) {
    auto pl = *poly.rbegin();
    for (auto& p : poly) {
      if (!(std::fabs(pl[0]-p[0])< dupe_threshold && std::fabs(pl[1]-p[1])< dupe_threshold && std::fabs(pl[2]-p[2])< dupe_threshold)) {
        new_ring.push_back(p);
      }
      pl=p;
    }
  }

  LinearRing fix_duplicates(LinearRing& poly, float& dupe_threshold) {
    LinearRing new_lr;
    fix_duplicates_ring(poly, new_lr, dupe_threshold);

    for (auto& ring : poly.interior_rings()) {
      vec3f new_ring;
      fix_duplicates_ring(ring, new_ring, dupe_threshold);
      new_lr.interior_rings().push_back(new_ring);
    }
    return new_lr;
  }
}