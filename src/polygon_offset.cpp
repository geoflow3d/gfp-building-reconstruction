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
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_2.h>
#include<CGAL/create_offset_polygons_2.h>
#include<CGAL/create_straight_skeleton_2.h>

#include "stepedge_nodes.hpp"


namespace geoflow::nodes::stepedge {
    
    void PolygonOffsetNode::process() {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;
        typedef K::FT FT;
        typedef K::Point_2                   Point ;
        typedef CGAL::Polygon_2<K>           Polygon_2 ;
        typedef CGAL::Straight_skeleton_2<K> Ss ;
        auto& ipolygons = vector_input("polygons");
        auto& opolygons = vector_output("offset_polygons");

        // offset 0 will crash the cgal algorithm
        if (offset==0) {
            opolygons = ipolygons.get_data_vec();
            return;
        }

        for (size_t i=0; i<ipolygons.size(); ++i) {
            auto& lr = ipolygons.get<LinearRing>(i);
            Polygon_2 poly;
            for (auto&p : lr) {
                poly.push_back(Point(p[0], p[1]));
            }
            if (poly.is_clockwise_oriented()){
                poly.reverse_orientation();
            }
            auto oss = CGAL::create_exterior_straight_skeleton_2(FT(offset), poly);
            auto offset_polygons = CGAL::create_offset_polygons_2<Polygon_2>(FT(offset), *oss);
            
            // std::cout << "there are " << offset_polygons.size() << " offset polygons\n";
            // for (auto offset_poly : offset_polygons) {
            // the first offset_polygons appears to be some kind of bounding box, the second (and last) is the one we are looking for (at least it was in my test data)
            auto offset_poly = *offset_polygons.rbegin();
            LinearRing lr_offset;
            for (auto& p : *offset_poly) {
                lr_offset.push_back(arr3f{float(p.x()), float(p.y()), 0});
            }
            // }
            opolygons.push_back(lr_offset);
        }
    }

} //namespace geoflow::nodes::stepedge