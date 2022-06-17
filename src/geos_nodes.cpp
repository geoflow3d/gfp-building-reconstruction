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
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/algorithm/Orientation.h>

#include <geos/simplify/DouglasPeuckerSimplifier.h>

#include "stepedge_nodes.hpp"


namespace geoflow::nodes::stepedge {

    geos::geom::GeometryFactory::Ptr geos_global_factory;

    std::unique_ptr<geos::geom::Polygon> to_geos_polygon(const LinearRing& lr) {
        auto coordSeq = std::make_unique<geos::geom::CoordinateArraySequence>();
        for (auto&p : lr) {
            coordSeq->add(geos::geom::Coordinate(p[0], p[1]));
        }
        coordSeq->add(geos::geom::Coordinate(lr[0][0], lr[0][1]));
        if (!geos::algorithm::Orientation::isCCW(coordSeq.get())){
            geos::geom::CoordinateArraySequence::reverse(coordSeq.get());
        }
        auto linearRing = geos_global_factory->createLinearRing(std::move(coordSeq));
        return geos_global_factory->createPolygon(std::move(linearRing));
    }

    void PolygonSimplifyGEOSNode::process() {
        auto& ipolygons = vector_input("polygons");
        auto& opolygons = vector_output("simplified_polygons");

        geos_global_factory = geos::geom::GeometryFactory::create();

        for (size_t i=0; i<ipolygons.size(); ++i) {
            auto& lr = ipolygons.get<LinearRing>(i);
            
            auto polygon = to_geos_polygon(lr);

            auto simplified_geom = geos::simplify::DouglasPeuckerSimplifier::simplify(polygon.get(), tolerance).release();
            // auto buf_geom = polygon->buffer(offset).release();

            if(!simplified_geom->isValid()) {
                std::cout << "feature not simplified\n";
                opolygons.push_back(lr);
            } else if (auto buf_poly = dynamic_cast<geos::geom::Polygon*>(simplified_geom)) {
                auto polygon_ring = buf_poly->getExteriorRing();

                LinearRing lr_offset;
                for (size_t i=0; i<polygon_ring->getNumPoints()-1; ++i) {
                    auto& p = polygon_ring->getCoordinateN(i);
                    lr_offset.push_back(arr3f{float(p.x), float(p.y), 0});
                }
                opolygons.push_back(lr_offset);                
            } else {
                std::cout << "feature not simplified\n";
                opolygons.push_back(lr);
            }

        }
    }

    void PolygonBufferGEOSNode::process() {
        auto& ipolygons = vector_input("polygons");
        auto& opolygons = vector_output("offset_polygons");

        geos_global_factory = geos::geom::GeometryFactory::create();

        for (size_t i=0; i<ipolygons.size(); ++i) {
            auto& lr = ipolygons.get<LinearRing>(i);
            
            auto polygon = to_geos_polygon(lr);

            auto buf_geom = polygon->buffer(offset).release();
            
            if(!buf_geom->isValid()) {
                std::cout << "feature not buffered\n";
                opolygons.push_back(lr);
            } else if (auto buf_poly = dynamic_cast<geos::geom::Polygon*>(buf_geom)) {
                auto polygon_ring = buf_poly->getExteriorRing();

                LinearRing lr_offset;
                for (size_t i=0; i<polygon_ring->getNumPoints()-1; ++i) {
                    auto& p = polygon_ring->getCoordinateN(i);
                    lr_offset.push_back(arr3f{float(p.x), float(p.y), 0});
                }
                // }
                opolygons.push_back(lr_offset);
            } else {
                std::cout << "feature not buffered\n";
                opolygons.push_back(lr);
            }
            delete buf_geom;
        }
    }

} //namespace geoflow::nodes::stepedge