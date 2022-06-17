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
#include "stepedge_nodes.hpp"
#include "polygon_util.hpp"
#include "cdt_util.hpp"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>

#include <CGAL/Arr_walk_along_line_point_location.h>



namespace geoflow::nodes::stepedge {

  typedef tri_util::CDT CDT;

  void SnapRoundNode::process() {
    typedef CGAL::Quotient<CGAL::MP_Float>           Number_type;
    typedef CGAL::Cartesian<Number_type>             QKernel;
    typedef CGAL::Snap_rounding_traits_2<QKernel>    QTraits;
    typedef std::list<QTraits::Segment_2>            Segment_list_2;
    typedef std::list<QTraits::Point_2>              Polyline_2;
    typedef std::list<Polyline_2>                    Polyline_list_2;

    typedef CGAL::Arr_segment_traits_2<QKernel>                      Traits_2;
    typedef Traits_2::X_monotone_curve_2                            ARSegment_2;
    typedef CGAL::Arrangement_2<Traits_2>                           QArrangement_2;

    typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;

    struct overlay_functor {
      FaceInfo operator()(const FaceInfo& a, FaceInfo& b) { 
        b=a;
        return a;
      }
    };
    typedef CGAL::Arr_face_overlay_traits<Arrangement_2,
                                          Arrangement_2,
                                          Arrangement_2,
                                          overlay_functor >  Overlay_traits;

    auto arr = input("arrangement").get<Arrangement_2>();
    
    Segment_list_2 seg_list;
    Polyline_list_2 output_list;
    for (auto he : arr.edge_handles()) {
      auto& s = he->source()->point();
      auto& t = he->target()->point();
      seg_list.push_back (
        QTraits::Segment_2(
          QTraits::Point_2(CGAL::to_double(s.x()), CGAL::to_double(s.y())), 
          QTraits::Point_2(CGAL::to_double(t.x()), CGAL::to_double(t.y())))
      );
    }

    CGAL::snap_rounding_2<QTraits,Segment_list_2::const_iterator,Polyline_list_2>
      (seg_list.begin(), seg_list.end(), output_list, pixel_size, true, false, 1);

    // QArrangement_2 arrSR;
    Arrangement_2 arrSR;
    // int counter = 0;
    Polyline_list_2::const_iterator polyline_it;
    for (polyline_it = output_list.begin(); polyline_it != output_list.end(); ++polyline_it) {
      // std::cout << "Polyline number " << ++counter << ":\n";

      if (polyline_it->size() > 1) {
        Polyline_2::const_iterator p = polyline_it->begin();
        Polyline_2::const_iterator p_prev = p;
        for (++p; p != polyline_it->end(); ++p) {
          // std::cout << "    (" << CGAL::to_double(p->x()) << ":" << CGAL::to_double(p->y()) << ")\n";
          insert(arrSR, Segment_2(
            Point_2(CGAL::to_double(p->x()), CGAL::to_double(p->y())), 
            Point_2(CGAL::to_double(p_prev->x()), CGAL::to_double(p_prev->y()))
          ));
          p_prev = p;
        }
      }
    }
    std::cout << "vcount before: " << arr.number_of_vertices() << "\n";
    std::cout << "vcount after:  " << arrSR.number_of_vertices() << "\n";

    Walk_pl walk_pl (arr);
    for (auto& arrFace : arrSR.face_handles()) {
      LinearRing poly;
      if(arrFace->is_fictitious()) continue;
      
      if (arrangementface_to_polygon(arrFace, poly)) {

        CDT cdt = tri_util::create_from_polygon(poly);
        for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
          fit != cdt.finite_faces_end(); ++fit) {

          if (!fit->info().in_domain()) continue;
          auto p = CGAL::centroid(cdt.triangle(fit));
          auto obj = walk_pl.locate( Walk_pl::Arrangement_2::Point_2(p.x(), p.y()) );

          // std::cout << "The point (" << p << ") is located ";
          if (auto f = boost::get<Face_const_handle>(&obj)) { // located inside a face
            // std::cout << "inside "
            //           << (((*f)->is_unbounded()) ? "the unbounded" : "a bounded")
            //           << " face." << std::endl;
            arrFace->data() = (*f)->data();
          }
          // else if (auto e = boost::get<Halfedge_const_handle>(&obj)) // located on an edge
          //   std::cout << "on an edge: " << (*e)->curve() << std::endl;
          // else if (auto v = boost::get<Vertex_const_handle>(&obj)) // located on a vertex
          //   std::cout << "on " << (((*v)->is_isolated()) ? "an isolated" : "a")
          //             << " vertex: " << (*v)->point() << std::endl;
          // else CGAL_error_msg("Invalid object.");
        }
      }

    }

    // Arrangement_2 arr_overlay;
    // Overlay_traits overlay_traits;
    // overlay(arr, arrSR, arr_overlay, overlay_traits);

    output("arrangement").set(arrSR);
    
  }
}