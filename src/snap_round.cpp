#include "stepedge_nodes.hpp"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>


namespace geoflow::nodes::stepedge {

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

    output("arrangement").set(arrSR);
    
  }
}