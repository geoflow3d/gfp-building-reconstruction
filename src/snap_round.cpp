#include "stepedge_nodes.hpp"

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>


namespace geoflow::nodes::stepedge {

  void SnapRoundNode::process() {
    typedef CGAL::Quotient<CGAL::MP_Float>           Number_type;
    typedef CGAL::Cartesian<Number_type>             Kernel;
    typedef CGAL::Snap_rounding_traits_2<Kernel>     Traits;
    typedef Kernel::Segment_2                        Segment_2;
    typedef Kernel::Point_2                          Point_2;
    typedef std::list<Segment_2>                     Segment_list_2;
    typedef std::list<Point_2>                       Polyline_2;
    typedef std::list<Polyline_2>                    Polyline_list_2;

    auto arr = input("arrangement").get<Arrangement_2>();
    
    Segment_list_2 seg_list;
    Polyline_list_2 output_list;
    for (auto he : arr.edge_handles()) {
      auto& s = he->source()->point();
      auto& t = he->target()->point();
      seg_list.push_back (
        Segment_2(
          Point_2(CGAL::to_double(s.x()), CGAL::to_double(s.y())), 
          Point_2(CGAL::to_double(t.x()), CGAL::to_double(t.y())))
      );
    }

    CGAL::snap_rounding_2<Traits,Segment_list_2::const_iterator,Polyline_list_2>
      (seg_list.begin(), seg_list.end(), output_list, pixel_size, true, false, 1);

    int counter = 0;
    Polyline_list_2::const_iterator iter1;
    for (iter1 = output_list.begin(); iter1 != output_list.end(); ++iter1) {
      std::cout << "Polyline number " << ++counter << ":\n";
      Polyline_2::const_iterator iter2;
      for (iter2 = iter1->begin(); iter2 != iter1->end(); ++iter2)
        std::cout << "    (" << iter2->x() << ":" << iter2->y() << ")\n";
    }
    
  }
}