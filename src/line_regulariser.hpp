#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <geoflow/geoflow.hpp>

namespace linereg {
  typedef CGAL::Exact_predicates_exact_constructions_kernel K;
  // typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  typedef K::Point_2 Point_2;
  typedef K::Point_3 Point_3;
  typedef K::Vector_2 Vector_2;
  typedef K::Line_2 Line_2;
  typedef K::Segment_2 Segment_2;
  typedef CGAL::Polygon_2<K> Polygon_2;
  typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

  struct linetype;

  template <typename T> struct Cluster {
    T value;
    bool has_intersection_line;
    std::vector<linetype*> lines;
    virtual double distance(Cluster<T>* other_cluster)=0;
    virtual void calc_mean_value()=0;
  };
  struct AngleCluster : public Cluster<double>{
    double distance(Cluster<double>* other_cluster);
    void calc_mean_value();
  };
  struct DistCluster : public Cluster<Segment_2>{
    double distance(Cluster<Segment_2>* other_cluster);
    void calc_mean_value();
  };

  template <typename ClusterH> struct DistanceTable {
    typedef std::pair<ClusterH, ClusterH> ClusterPair;

    // define has function such that the same hash results regardless of the order of cluster handles in the pair
    struct KeyHash {
      size_t operator()(const ClusterPair& key) const {
        if (key.first.get() < key.second.get())
          return std::hash<ClusterH>()(key.first) ^
            (std::hash<ClusterH>()(key.second) << 1);
        else
          return std::hash<ClusterH>()(key.second) ^
            (std::hash<ClusterH>()(key.first) << 1);

      }
    };
    // True equality function is needed to deal with collisions
    struct KeyEqual {
      bool operator()(const ClusterPair& lhs, const ClusterPair& rhs) const {
        return 
          ((lhs.first==rhs.first) && (lhs.second==rhs.second)) 
          ||
          ((lhs.first==rhs.second) && (lhs.second==rhs.first));
      }
    };
    typedef std::unordered_map<ClusterPair, double, KeyHash, KeyEqual> DistanceMap;

    DistanceMap distances;
    std::set<ClusterH>& clusters;

    DistanceTable(std::set<ClusterH>& clusters); //computes initial distances
    void merge(ClusterH lhs, ClusterH rhs); // merges two clusters, then removes one from the distances map and update the affected distances
    typename DistanceMap::iterator get_closest_pair(); //returns the cluster pair with the smallest distance
  };
  
  // extern template class Cluster<double>;
  // extern template class Cluster<Segment_2>;

  typedef std::shared_ptr<AngleCluster> AngleClusterH;
  typedef std::shared_ptr<DistCluster> DistClusterH;

  extern template class DistanceTable<AngleClusterH>;
  extern template class DistanceTable<DistClusterH>;

  double calc_mean_angle(const std::vector<linetype*>& lines);
  Point_2 calc_centroid(const std::vector<linetype*>& lines);
  Segment_2 calc_segment(Point_2 centroid, double mean_angle, const std::vector<linetype*>& lines, double extension=0);

  struct linetype {
    linetype(Segment_2 segment_, double angle_, Point_2 midpoint_, double dist_in_ang_cluster_, size_t priority_, size_t segment_id_, double sqlength_) :
    segment(segment_), angle(angle_), midpoint(midpoint_), priority(priority_), segment_id(segment_id_), sqlength(sqlength_) {};
    
    Segment_2 segment;
    double angle;
    Point_2 midpoint;
    Vector_2 direction;
    double sqlength;
    // double dist_in_ang_cluster;
    size_t priority;
    size_t segment_id;

    Segment_2 reg_segment;
    AngleClusterH angle_cluster;
    DistClusterH dist_cluster;
    size_t angle_cluster_id, dist_cluster_id;
  };

  static constexpr double pi = 3.14159265358979323846;
  class LineRegulariser {

    typedef std::vector<Segment_2> SegmentVec;

    // geoflow::SegmentCollection& input_segments;
    public:
    std::vector<linetype> lines;
    // SegmentVec input_reg_exact;
    double angle_threshold, dist_threshold;

    std::unordered_map<size_t, SegmentVec> segments;
    std::set<AngleClusterH> angle_clusters;
    std::set<DistClusterH> dist_clusters;

    LineRegulariser() {};

    void add_segments(size_t priority, const Polygon_2& polygon, double offset) {
      size_t i;
      if(segments.find(priority) == segments.end()) {
        i=0;
      } else {
        i=segments.size();
      }
      auto orientation = polygon.orientation();
      for(auto edge = polygon.edges_begin(); edge != polygon.edges_end(); ++edge) {
        auto source = edge->source();
        auto target = edge->target();
        auto perp = (target-source).perpendicular(orientation);
        auto len = CGAL::sqrt(CGAL::to_double(perp.squared_length()));
        // std::cout << "len: " << len << "\n"; 
        perp = offset * (perp/len);
        target -= perp;
        source -= perp;
        auto v = target-source;
        auto p_ = source + v/2;
        auto p = Point_2(CGAL::to_double(p_.x()),CGAL::to_double(p_.y()));
        auto l = CGAL::to_double(v.squared_length());
        auto angle = std::atan2(CGAL::to_double(v.x()), CGAL::to_double(v.y()));
        if (angle < 0) angle += pi;
        lines.push_back(linetype(*edge, angle,p,0,priority,i++,l));
        segments[priority].push_back(Segment_2(target,source));
      }
    }

    void add_segments(size_t priority, geoflow::SegmentCollection& segs) {
      if (segs.size()==0) return;
      size_t i;
      if(segments.find(priority) == segments.end()) {
        i=0;
      } else {
        i=segments.size();
      }

      for(auto& edge : segs) {
        auto source = Point_2(edge[0][0], edge[0][1]);
        auto target = Point_2(edge[1][0], edge[1][1]);
        auto v = target-source;
        auto p_ = source + v/2;
        auto p = Point_2(CGAL::to_double(p_.x()),CGAL::to_double(p_.y()));
        auto l = CGAL::to_double(v.squared_length());
        auto angle = std::atan2(CGAL::to_double(v.x()), CGAL::to_double(v.y()));
        if (angle < 0) angle += pi;
        lines.push_back(linetype(Segment_2(source, target), angle,p,0,priority,i++,l));
        segments[priority].push_back(Segment_2(source, target));
      }
    }

    SegmentVec& get_segments(const size_t& priority) {
      return segments[priority];
    }

    void perform_angle_clustering();
    void perform_distance_clustering();
  };
  template<class Kernel> void 
  chain(
    const typename Kernel::Segment_2& a, 
    const typename Kernel::Segment_2& b, 
    typename CGAL::Polygon_2<Kernel>& ring, 
    const float& snap_threshold) {

    auto l_a = a.supporting_line();
    auto l_b = b.supporting_line();
    typename Kernel::Segment_2 s(a.target(), b.source());
    auto result = CGAL::intersection(l_a, l_b);
    if (result) {
      if (auto p = boost::get<typename Kernel::Point_2>(&*result)) {
        if (CGAL::squared_distance(*p, s) < snap_threshold) {
          ring.push_back(*p);
        } else {
          ring.push_back(a.target());
          ring.push_back(b.source());
        }
      }
    // } else if (auto l = boost::get<K::Line_2>(&*result)) {

    // }
    } else { // there is no intersection
      ring.push_back(a.target());
      ring.push_back(b.source());
    }
  }

  // void chain(Segment& a, Segment& b, LinearRing& ring, const float& snap_threshold) {
  template <class Kernel> inline void check_dist(const CGAL::Polygon_2<Kernel>& pos, CGAL::Polygon_2<Kernel>& pot, const size_t a, const size_t b) {
    auto d = CGAL::squared_distance(pos.vertex(a), pos.vertex(b));
    if (d > 1E-6) pot.push_back(pos.vertex(a));
  }
  
  template<class Kernel> CGAL::Polygon_2<Kernel> 
  chain_ring(
    const std::vector<size_t>& idx, 
    const std::vector<typename Kernel::Segment_2>& segments, 
    const float& snap_threshold) {

    typename CGAL::Polygon_2<Kernel>  ring, fixed_ring;

    if (idx.size()>1) { // we need at least 2 segments
      for (size_t i=1; i<idx.size(); ++i) {
        chain<Kernel>(segments[idx[i-1]], segments[idx[i]], ring, snap_threshold);
      }
      chain<Kernel>(segments[idx[idx.size()-1]], segments[idx[0]], ring, snap_threshold);

      // get rid of segments with zero length (ie duplicate points)
      // check again the size, to ignore degenerate case of input ring that consists of 3 co-linear segments (would get chained to eg 0 vertices)
      if (ring.size()>2) {
        for (size_t i=1; i<ring.size(); ++i) {
          check_dist<Kernel>(ring, fixed_ring, i-1, i);
        }
        check_dist<Kernel>(ring, fixed_ring, ring.size()-1, 0);
      }

      // NOTE: at this point there can still be vertices between co-linear segments (ie 3 consecutive points on the same line)
    }

    return fixed_ring;
  }
}