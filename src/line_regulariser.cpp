#include "line_regulariser.hpp"
#include <iterator>

namespace linereg {

  double calc_mean_angle(const std::vector<linetype*>& lines) {
    // length-weighted mean angle
    double angle_sum=0, sqlenth_sum=0;
    for (auto& line : lines) {
      sqlenth_sum += line->sqlength;
    }
    for (auto& line : lines) {
      angle_sum += line->angle * line->sqlength;
    }
    return angle_sum/sqlenth_sum;
  }

  Point_2 calc_centroid(const std::vector<linetype*>& lines) {
    double cx=0, cy=0;
    for (auto& line : lines) {
      auto p = line->segment.source();
      cx += CGAL::to_double(p.x());
      cy += CGAL::to_double(p.y());
      p = line->segment.target();
      cx += CGAL::to_double(p.x());
      cy += CGAL::to_double(p.y());
    }
    size_t np = 2*lines.size();
    return Point_2(cx/np, cy/np);
  }

  // construct a line L with centroid and mean_angle, then project all the segments from lines on L to bound it to a segment
  Segment_2 calc_segment(Point_2 centroid, double mean_angle, const std::vector<linetype*>& lines, double extension) {
    auto lv = Vector_2(std::tan(mean_angle), 1.0);
    lv = lv / CGAL::sqrt(CGAL::to_double(lv.squared_length()));

    bool setminmax=false;
    Point_2 pmin, pmax;
    double dmin, dmax;
    for (auto& line : lines) {
      auto p = line->segment.source();
      auto d = CGAL::to_double(Vector_2(p.x(),p.y())*lv);
      if (!setminmax) {
        setminmax=true;
        dmin=dmax=d;
        pmin=pmax=p;
      }
      if (d < dmin){
        dmin = d;
        pmin = p;
      }
      if (d > dmax) {
        dmax = d;
        pmax = p;
      }

      p = line->segment.target();
      d = CGAL::to_double(Vector_2(p.x(),p.y())*lv);
      if (d < dmin){
        dmin = d;
        pmin = p;
      }
      if (d > dmax) {
        dmax = d;
        pmax = p;
      }
    }
    Line_2 l(centroid, lv);
    pmin = l.projection(pmin) - lv*extension;
    pmax = l.projection(pmax) + lv*extension;
    return Segment_2(pmin, pmax);
  }

  double AngleCluster::distance(Cluster<double>* other_cluster) {
    return std::fabs(value - other_cluster->value);
  }
  void AngleCluster::calc_mean_value() {
    value = calc_mean_angle(lines);
  }

  double DistCluster::distance(Cluster<Segment_2>* other_cluster) {
    return CGAL::to_double(CGAL::squared_distance(value, other_cluster->value));
  }
  void DistCluster::calc_mean_value() {
    // a segment through the length-weighted mean centroid and elongated to 'cover' all the segments
    double mean_angle = calc_mean_angle(lines);
    Point_2 centroid = calc_centroid(lines);
    value = calc_segment(centroid, mean_angle, lines);
  }

  template <typename ClusterH> DistanceTable<ClusterH>::DistanceTable(std::set<ClusterH>& clusters) : clusters(clusters) {
    // compute only half of the distance table, since the other half will be exactly the same
    for(auto cit_a = clusters.begin(); cit_a != clusters.end(); ++cit_a) {
      for(auto cit_b = std::next(cit_a); cit_b != clusters.end(); ++cit_b) {
        auto cluster_a = *cit_a;
        auto cluster_b = *cit_b;
        distances[std::make_pair(cluster_a, cluster_b)] = cluster_a->distance(cluster_b.get());
      }
    }
  }
  template <typename ClusterH> void DistanceTable<ClusterH>::merge(ClusterH lhs, ClusterH rhs) {
    // merges two clusters, then removes one from the distances map and update the affected distances
    // merge
    lhs->lines.insert(lhs->lines.end(), rhs->lines.begin(), rhs->lines.end());
    lhs->calc_mean_value();

    clusters.erase(rhs);
    // iterate distances
    // if rhs in pair remove element
    // if lhs in pair update dist
    std::vector<typename DistanceMap::iterator> to_remove;
    for (auto it = distances.begin(); it != distances.end(); ++it) {
      if (it->first.first == rhs || it->first.second == rhs) {
        to_remove.push_back(it);
      }
    }
    for (auto& it : to_remove) {
      distances.erase(it);
    }
    for (auto it = distances.begin(); it != distances.end(); ++it) {
      if (it->first.first == lhs || it->first.second == lhs) {
        it->second = it->first.first->distance(it->first.second.get());
      }
    }
  }
  //NB: this will crash if there is only one cluster!!
  template <typename ClusterH> std::pair<typename DistanceTable<ClusterH>::ClusterPair, double>  DistanceTable<ClusterH>::get_closest_pair() {
    auto min_it = distances.begin();
    for(auto it = distances.begin(); it != distances.end(); ++it) {
      if(it->second < min_it->second) {
        min_it = it;
      }
    }
    return *min_it;
  }

  template class DistanceTable<AngleClusterH>;
  template class DistanceTable<DistClusterH>;

  void LineRegulariser::perform_angle_clustering() {
    //make clusters
    angle_clusters.clear();
    for(auto& line : lines) {
      auto aclusterh = std::make_shared<AngleCluster>();
      aclusterh->value = line.angle;
      aclusterh->lines.push_back(&line);
      angle_clusters.insert(aclusterh);
    }
    // make distance table
    DistanceTable adt(angle_clusters);

    if (angle_clusters.size()>1) {
      auto apair = adt.get_closest_pair();
      while (apair.second < angle_threshold) {
        adt.merge(apair.first.first, apair.first.second);
        if (angle_clusters.size()==1) break;
        apair = adt.get_closest_pair();
      }
    }

    size_t id_cntr=0;
    for(auto& aclusterh : angle_clusters) {
      for(auto line : aclusterh->lines){
        line->angle_cluster = aclusterh;
        line->angle_cluster_id = id_cntr;
      }
      ++id_cntr;
    }
  }
  void LineRegulariser::perform_distance_clustering() {
    dist_clusters.clear();

    // perform distance clustering for each angle cluster
    for(auto& aclusterh : angle_clusters) {
      std::set<DistClusterH> dclusters;
      for(auto& line : aclusterh->lines) {
        auto dclusterh = std::make_shared<DistCluster>();
        dclusterh->value = line->segment;
        dclusterh->lines.push_back(line);
        dclusters.insert(dclusterh);
      }
      
      if (dclusters.size()>1) {
        // make distance table
        DistanceTable ddt(dclusters);

        // do clustering
        auto dpair = ddt.get_closest_pair();
        while (dpair.second < dist_threshold) {
          ddt.merge(dpair.first.first, dpair.first.second);
          if (dclusters.size()==1) break;
          dpair = ddt.get_closest_pair();
        }
      }
      dist_clusters.insert(dclusters.begin(), dclusters.end());
    }

    size_t id_cntr=0;
    for(auto& dclusterh : dist_clusters) {
      for(auto line : dclusterh->lines) {
        line->dist_cluster = dclusterh;
        line->dist_cluster_id = id_cntr;
      }
      ++id_cntr;
    }
  }
}