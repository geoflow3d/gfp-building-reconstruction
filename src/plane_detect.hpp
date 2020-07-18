#pragma once

#include <geoflow/common.hpp>

#include <region_grower_DS_CGAL.hpp>
#include <region_grower.hpp>

#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Plane_3.h>

#include <deque>

namespace planedect {

typedef CGAL::Exact_predicates_inexact_constructions_kernel cgal_kernel;
typedef cgal_kernel::Point_3 Point;
typedef cgal_kernel::Vector_3 Vector;
typedef cgal_kernel::Plane_3 Plane;

class PlaneDS : public regiongrower::CGAL_RegionGrowerDS {
  public:
  geoflow::vec3f& normals;
  std::vector<Plane> seed_planes;
  
  PlaneDS(geoflow::PointCollection& points, geoflow::vec3f& normals, size_t N=15) 
    : CGAL_RegionGrowerDS(points, N), normals(normals) {};

  // Note this crashes when idx.size()==1;
  inline double fit_plane(std::vector<size_t>& idx, Plane& plane){
    std::vector<Point> neighbor_points;
    for (auto i: idx)
      neighbor_points.push_back(Point(points[i][0], points[i][1], points[i][2]));
    double quality = linear_least_squares_fitting_3(neighbor_points.begin(),neighbor_points.end(),plane,CGAL::Dimension_tag<0>());
    return quality;
  }

  virtual std::deque<size_t> get_seeds() override {
    // seed generation
    typedef std::pair<size_t,double> index_dist_pair;
    auto cmp = [](index_dist_pair left, index_dist_pair right) {return left.second < right.second;};
    std::priority_queue<index_dist_pair, std::vector<index_dist_pair>, decltype(cmp)> pq(cmp);

    size_t i=0;
    Plane plane;
    for(size_t pi=0; pi<size; ++pi){
      auto neighbours = get_neighbours(pi);
      auto quality = fit_plane(neighbours, plane);
      pq.push(index_dist_pair(i++, quality));
    }

    std::deque<size_t> seed_order;
    while (pq.size()>0) {
      seed_order.push_back(pq.top().first); pq.pop();
    }
    return seed_order;
  }
};

class PlaneRegion : public regiongrower::Region {
  public:
  using regiongrower::Region::Region;
  Plane plane;
  std::vector<size_t> inliers;
};

class DistAndNormalTester {
  public:
  float dist_thres;
  float normal_thres;
  size_t n_refit, refit_counter=0;

  DistAndNormalTester(float dist_thres=0.04, float normal_thres=0.9, size_t n_refit=5) : 
  dist_thres(dist_thres), normal_thres(normal_thres), n_refit(n_refit) {};

  bool is_valid(PlaneDS& cds, size_t candidate, size_t neighbour, PlaneRegion& shape) {
    Point p_c = Point(cds.points[candidate][0], cds.points[candidate][1], cds.points[candidate][2]);
    Vector n_c = Vector(cds.normals[candidate][0], cds.normals[candidate][1], cds.normals[candidate][2]);
    Point p = Point(cds.points[neighbour][0], cds.points[neighbour][1], cds.points[neighbour][2]);
    Vector n = Vector(cds.normals[neighbour][0], cds.normals[neighbour][1], cds.normals[neighbour][2]);

    if (shape.inliers.size()==0) {
      shape.plane = Plane(p_c,n_c);
      shape.inliers.push_back(candidate);
    }

    bool valid = (CGAL::squared_distance(shape.plane, p) < dist_thres) && (std::abs(shape.plane.orthogonal_vector()*n) > normal_thres);
    if (valid) {
      shape.inliers.push_back(neighbour);
      if (shape.inliers.size() % n_refit ==0) {
        cds.fit_plane(shape.inliers, shape.plane);
      }
    }
    return valid;
  }
};
};