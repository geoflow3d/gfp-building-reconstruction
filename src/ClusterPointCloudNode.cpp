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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/cluster_point_set.h>
#include <CGAL/compute_average_spacing.h>

#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

#include <fstream>
#include <iostream>

namespace geoflow::nodes::stepedge {

  void ClusterPointCloudNode::process() {
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point_3 = Kernel::Point_3;
    using Point_set = CGAL::Point_set_3<Point_3>;


    auto& points = input("points").get<PointCollection&>();

    Point_set point_set;

    for(auto& p : points) {
      if (flatten) {
        point_set.insert(Point_3(p[0], p[1], 0));
      } else { 
        point_set.insert(Point_3(p[0], p[1], p[2]));
      }
    }

    // Add a cluster map
    Point_set::Property_map<int> cluster_map = point_set.add_property_map<int>("cluster", -1).first;

    // Compute average spacing
    // double spacing = CGAL::compute_average_spacing<CGAL::Parallel_if_available_tag> (point_set, 12);
    // std::cerr << "Spacing = " << spacing << std::endl;

    // Adjacencies stored in vector
    std::vector<std::pair<std::size_t, std::size_t> > adjacencies;
    // Compute clusters
    CGAL::Real_timer t;
    t.start();
    // unsigned int k=15;
    std::size_t nb_clusters
      = CGAL::cluster_point_set(point_set, cluster_map,
                                point_set.parameters().neighbor_radius(spacing)
                                                  .adjacencies(std::back_inserter(adjacencies)));
    t.stop();
    std::cerr << "Found " << nb_clusters << " clusters with " << adjacencies.size()
              << " adjacencies in " << t.time() << " seconds" << std::endl;

    // Output a colored PLY file
    Point_set::Property_map<unsigned char> red = point_set.add_property_map<unsigned char>("red", 0).first;
    Point_set::Property_map<unsigned char> green = point_set.add_property_map<unsigned char>("green", 0).first;
    Point_set::Property_map<unsigned char> blue = point_set.add_property_map<unsigned char>("blue", 0).first;
    vec1i cluster_idx(points.size());
    size_t i=0;
    IndexedPlanesWithPoints pts_per_roofplane;
    for(Point_set::const_iterator it = point_set.begin();
       it != point_set.end(); ++ it)
    {
      pts_per_roofplane[cluster_map[i]].second.push_back(
        // point_set.point(*it)
        Point_3(points[i][0], points[i][1], points[i][2])
        );
      // One color per cluster
      cluster_idx[i] = cluster_map[i];
      CGAL::Random rand (cluster_map[i]);
      red[i] = rand.get_int(64, 192);
      green[i] = rand.get_int(64, 192);
      blue[i] = rand.get_int(64, 192);
      ++i;
    }
    std::ofstream ofile("out.ply", std::ios_base::binary);
    CGAL::IO::set_binary_mode(ofile);
    ofile << point_set;

    output("cluster_id").set(cluster_idx);
    output("pts_per_roofplane").set(pts_per_roofplane);

  }
}