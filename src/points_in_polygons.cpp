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

// #include <sstream>
#ifdef GFP_WITH_PDAL
  #include <pdal/PointTable.hpp>
  #include <pdal/PointView.hpp>
  #include <pdal/Options.hpp>
  #include <pdal/io/EptReader.hpp>
  #include <pdal/filters/RangeFilter.hpp>
  #include <pdal/filters/CropFilter.hpp>
#endif

#include <lasreader.hpp>
#include "pip_util.hpp"

namespace geoflow::nodes::stepedge {

class PointsInPolygonsCollector  {
  gfSingleFeatureInputTerminal& polygons;
  gfSingleFeatureOutputTerminal& point_clouds;
  gfSingleFeatureOutputTerminal& ground_elevations;

  // ground elevations
  std::vector<std::vector<arr3f>> ground_buffer_points;
  RasterTools::Raster pindex;
  std::vector<std::vector<size_t>> pindex_vals;
  std::vector<pGridSet> poly_grids, buf_poly_grids;
  std::vector<vec1f> z_ground;
  std::unordered_map<std::unique_ptr<arr3f>, std::vector<size_t>> points_overlap; // point, [poly id's], these are points that intersect with multiple polygons

  int ground_class, building_class;
  
  public:
  Box completearea_bb;
  float min_ground_elevation = std::numeric_limits<float>::max();

  PointsInPolygonsCollector(
    gfSingleFeatureInputTerminal& polygons, 
    gfSingleFeatureInputTerminal& buf_polygons, 
    gfSingleFeatureOutputTerminal& point_clouds, 
    gfSingleFeatureOutputTerminal& ground_elevations, 
    float& cellsize, 
    float& buffer,
    int ground_class=2,
    int building_class=6
    )
    : polygons(polygons), point_clouds(point_clouds), ground_elevations(ground_elevations), ground_class(ground_class), building_class(building_class) {
    // point_clouds_ground.resize(polygons.size());
    point_clouds.resize<PointCollection>(polygons.size());
    z_ground.resize(polygons.size());
    ground_buffer_points.resize(polygons.size());

    for (size_t i=0; i< point_clouds.size(); ++i) {
      point_clouds.get<PointCollection&>(i).add_attribute_vec1i("classification");
    }

    // make a vector of BOX2D for the set of input polygons
    // build point in polygon grids
    for (size_t i=0; i<polygons.size(); ++i) {
      auto ring = polygons.get<LinearRing>(i);
      auto buf_ring = buf_polygons.get<LinearRing>(i);
      poly_grids.push_back(build_grid(ring));
      buf_poly_grids.push_back(build_grid(buf_ring));
      completearea_bb.add(buf_ring.box());
    }

    // build an index grid for the polygons

    // build raster index (pindex) and store for each raster cell all the polygons with intersecting bbox (pindex_vals)
    float minx = completearea_bb.min()[0]-buffer; 
    float miny = completearea_bb.min()[1]-buffer;
    float maxx = completearea_bb.max()[0]+buffer;
    float maxy = completearea_bb.max()[1]+buffer;
    pindex =  RasterTools::Raster(cellsize, minx, maxx, miny, maxy);
    pindex_vals.resize(pindex.dimx_*pindex.dimy_);

    // populate pindex_vals
    for (size_t i=0; i<buf_polygons.size(); ++i) {
      auto ring = buf_polygons.get<LinearRing>(i);
      auto& b = ring.box();
      size_t r_min = pindex.getRow(b.min()[0], b.min()[1]);
      size_t c_min = pindex.getCol(b.min()[0], b.min()[1]);
      size_t r_max = static_cast<size_t>( ceil(b.size_y() / cellsize) ) + r_min;
      size_t c_max = static_cast<size_t>( ceil(b.size_x() / cellsize) ) + c_min;
      if (r_max >= pindex.dimy_)
        r_max = pindex.dimy_-1;
      if (c_max >= pindex.dimx_)
        c_max = pindex.dimx_-1;
      for(size_t r = r_min; r <= r_max; ++r) {
        for(size_t c = c_min; c <= c_max; ++c) {
          pindex_vals[r*pindex.dimx_ + c].push_back(i);
        }
      }
    }

  }

  ~PointsInPolygonsCollector() {
    for (int i=0; i<poly_grids.size(); i++) {
      delete poly_grids[i];
      delete buf_poly_grids[i];
    }
  }

  void add_point(arr3f point, int point_class) {
    // look up grid index cell and do pip for all polygons retreived from that cell
    size_t lincoord = pindex.getLinearCoord(point[0],point[1]);
    if (lincoord >= pindex_vals.size() || lincoord < 0) {
      // std::cout << "Point (" << point[0] << ", " <<point[1] << ", "  << point[2] << ") is not in the polygon bbox.\n";
      return;
    }
    // For the ground points we only test if the point is within the grid
    // index cell, but we do not do a pip for the footprint itself.
    // The reasons for this:
    //  - If we would do a pip for with the ground points, we would need to
    //    have a separate list of buffered footprints that we use for the
    //    ground pip. Still it would be often the case that there are no
    //    ground points found for a buffered footprint. Thus we need a larger
    //    area in which we can guarantee that we find at least a couple of
    //    ground points.
    //  - Because we work with Netherlands data, the ground relief is small.
    //    Thus a single ground height value per grid cell is good enough for
    //    representing the ground/floor elevation of the buildings in that
    //    grid cell.
    pPipoint pipoint = new Pipoint{point[0],point[1]};
    std::vector<size_t> poly_intersect;
    for(size_t& poly_i : pindex_vals[lincoord]) {
      if (GridTest(buf_poly_grids[poly_i], pipoint)) {
        auto& point_cloud = point_clouds.get<PointCollection&>(poly_i);
        auto classification = point_cloud.get_attribute_vec1i("classification");
        
        if (point_class == ground_class) {
          min_ground_elevation = std::min(min_ground_elevation, point[2]);
          z_ground[poly_i].push_back(point[2]);
        }

        if (GridTest(poly_grids[poly_i], pipoint)) {
          if (point_class == ground_class) {            
            point_cloud.push_back(point);
            (*classification).push_back(2);
          } else if (point_class == building_class) {
              poly_intersect.push_back(poly_i);
          }
        } else if (point_class == ground_class) {
          ground_buffer_points[poly_i].push_back( point );
        }
      }
    }

    if (point_class == building_class) {
      if (poly_intersect.size() == 1) {
        auto& point_cloud = point_clouds.get<PointCollection&>(poly_intersect[0]);
        auto classification = point_cloud.get_attribute_vec1i("classification");
        point_cloud.push_back(point);
        (*classification).push_back(6);
      } else if (poly_intersect.size() > 1) {
        points_overlap[std::make_unique<arr3f>(point)] = poly_intersect;
      }
    }
    delete pipoint;
  }

  void do_post_process(
      float& ground_percentile, 
      float& max_density_delta, 
      float& coverage_threshold,
      bool& clear_if_insufficient,
      gfSingleFeatureOutputTerminal& poly_areas,
      gfSingleFeatureOutputTerminal& poly_pt_counts_bld,
      gfSingleFeatureOutputTerminal& poly_pt_counts_grd,
      gfSingleFeatureOutputTerminal& poly_ptcoverage_class,
      gfSingleFeatureOutputTerminal& poly_densities
    ) {
    
    // compute poly properties
    struct PolyInfo { size_t pt_count_bld; size_t pt_count_grd; size_t pt_count_bld_overlap{0}; float avg_elevation; float area; };
    std::unordered_map<size_t, PolyInfo> poly_info;
    
    for (size_t poly_i=0; poly_i < polygons.size(); poly_i++) {
      auto& polygon = polygons.get<LinearRing>(poly_i);
      auto& point_cloud = point_clouds.get<PointCollection&>(poly_i);
      auto classification = point_cloud.get_attribute_vec1i("classification");
      PolyInfo info;
      
      info.area = polygon.signed_area();
      size_t pt_cnt_bld = 0;
      size_t pt_cnt_grd = 0;
      float z_sum = 0;
      for(size_t pi=0; pi < point_cloud.size(); ++pi) {
        if((*classification)[pi] == 6) {
          ++pt_cnt_bld;
          z_sum += point_cloud[pi][2];
        } else if ((*classification)[pi] == 2) {
          ++pt_cnt_grd;
        }
      }
      info.pt_count_bld = pt_cnt_bld;
      info.pt_count_grd = pt_cnt_grd;
      if (info.pt_count_bld > 0) {
        info.avg_elevation = z_sum/pt_cnt_bld;
      }
      poly_info.insert({poly_i, info});
      
    }

    // merge buffer ground points into regular point_clouds now that the proper counts have been established
    for (size_t poly_i; poly_i < polygons.size(); poly_i++) {
      auto& point_cloud = point_clouds.get<PointCollection&>(poly_i);
      auto classification = point_cloud.get_attribute_vec1i("classification");
      for (auto& p : ground_buffer_points[poly_i]) {
        point_cloud.push_back(p);
        (*classification).push_back(2);
      }
    }
    ground_buffer_points.clear();

    // assign points_overlap
    for(auto& [p, polylist] : points_overlap) {
      for( auto& poly_i : polylist ) {
        poly_info[poly_i].pt_count_bld_overlap++;
      }
    }
    for(auto& [p, polylist] : points_overlap) {
      // find best polygon to assign this point to
      std::sort(polylist.begin(), polylist.end(), [&max_density_delta, &poly_info, this](auto& d1, auto& d2) {
        // we look at the maximim possible point density (proxy for point coverage) and the average elevation
        // compute poitncloud density for both polygons
        float pd1 = (poly_info[d1].pt_count_bld + poly_info[d1].pt_count_bld_overlap) / poly_info[d1].area;
        float pd2 = (poly_info[d2].pt_count_bld + poly_info[d2].pt_count_bld_overlap) / poly_info[d2].area;

        // check if the difference in point densities is less than 5%
        if (std::abs(1 - pd1/pd2) < max_density_delta) {
          // if true, then look at the polygon with the highest elevation point cloud
          return poly_info[d1].avg_elevation < poly_info[d2].avg_elevation;
        } else {
          // otherwise decide based on the density values
          return pd1 < pd2;
        }
      });
      
      // now the most suitable polygon (footprint) is the last in the list. We will assign this point to that footprint.
      auto& point_cloud = point_clouds.get<PointCollection&>( polylist.back() );
      auto classification = point_cloud.get_attribute_vec1i("classification");
      point_cloud.push_back(*p);
      (*classification).push_back(6);
      poly_info[ polylist.back() ].pt_count_bld++;
    }

    // Compute average elevation per polygon
    std::cout <<"Computing the average ground elevation per polygon..." << std::endl;
    for (size_t i=0; i<z_ground.size(); ++i) {
      float ground_ele = min_ground_elevation;
      if (z_ground[i].size()!=0) {
        std::sort(z_ground[i].begin(), z_ground[i].end(), [](auto& z1, auto& z2) {
          return z1 < z2;
        });
        int elevation_id = std::floor(ground_percentile*float(z_ground[i].size()-1));
        ground_ele = z_ground[i][elevation_id];
      } else {
        std::cout << "no ground pts found for polygon\n";
      }
      // Assign the median ground elevation to each polygon
      ground_elevations.push_back(ground_ele);
    }

    // clear footprints with very low coverage (ie. underground footprints)
    // TODO: improve method for computing mean_density
    float total_cnt=0, total_area=0;
    for( auto& [poly_i, info] : poly_info ) {
      total_cnt += info.pt_count_bld + info.pt_count_grd;
      total_area += info.area;
    }
    float mean_density = total_cnt/total_area;
    float diff_sum = 0;
    for( auto& [poly_i, info] : poly_info ) {
      diff_sum += std::pow(mean_density - (info.pt_count_bld / info.area), 2);
    }
    float std_dev_density = std::sqrt(diff_sum / poly_info.size());
    std::cout << "Mean point density = " << mean_density << std::endl;
    std::cout << "\t standard deviation = " << std_dev_density << std::endl;

    float cov_thres = mean_density - coverage_threshold * std_dev_density;
    for (size_t poly_i=0; poly_i < polygons.size(); ++poly_i) {
      auto& info = poly_info[poly_i];

      if ( (info.pt_count_bld / info.area) < cov_thres ) {
        if (clear_if_insufficient) {
          auto& point_cloud = point_clouds.get<PointCollection&>( poly_i );
          point_cloud.clear();
        }
        poly_ptcoverage_class.push_back(std::string("insufficient"));
      } else {
        poly_ptcoverage_class.push_back(std::string("sufficient"));
      }
      // info.pt_count = point_cloud.size();
      poly_areas.push_back(float(info.area));
      poly_pt_counts_bld.push_back(int(info.pt_count_bld));
      poly_pt_counts_grd.push_back(int(info.pt_count_grd));
      poly_densities.push_back(float(info.pt_count_bld / info.area));
    }
  }
};

void LASInPolygonsNode::process() {
  auto& polygons = vector_input("polygons");
  auto& buf_polygons = vector_input("buf_polygons");

  auto& point_clouds = vector_output("point_clouds");
  auto& ground_elevations = vector_output("ground_elevations");
  auto& poly_areas = vector_output("poly_areas");
  auto& poly_pt_counts_bld = vector_output("poly_pt_counts_bld");
  auto& poly_pt_counts_grd = vector_output("poly_pt_counts_grd");
  auto& poly_ptcoverage_class = vector_output("poly_ptcoverage_class");
  auto& poly_densities = vector_output("poly_densities");

  PointsInPolygonsCollector pip_collector{
    polygons, 
    buf_polygons, 
    point_clouds, 
    ground_elevations, 
    cellsize, 
    buffer, 
    ground_class, 
    building_class
  };

  auto aoi_min = pip_collector.completearea_bb.min();
  auto aoi_max = pip_collector.completearea_bb.max();

  for (auto filepath : split_string(manager.substitute_globals(filepaths), " "))
  {
    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(filepath.c_str());
    LASreader* lasreader = lasreadopener.open();
    
    if (!lasreader){
      std::cout << "cannot read las file: " << filepath << "\n";
      continue;
    }

    Box file_bbox;
    file_bbox.add(arr3f{
      float(lasreader->get_min_x()-(*manager.data_offset)[0]),
      float(lasreader->get_min_y()-(*manager.data_offset)[1]),
      float(lasreader->get_min_z()-(*manager.data_offset)[2])
    });
    file_bbox.add(arr3f{
      float(lasreader->get_max_x()-(*manager.data_offset)[0]),
      float(lasreader->get_max_y()-(*manager.data_offset)[1]),
      float(lasreader->get_max_z()-(*manager.data_offset)[2])
    });

    if(!file_bbox.intersects(pip_collector.completearea_bb)){
      std::cout << "no intersection footprints with las file: " << filepath << "\n";
      continue;
    }

    // tell lasreader our area of interest. It will then use quadtree indexing if available (.lax file created with lasindex)
    lasreader->inside_rectangle(
      aoi_min[0] + (*manager.data_offset)[0], 
      aoi_min[1] + (*manager.data_offset)[1], 
      aoi_max[0] + (*manager.data_offset)[0], 
      aoi_max[1] + (*manager.data_offset)[1]
    );

    while (lasreader->read_point()) {
      pip_collector.add_point(
        {
        float(lasreader->point.get_x()-(*manager.data_offset)[0]), 
        float(lasreader->point.get_y()-(*manager.data_offset)[1]),
        float(lasreader->point.get_z()-(*manager.data_offset)[2]) 
        }, 
        lasreader->point.get_classification()
      );
    }
    lasreader->close();
    delete lasreader;
  }

  pip_collector.do_post_process(
    ground_percentile, 
    max_density_delta, 
    coverage_threshold, 
    clear_if_insufficient,
    poly_areas, 
    poly_pt_counts_bld, 
    poly_pt_counts_grd, 
    poly_ptcoverage_class, 
    poly_densities
  );
}

#ifdef GFP_WITH_PDAL
void EptInPolygonsNode::process()
{
  // Prepare inputs
  auto& polygons = vector_input("polygons");
  auto& buf_polygons = vector_input("buf_polygons");

  // Prepare outputs
  auto& point_clouds = vector_output("point_clouds");
  auto& ground_point_clouds = vector_output("ground_point_clouds");
  auto& ground_elevations = vector_output("ground_elevations");

  PointsInPolygonsCollector pip_collector{polygons, buf_polygons, point_clouds, ground_point_clouds, ground_elevations, cellsize, buffer};

  std::string ept_path = "ept://" + manager.substitute_globals(dirpath);
  try {
    // This scope is for debugging
    pdal::EptReader reader;
    pdal::Options   options;
    options.add("filename", ept_path);
    reader.setOptions(options);
    // .preview() only reads the metadata
    const pdal::QuickInfo qi(reader.preview());
    std::cout << std::endl << "EPT Bounds:\t" << qi.m_bounds << std::endl;
    std::cout << "EPT Point Count:\t" << qi.m_pointCount << std::endl;
    std::cout << "EPT WKT:\t" << qi.m_srs.getWKT() << std::endl;
  } catch (pdal::pdal_error& e) {
    throw(gfException(e.what()));
  }

  // TODO: The reader and filter init can go to a function that returns a
  //  const PointViewSet.
  auto pmin = pip_collector.completearea_bb.min();
  auto pmax = pip_collector.completearea_bb.max();
  pdal::BOX2D poly_bbox(
    pmin[0] - buffer + (*manager.data_offset)[0],
    pmin[1] - buffer + (*manager.data_offset)[1],
    pmax[0] + buffer + (*manager.data_offset)[0],
    pmax[1] + buffer + (*manager.data_offset)[1]
  );
  pdal::EptReader reader;
  {
    pdal::Options options;
    options.add("filename", ept_path);
    options.add("bounds", poly_bbox);
//      options.add("polygon", wkt.str());
    reader.setOptions(options);
  }
  pdal::RangeFilter range;
  {
    pdal::Options options;
    if (filter_limits.length() == 0) {
      std::cout << "Uh oh, PDAL Range filter cannot be empty. GOING TO CRASH!"
                << std::endl;
    }
    options.add("limits", filter_limits);
    range.setOptions(options);
    range.setInput(reader);
  }
  pdal::PointTable ept_table;
  range.prepare(ept_table);
  const pdal::PointViewSet pw_set(range.execute(ept_table));

  std::cout << "...PDAL Pipeline executed\n Continuing with point-in-polygon tests\n";

  for (const pdal::PointViewPtr& view : pw_set) {

    if (view->size() == 0) {
      std::cout << "0 points were found" << std::endl;
    }
    for (pdal::point_count_t p(0); p < view->size(); ++p) {
      float px = view->getFieldAs<float>(pdal::Dimension::Id::X, p) - (*manager.data_offset)[0];
      float py = view->getFieldAs<float>(pdal::Dimension::Id::Y, p) - (*manager.data_offset)[1];
      float pz = view->getFieldAs<float>(pdal::Dimension::Id::Z, p) - (*manager.data_offset)[2];
      int pclass = view->getFieldAs<int>(pdal::Dimension::Id::Classification, p);
      
      pip_collector.add_point({px, py, pz}, pclass);
    }
  }

  pip_collector.do_post_process(ground_percentile, 0.05, 0.5);
}
#endif

}