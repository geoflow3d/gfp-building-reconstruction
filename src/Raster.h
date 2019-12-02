// Copyright (c) 2012, 2013
// Ravi Peters -- r.y.peters@tudelft.nl
// All rights reserved
// 
// This file is part of Surfonoi.
// 
// Surfonoi is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Surfonoi is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Surfonoi.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <iostream>
#include <cmath>
#include <algorithm>
#include <cfloat>
#include <fstream>
#include <array>
#include <vector>

// #include <gdal_priv.h>
// #include <cpl_string.h>
// #include <cpl_conv.h>
// #include <ogr_spatialref.h>
namespace RasterTools {
    enum alg {AVG,MIN,MAX,CNT};
    class Raster
    {
    public:
        typedef std::array<float,3> point3d;
        typedef std::array<float,2> point2d;
        Raster(double cellsize, double min_x, double max_x, double min_y, double max_y);
        Raster(){};
        ~Raster(){};
        void prefill_arrays(alg a);
        void add_point(double x, double y, double z, alg a);
        void add_raster(double x, double y, double z, alg a);
        size_t getLinearCoord(double x, double y) const;
        std::array<double,2> getColRowCoord(double x, double y) const;
        point3d getPointFromRasterCoords(size_t col, size_t row) const;
        double getNoDataVal() {return noDataVal_;};
        double sample(double &x, double &y);
        void set_val(size_t col, size_t row, double val);
        // void write(const char* WKGCS, alg a, void * dataPtr, const char* outFile);

        // rasterise a polygon and return a list with points - one in the center of each pixel inside the polygon
        // in the polygon first point is *not* repeated as last
        // T should be a vector of arr<float,2> or arr<float,3>
        template<typename T> std::vector<point3d> rasterise_polygon(T& polygon, bool returnNoData=true) const {
          // code adapted from http://alienryderflex.com/polygon_fill/
          int n_nodes, pixelX, pixelY, i, j, swap ;
          int n_vertices = polygon.size();
          std::vector<point3d> result;

          // perhaps we can specialise these to the bounding box of the polygon
          int IMAGE_TOP = 0, IMAGE_BOT = dimy_, IMAGE_LEFT = 0, IMAGE_RIGHT=dimx_;

          // Loop through the rows of the image.
          for (pixelY=IMAGE_TOP; pixelY<IMAGE_BOT; pixelY++) {
            std::vector<int> intersect_x; // vector to hold the x-coordinates where the scanline intersects the polygon

            // Build a list of nodes.
            n_nodes=0; j=n_vertices-1;
            for (i=0; i<n_vertices; i++) {
              auto pi = getColRowCoord((double)polygon[i][0], (double)polygon[i][1]);
              auto pj = getColRowCoord((double)polygon[j][0], (double)polygon[j][1]);
              // std::cerr << pi[0] << " " << pi[1] << "\n";
              // std::cerr << pj[0] << " " << pj[1] << "\n";
              if ( (pi[1]<(double) pixelY && pj[1]>=(double) pixelY)
              || (pj[1]<(double) pixelY && pi[1]>=(double) pixelY)) {
                intersect_x.push_back((int) (pi[0]+(pixelY-pi[1])/(pj[1]-pi[1])
                *(pj[0]-pi[0])));
                ++n_nodes;
              }
              j=i; 
            }

            // Sort the nodes, via a simple “Bubble” sort.
            i=0;
            while (i<n_nodes-1) {
              if (intersect_x[i]>intersect_x[i+1]) {
                swap=intersect_x[i]; intersect_x[i]=intersect_x[i+1]; intersect_x[i+1]=swap; if (i) i--; 
              } else {
                i++; 
              }
            }

            // Fill the pixels between node pairs.
            for (i=0; i<n_nodes; i+=2) {
              if  (intersect_x[i ]>=IMAGE_RIGHT) 
                break;
              if  (intersect_x[i+1]> IMAGE_LEFT ) {
                if (intersect_x[i ]< IMAGE_LEFT ) 
                  intersect_x[i ]=IMAGE_LEFT ;
                if (intersect_x[i+1]> IMAGE_RIGHT) 
                  intersect_x[i+1]=IMAGE_RIGHT;
                for (pixelX=intersect_x[i]; pixelX<intersect_x[i+1]; pixelX++) 
                  if (returnNoData) {
                    result.push_back(getPointFromRasterCoords(pixelX,pixelY)); 
                  } else if((vals_[pixelX+pixelY*dimx_] != noDataVal_) ) {
                    result.push_back(getPointFromRasterCoords(pixelX,pixelY)); 
                  }
              }
            }
          }

          return result;
        };
        // template<> std::vector<point3d> rasterise_polygon(std::vector<point2d>& polygon) const;
        // template<> std::vector<point3d> rasterise_polygon(std::vector<point3d>& polygon) const;

        double cellSize_, minx_, miny_, maxx_, maxy_;
        int dimx_, dimy_;
        double noDataVal_;
        std::vector<int16_t> counts_;
        std::vector<double> vals_;
    private:
        void avg(double &x, double &y, double &val);
        void min(double &x, double &y, double &val);
        void max(double &x, double &y, double &val);
        void cnt(double &x, double &y);
        // OGRSpatialReference oSRS;
    };
}