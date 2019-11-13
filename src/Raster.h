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
        Raster(double cellsize, double min_x, double max_x, double min_y, double max_y);
        ~Raster();
        void prefill_arrays(alg a);
        void add_point(double x, double y, double z, alg a);
        size_t getLinearCoord(double &x, double &y);
        inline std::array<double,2> getColRowCoord(double x, double y);
        point3d getPointFromRasterCoords(size_t col, size_t row);
        double sample(double &x, double &y);
        // void write(const char* WKGCS, alg a, void * dataPtr, const char* outFile);

        // rasterise a polygon and return a list with points - one in the center of each pixel inside the polygon
        // in the polygon first point is *not* repeated as last
        std::vector<point3d> rasterise_polygon(std::vector<point3d>& polygon);

        double cellSize_, minx_, miny_, maxx_, maxy_;
        int dimx_, dimy_;
    private:
        void avg(double &x, double &y, double &val);
        void min(double &x, double &y, double &val);
        void max(double &x, double &y, double &val);
        void cnt(double &x, double &y);
        // OGRSpatialReference oSRS;
        double noDataVal_;
        int16_t *counts_;
        double *vals_;
    };
}