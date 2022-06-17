# gfp-building-reconstruction
Geoflow plugin for building LoD2 reconstruction from a point cloud

## Installation

### As part of geoflow-bundle
This the recommended way since it will also include the commonly used plugins and flowcharts to get you started quickly. Also binary pacakges are available.

### Building from source
Need to install first [geoflow](https://github.com/geoflow3d/geoflow).
Requires compiler with c++17 support (see https://en.cppreference.com/w/cpp/compiler_support).

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --build . --parallel 4 --config Release
cmake --build . --target install
```

Dependencies:

* [LASlib](https://github.com/LAStools/LAStools)
* [CGAL](https://github.com/CGAL/cgal)
