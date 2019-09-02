#include <array>
#include <vector>
#include <map>
#include <unordered_set>
#include <iostream>

class IntervalList {
  // list of all vertices: distance, segment_id
  typedef std::multimap<double, size_t> ValueMap;
  typedef std::vector<ValueMap::iterator> Component;
  typedef std::array<double,2> Interval;
  ValueMap values;
  size_t id = 1; //id of 0 is reserved for clip_segments

  public:

  void insert(Interval i) {
    values.insert({i[0], id});
    values.insert({i[1], id++});
  }

  size_t size() {
    return values.size()/2;
  }

  std::vector<Interval> merge_overlap() {
    std::vector<Interval> merged_intervals;
    if (size()==0)
      return merged_intervals;
    else if (size()==1)
      return {{values.begin()->first, values.rbegin()->first}};
    // find connected components
    std::unordered_set<size_t> current;
    std::vector< Component > components;

    // iterate over all values (in sorted order)
    for (auto it=values.begin(); it!=values.end(); ++it) {
      // keep track of current segments
      if (current.find(it->second) == current.end()) {
        current.emplace(it->second);
      } else {
        current.erase(it->second);
      }

      // make a new component if current is empty, store it in current component otherwise
      if (current.empty()) {
        components.resize(components.size()+1);
      } else {
        auto x = it->first;
        std::cout << it-> first << "\n";
        // CRASH!! push_back is not working..
        components.back().push_back(it);
      }
    }
    // remove al values in each component and replace by just the first and last vertex in the component
    
    for (auto& comp : components) {
      double new_source = (*comp.begin())->first;
      // size_t new_id = (*comp.begin())->second;
      double new_target = (*comp.rbegin())->first;
      for (auto it : comp) {
        values.erase(it);
      }
      insert({new_source, new_target});
      merged_intervals.push_back({new_source, new_target});
    }
    
    return merged_intervals;
  }

  // void clip(std::vector<std::array<double,2>> clip_segments) {
  //   // ensures that there are no intervals that overlap with source-target
  //   // we assume the clip segments do not overlap (could be guaranteed by first doing the merge_connected_component() function above)
  //   std::unordered_set<size_t> current;

  //   for (auto& clipseg : clip_segments) {
  //     values.insert({clipseg[0], 0});
  //     values.insert({clipseg[1], 0});
  //   }

  //   // iterate over all values (in sorted order)
  //   bool inside_clipseg = false;
  //   for (auto it=values.begin(); it!=values.end(); ++it) {
  //     if (it->second == 0)
  //       inside_clipseg = !inside_clipseg;

  //     if (current.find(it->second) == current.end()) {
  //       current.emplace(it->second);
  //     } else {
  //       current.erase(it->second);
  //     }
  //   }
  // }
};