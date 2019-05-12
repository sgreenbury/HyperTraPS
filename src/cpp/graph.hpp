#ifndef __GRAPH__HPP
#define __GRAPH__HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <cmath>
#include <map>
#include <deque>
#include <limits>

class Graph
{
public:
  // Data structures for working with graph
  std::map <int, std::vector <int>> out_edges_;
  std::map <int, std::vector <int>> in_edges_;
  std::vector <int> roots_;
  std::vector <int> leaves_;
  std::vector <std::pair<int, int>> edge_list_;
  std::map <int, std::vector <int>> node_labels_;
  std::map <int, int> one_counter_;
  std::map <int, std::deque <int>> path_to_leaves_;
  
  // Methods
  void LoadEdgelist_(const std::string file_name);
  void LoadNodelist_(const std::string file_name);
  void FindRoots_();
  void FindLeaves_();
  void FindPathToLeaf_(int leaf);
  void FindPathToLeaves_();
  void PrintGraph_();
  void Clear_();

  // For copy construction and data structure usage
  Graph();
  Graph(const std::string file_name);
  Graph(const Graph& G);
  Graph& operator=(const Graph& G);

  void Init_(const Graph& G);
  void Init_(const std::string file_name);
};

#endif
