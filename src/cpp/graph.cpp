#include "graph.hpp"

Graph::Graph()
{
}

Graph::Graph(const std::string file_name)
{
  Init_(file_name);
}

void Graph::Init_(const std::string file_name)
{
  std::string base = file_name.substr(0, file_name.size()-4);
  LoadNodelist_(base + ".nodelist");
  LoadEdgelist_(base + ".edgelist");
  FindRoots_();
  FindLeaves_();
  FindPathToLeaves_();
}

void Graph::Clear_()
{
  node_labels_.clear();
  in_edges_.clear();
  out_edges_.clear();
  roots_.clear();
  leaves_.clear();

  edge_list_.clear();
  one_counter_.clear();

  path_to_leaves_.clear();
}

void Graph::FindPathToLeaf_(int leaf)
{
  std::deque <int> path;
  int current_node = leaf;
  path.push_front(current_node);
  while(in_edges_[current_node].size() != 0) {
    current_node = in_edges_[current_node][0];
    path.push_front(current_node);
  }
  path_to_leaves_[leaf] = path;
}

void Graph::FindPathToLeaves_()
{
  for(unsigned int i=0; i<leaves_.size(); i++)
    FindPathToLeaf_(leaves_[i]);
}

void Graph::Init_(const Graph& G)
{
  node_labels_.insert(G.node_labels_.begin(), G.node_labels_.end());
  in_edges_.insert(G.in_edges_.begin(), G.in_edges_.end());
  out_edges_.insert(G.out_edges_.begin(), G.out_edges_.end());
  one_counter_.insert(G.one_counter_.begin(), G.one_counter_.end());

  for(unsigned int j=0; j<G.edge_list_.size(); j++)
    edge_list_.push_back(G.edge_list_[j]);

  FindRoots_();
  FindLeaves_();
  FindPathToLeaves_();
}

Graph& Graph::operator=(const Graph &G)
{
  if( &G != this ) {
    Clear_();
    Init_(G);
  }
  return *this;
}


void Graph::LoadEdgelist_(const std::string file_name)
{
  for(std::map<int, std::vector<int>>::iterator it=node_labels_.begin(); it!=node_labels_.end(); it++) {
    out_edges_[it->first] = std::vector <int> ();
    in_edges_[it->first]  = std::vector <int> ();
  }
  std::string line;
  std::string start_string;
  std::string end_string;
  std::string line_label;
  std::ifstream edge_file(file_name);
  if(!edge_file.is_open()){
    std::cout << "Graph file is not open..." << std::endl;
    exit(1);
  }
  edge_list_.clear();
  while(std::getline(edge_file, line) ) {
    std::vector <int> temp;
    std::stringstream linestream(line);
    std::getline(linestream, start_string, '\t');
    std::getline(linestream, end_string, '\t');
    int source = std::stoi(start_string);
    int target = std::stoi(end_string);
    out_edges_[source].push_back(target);
    in_edges_[target].push_back(source);
    edge_list_.push_back(std::make_pair(source, target));
  }
  edge_file.close();
}

void Graph::LoadNodelist_(const std::string file_name)
{
  std::string line;
  std::string line_label;
  std::string label_string;
  std::string temporary_string;
  std::ifstream node_file(file_name);
  int d=0;
  while(std::getline(node_file, line)) {
    std::vector <int> temp;
    std::stringstream linestream(line);
    std::getline(linestream, line_label, '\t');
    std::getline(linestream, label_string, '\t');
    std::stringstream labelstream(label_string);
    while(labelstream >> d) {
      temp.push_back(d);
    }
    node_labels_[std::stoi(line_label)] = temp;
    int counter = 0;
    for(unsigned int i=0; i<temp.size(); i++)
      counter += (temp[i] == 1 ? 1 : 0);
    one_counter_[std::stoi(line_label)] = counter;
  }
  node_file.close();
}

void Graph::FindRoots_()
{
  for(std::map<int, std::vector<int>>::iterator it = in_edges_.begin();
      it != in_edges_.end();
      it++) {
    if(it->second.size() == 0)
      roots_.push_back(it->first);
  }
}

void Graph::FindLeaves_()
{
  for(std::map<int, std::vector<int>>::iterator it = out_edges_.begin();
      it != out_edges_.end();
      it++) {
    if(it->second.size() == 0)
      leaves_.push_back(it->first);
  }
}

void Graph::PrintGraph_()
{
  for(unsigned int i=0; i<edge_list_.size(); i++) {
    std::cout << edge_list_[i].first << "\t"
	      << edge_list_[i].second << std::endl;
  }
  
  for(std::map<int, std::vector<int>>::iterator it=node_labels_.begin(); it != node_labels_.end(); it++) {
    std::cout << "Node: " << it->first << std::endl;
    for(unsigned int i=0; i<it->second.size(); i++)
      std::cout << it->second[i] << " ";
    std::cout << "\t" << one_counter_[it->first] << std::endl;

    std::cout << "Out edges: " << std::endl;
    for(unsigned int i=0; i<out_edges_[it->first].size(); i++) {
      std::cout << out_edges_[it->first][i] << " ";
    }
    std::cout << std::endl;
    std::cout << "In edges: " << std::endl;
    for(unsigned int i=0; i<in_edges_[it->first].size(); i++) {
      std::cout << in_edges_[it->first][i] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "Roots: " << roots_.size() << std::endl;
  for(unsigned int i=0; i<roots_.size(); i++) {
    std::cout << roots_[i] << std::endl;
  }

  std::cout << "Leaves: " << leaves_.size() << std::endl;
  for(unsigned int i=0; i<leaves_.size(); i++) {
    std::cout << leaves_[i] << std::endl;
  }

  std::cout << "Paths to leaves: " << std::endl;
  for(std::map <int, std::deque <int>>::iterator it=path_to_leaves_.begin();
      it != path_to_leaves_.end();
      it++) {
    std::cout << it->first << "\t";
    for(unsigned int j=0; j<it->second.size(); j++) {
      std::cout << it->second[j];
      if(j < it->second.size()-1)
	std::cout << "-";
    }
    std::cout << std::endl;
  }
}
