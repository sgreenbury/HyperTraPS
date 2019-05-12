#include "updater_fast.hpp"
#include "graph.hpp"
#include <functional>
#include <unistd.h>

void PrintUsageMessage()
{
  std::string s(100, '-');
  s = s + "\n";
  std::cout << "\n";
  std::cout << s;
  std::cout << "\nHyperTraPS Discrete Time - Random Walk Sampling\n\n";
  std::cout << s
	    << "\nUsage options:\n";
  std::cout << "  -f <input data file = \"transitions.txt\">\n";
  std::cout << "  -L <model type = \"forwards.txt\">\n";
  std::cout << "  -w <walk simulation type  = \"zero-one\" =: \"WS1\">\n";
  std::cout << "     Can take: \"zero-one\" =: \"WS1\" or \"match-data\" =: \"WS2\"\n";
  std::cout << "  -b <burn in from MCMC sampler (assumed posterior sampling after) = 50000>\n";
  
  std::cout << "  -e <limit of posterior samples = -1 (unlimited by default)>\n";
  std::cout << "  -i <iid steps between posterior sample draws  = 100>\n";
  std::cout << "  -g <desired sample steps between sample draws = 100>\n";
  std::cout << "  -R <Number of times dataset is sampled for random walks = 10>\n";

  std::cout << "  -R <Generate resamples of dataset based on target counts = \"no\">\n";

  std::cout << "  -S <seed = 0>\n";
  std::cout << "  -v <verbosity = 0>\n";
  std::cout << "  -h or -? <print this message>\n";
  std::cout << s;
}

void WriteTrajectoryToFile(std::ofstream& tf, const Graph& G, Trajectory& t)
{
  int start_ones = t.start_.CountOnes_();
  if(start_ones == 0)
    tf << "root-";
  for(int i=0; i<t.n_visited_-1; i++) {
    int diff = t.visited_[i+1].Difference(t.visited_[i]);
    if(diff == -1)
      continue;
    tf << diff << "-";
  }
  int final_diff = t.end_.Difference(t.visited_[t.n_visited_-1]);
  if(final_diff == -1)
    tf << "end" << std::endl;
  else
    tf << final_diff << "-end" << std::endl;
  return;
}

void WriteTrajectoryToFileMatchData(std::ofstream& tf, Graph& G, std::map <std::pair<int, int>, Trajectory >& trajectories)
{
  for(std::map <int, std::deque<int>>::iterator it=G.path_to_leaves_.begin();
      it != G.path_to_leaves_.end(); it++) {
    int STARTED = 0;
    std::deque <int>& ptl = std::ref(it->second);
    for(unsigned int j=0; j<ptl.size()-1; j++) {
      int start_node = ptl[j];
      int end_node   = ptl[j+1];
      Trajectory& t   = std::ref(trajectories[std::make_pair(start_node, end_node)]);
      if(t.n_visited_ == 0) {
	if(j+1 == ptl.size()-1 and STARTED == 1)
	  tf << "end" << std::endl;
	continue;
      }
      for(int k=0; k<t.n_visited_; k++) {
	int diff;
	if(k == t.n_visited_ - 1)
	  diff = t.end_.Difference(t.visited_[k]);
	else
	  diff = t.visited_[k+1].Difference(t.visited_[k]);
	
	if(t.visited_[k].CountOnes_() == 0) {
	  STARTED = 1;
	  tf << "root-";
	}
	if(diff == -1)
	  continue;
	tf << diff << "-";
      }
      if(j+1 == ptl.size()-1 and STARTED == 1)
	tf << "end" << std::endl;
    }
  }
  return;
}

int LoadTransitions(const std::string file_name, const std::string load_type, std::vector <Transition>& pis, const Transition pi)
{
  // Loads each of the different possible types of matrix
  std::ifstream in_file(file_name);
  std::string line;
  std::string line_label;
  std::string row_string;
  int row = 0;
  int col = 0;
  double d;
  Transition pi_t = pi;
  pis.clear();
  while(std::getline(in_file, line)) {
    std::stringstream linestream(line);
    std::getline(linestream, line_label, '\t');
    for(row = 0; row < pi_t.L_; row++) {
      std::getline(linestream, row_string, '\t');
      std::stringstream rowstream(row_string);
      col = 0;
      while(rowstream >> d) {
	if((col >= pi_t.L_ and load_type != "other_states") or (col > 2 and load_type == "other_states"))
	    throw std::invalid_argument("Too many columns in data file.");
	if(load_type == "forwards")
	  pi_t.F_[row][col] = d;
	if(load_type == "backwards")
	  pi_t.B_[row][col] = d;
	if(load_type == "other_states")
	  pi_t.O_[row][col] = d;
	col++;
      }
    }
    pis.push_back(pi_t);
  }
  in_file.close();
  return 0;
}

int main(int argc, char *argv[])
{
  uint32_t seed          =  0;
  std::string file_label = "transitions.txt";
  std::string file_name  = "forwards.txt";
  std::string match_type = "zero-one";
  int start_test         = 50000;
  int end_test           = -1;
  int gap                = 100;
  int sample_gap         = 100;
  int M                  = 10;
  std::string resamples  = "no";
  int verbose            = 0;

  int opt = 0;
  static const char *optString = "f:L:w:b:e:i:g:R:S:D:vh?";
  opt = getopt(argc, argv, optString);
  while(opt!=-1) {
    switch(opt) {
    case 'f':
      file_label  = optarg;
      break;
    case 'L':
      file_name   = optarg;
      break;
    case 'w':
      match_type  = optarg;
      break;
    case 'b':
      start_test  = atoi(optarg);
      break;
    case 'e':
      end_test    = atoi(optarg);
      break;
    case 'i':
      gap         = atoi(optarg);
      break;
    case 'g':
      sample_gap  = atoi(optarg);
      break;
    case 'R':
      M           = atoi(optarg);
      break;
    case 'S':
      seed = atoi(optarg);
      break;
    case 'D':
      resamples = optarg;
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h':
    case '?':
      PrintUsageMessage();
      exit(1);
      break;
    default:
      break;
    }
    opt = getopt(argc, argv, optString);
  }

  
  RandomInitialise(seed);
  
  // 1. Find L
  int L=0;
  double sigma=0.1;
  Transition X(L, sigma);
  int a = X.FindL_(file_name, "forwards");
  if(a == 1) {
    std::cout << "Load unsuccessful..." << std::endl;
    exit(1);
  }
  L = a;

  Transition pi(L, sigma);
  Updater U(pi);
  State start(std::vector <int> (pi.L_, 0), 0);
  State end(std::vector <int> (pi.L_, 1), 0);

  
  // 2. Load data of start and end points
  std::vector <std::vector <int>> d;
  LoadData(file_label, d);
  int N = (int)((double)d.size()/2);
  std::vector <int> start_counts;
  std::vector <int> end_counts;
  
  std::vector <State> data_start(N, State(std::vector<int> (L,0),0));
  for(int i=0; i<N; i++) {
    int counter = 0;
    if(verbose)
      std::cout << "Data item end: " << (int)((double)i/2) << std::endl;
    for(unsigned int j=0; j<d[2*i].size(); j++) {
      data_start[i].position_[j] = d[2*i][j];
      if(d[2*i][j] == 1)
	counter++;
      if(verbose)
	std::cout << data_start[i].position_[j] << " ";
    }
    if(verbose)
      std::cout << std::endl;
    start_counts.push_back(counter);
  }
  
  std::vector <State> data_end(N, State(std::vector<int> (L,0),0));
  for(int i=0; i<N; i++) {
    int counter = 0;
    if(verbose)
      std::cout << "Data item end: " << (int)((double)i/2) << std::endl;
    for(unsigned int j=0; j<d[2*i+1].size(); j++) {
      data_end[i].position_[j] = d[2*i+1][j];
      if(d[2*i+1][j] == 1)
	counter++;
      if(verbose)
	std::cout << data_end[i].position_[j] << " ";
    }
    if(verbose)
      std::cout << std::endl;
    end_counts.push_back(counter);
  }

  // Make zero and one states:
  State zero_state(std::vector<int> (L,0),0);
  State one_state(std::vector<int> (L,1),0);
  
  // 3. Set the sample framework up
  int start_index = (int)((double)start_test/gap);

  // Load all the samples
  std::vector <Transition> pis;
  LoadTransitions(file_name, "forwards", pis, pi);

  int end_index = (end_test == -1 ? pis.size() : (int)((double)end_test/gap));
  
  std::string tag = "forwards_";
  std::size_t pos1 = file_name.find(tag);
  std::size_t pos2 = pos1 + tag.size();
  
  std::string out_file = file_name.substr(0, file_name.size()-4);
  if(match_type != "zero-one" and 0)
    out_file += "_list-pord.csv";
  else
    out_file += "_list-pord-" + match_type + ".csv";

  // Load graph
  Graph G;
  if(match_type == "match-data")
    G.Init_(file_label);
  std::map <std::pair<long long, long long>, int> transition_counter;
  std::map <long long, int> state_counter;
  
  std::cout << out_file << std::endl;
  std::ofstream step_list(out_file);
  step_list << "step,feature\n";

  std::vector <std::vector <double>> ordering(L, std::vector <double> (L, 0.0));

  // Trajectory file list
  std::ofstream traj_file(file_name.substr(0, file_name.size()-4) + "_list-pord-trajectory-" + match_type + ".txt");

  // Make match-data trajectory map
  std::map <std::pair<int, int>, Trajectory> trajectories;

  // TODO: -w "match-summary" option allows walks that result in the same number
  //       of acquisitions between source and target and subsequent regeneration of new
  //       "matched summary" sample data. However, this is currently non-functional.
  
  // Begin the looping 
  for(int m=0; m<M; m++) {
    // Determine whether to write to outfile
    std::ofstream rs;
    if(resamples == "yes" and match_type == "match-summary") {
      std::string resample_file = file_name.substr(0,pos1-1) + "_resamples-m_" + std::to_string(m) + "_" + file_name.substr(pos2, file_name.size());
      std::cout << resample_file << std::endl;
      rs.open(resample_file);
    }

    int data_limit = (match_type != "match-data" ? N : G.edge_list_.size());
    for(int n=0; n<data_limit; n++) {
      // Draw a random number for the sample to be used for the walk
      int k = end_index - start_index;
      double spacer = (double)sample_gap/gap;
      k = (int)(k / spacer);
      k = UniformRealRange(0,1)* k;
      k = k * ((double)sample_gap/gap) + start_index;
      pi = pis[k];
      
      State start_p;
      State end_p;

      // If match type is match_summary, walk to data_start[n]
      if(match_type == "match-summary") {
	Trajectory t(start, end);
	start_p = U.PerformWalk(t, pi, ordering, start, end, step_list,
				0, start_counts[n], "no");
	Trajectory t1(start_p, end);
	end_p   = U.PerformWalk(t1, pi, ordering, start_p, end, step_list,
				start_counts[n], end_counts[n], "yes");


	// Record new start and end data items as resampled data points
	if(resamples == "yes" and match_type == "match-summary") {
	  for(unsigned int j=0; j<start_p.position_.size(); j++) {
	    if(j < start_p.position_.size()-1)
	      rs << start_p.position_[j] << " ";
	    else
	      rs << start_p.position_[j] << std::endl;
	  }
	  for(unsigned int j=0; j<end_p.position_.size(); j++) {
	    if(j < end_p.position_.size()-1)
	      rs << end_p.position_[j] << " ";
	    else
	      rs << end_p.position_[j] << std::endl;
	  }
	}
      }

      if(match_type == "match-data") {
	start_p = State(G.node_labels_[G.edge_list_[n].first], 0);
	end_p   = State(G.node_labels_[G.edge_list_[n].second], 0);
	Trajectory t(start_p, end_p);
	end_p   = U.PerformWalk(t, pi, ordering, start_p, end_p, step_list,
				G.one_counter_[G.edge_list_[n].first],
				G.one_counter_[G.edge_list_[n].second],
				transition_counter,
				state_counter,
				"yes");
	trajectories[std::make_pair(G.edge_list_[n].first, G.edge_list_[n].second)] = t;
      }
      
      if(match_type == "zero-one") {
	start_p = zero_state;
	end_p   = one_state;
	Trajectory t(start_p, end_p);
	end_p   = U.PerformWalk(t, pi, ordering, start_p, end_p, step_list,
				0,
				L,
				transition_counter,
				state_counter,
				"yes");
	WriteTrajectoryToFile(traj_file, G, t);
      }	
    }
    if(match_type == "match-data") {
      WriteTrajectoryToFileMatchData(traj_file, G, trajectories);
    }

    // Close new dataset if resampling
    if(resamples == "yes" and match_type == "match-summary") {
      rs.close();
    }
  }
  step_list.close();
  traj_file.close();

  std::ofstream tcf(file_name.substr(0, file_name.size()-4) + "_list-pord-edge-list-long-" + match_type + ".txt");
  tcf << "from to weight\n";
  for(std::map<std::pair<long long, long long>, int>::iterator it=transition_counter.begin();
      it!=transition_counter.end();
      it++) {
    tcf << it->first.first << " "
	<< it->first.second << " "
	<< it->second << std::endl;
  }
  tcf.close();

  std::ofstream scf(file_name.substr(0, file_name.size()-4) + "_list-pord-state-list-long-" + match_type + ".txt");
  scf << "state count\n";
  for(std::map<long long, int>::iterator it=state_counter.begin();
      it!=state_counter.end();
      it++) {
    scf << it->first << " "
	<< it->second << std::endl;
  }
  scf.close();

  return 0;
}
