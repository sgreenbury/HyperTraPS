#include "likelihood_sampler_thread.hpp"
#include <unistd.h>

void PrintUsageMessage()
{
  std::string s(100, '-');
  s = s + "\n";
  std::cout << "\n";
  std::cout << s;
  std::cout << "\nHyperTraPS Discrete Time - regularisation and model selection\n\n";
  std::cout << s;
  std::cout << "\nDescription:\n";
  std::cout << "    Program for the purpose of perfoming model selection\n"
	    << "    on the basis of:\n"
	    << "      - A zero order model has been run\n"
	    << "      - A first order model has been run\n"
	    << "\n"
	    << "    Program looks up the max likelihood (in 100 steps) across MC run.\n"
	    << "\n"
	    << "    Loads matrix in each case.\n"
	    << "\n"
	    << "    Strips away elements one by one calculating the BIC and AIC for each new parameterisation.\n"
	    << "\n"
	    << "    Outputs the:\n"
	    << "    Parameters,DPs,log_lik,AIC,BIC to a file for each model type\n"
	    << "\n"
	    << "    Outputs the matrix that is optimal for each model type.\n"
	    << "\n"
	    << "\nUsage options:\n";
  std::cout << "  -f <input data file = \"transitions.txt\">\n";
  std::cout << "  -L <model type = \"forwards.txt\">\n";
  std::cout << "  -t <number of threads  = hardware_concurrency - 1>\n";
  std::cout << "  -T <number of tests for greedy backward selection = 100>\n";
  std::cout << "  -G <number of random pruning steps before greedy process = 0>\n";
  std::cout << "  -r <number of runs for likelihood calculation = 200>\n";
  std::cout << "  -I <the information criterion measure for optimisation = \"AIC\" (or \"BIC\")>\n";
  std::cout << "  -S <seed = 0>\n";
  std::cout << "  -v <verbosity = 0>\n";
  std::cout << "  -h or -? <print this message>\n";
  std::cout << s;
}

double BIC(const int n, const int k, const double max_ll)
{
  return k*log(n) - 2*max_ll;
}

double AIC(const int k, const double max_ll)
{
  return 2*(k-max_ll);
}

double IC(const int n, const int k, const double max_ll, const std::string ic_type)
{
  if(ic_type == "AIC")
    return AIC(k, max_ll);
  if(ic_type == "BIC")
    return BIC(n, k, max_ll);
  throw std::invalid_argument("Invalid argument for ic_type...");
}

double MaxLL(Likelihood& l, int TESTS,
	     std::vector <Transition>& pis,
	     std::vector <Updater>& us,
	     std::vector <State>& data_end,
	     std::vector <State>& data_start,
	     std::vector <Trajectory>& ts)
{
  std::vector <double> results;
  for(int i=0; i<TESTS; i++) {
    l.RunLikelihoodCalculation(data_end, ts, pis, us, data_start);
    results.push_back(l.total_log_likelihood_);
  }
  return *std::max_element(results.begin(), results.end());
}

void MaxLikelihoodFromMC(const std::string file_name,
			 std::vector <std::pair <double, std::string>>& lls)
{
  // Use the forwards file-name passed to the program, but
  // substitute stats for forwards in the file string.
  std::string base      = file_name;
  std::string str_stats = "stats";
  std::size_t found     = base.find("forwards");
  base = base.replace(found, 8, str_stats);

  base = "stats.csv";
  // Look through the stats file and compile lls and string
  std::ifstream in_file(base);
  lls.clear();
  std::string line;
  std::string el;
  std::vector <std::string> ds;
  std::getline(in_file, line);
  int row = 0;
  while(std::getline(in_file, line)) {
    // Assumed that there are 100 mcmc iterations for each output parameterisation
    if(row % 100 != 0) {
      row++;
      continue;
    }
    row++;
    std::stringstream linestream(line);
    ds.clear();
    while(getline(linestream, el, ',')){
      ds.push_back(el);
    }
    double ll = std::stod(ds[1]);
    std::string lab = ds[0];
    std::pair <double, std::string> p(ll, lab);
    lls.push_back(p);
  }
  in_file.close();

  // Order lls by size of first el
  std::sort(lls.rbegin(), lls.rend());
}

void MakeBestICFromOutputs(const std::string file_name,
			   const std::string file_name_matrix,
			   std::vector <std::pair <double, std::string>>& ics,
			   int L_in,
			   int out_limit,
			   const std::string file_name_original)
{
  // Look through the stats file and compile lls and string
  std::ifstream in_file(file_name);
  ics.clear();
  std::string line;
  std::string el;
  std::vector <std::string> ds;
  std::getline(in_file, line);
  int row = 0;
  std::getline(in_file, line);
  while(std::getline(in_file, line)) {
    std::stringstream linestream(line);
    ds.clear();
    while(getline(linestream, el, ',')){
      ds.push_back(el);
    }
    double ic = std::stod(ds[4]);
    std::string label = ds[6];
    
    std::pair <double, std::string> p(ic, label);
    ics.push_back(p);
    row++;
  }
  in_file.close();

  // Order ics by size of first el
  std::sort(ics.begin(), ics.end());

  // Output top number_of_tests as regularised models
  Transition pi(L_in, 0.01);
  std::ofstream out_file_matrix(file_name_original.substr(0,file_name_original.size()-4) + "_regularised.txt");
  for(unsigned int i=0; i<out_limit; i++) {
    std::cout << ics[i].first << "\t" << ics[i].second << std::endl;
    pi.Load_(file_name_matrix,
	     "forwards",
	     ics[i].second);
    out_file_matrix << ics[i].second << "\t";
    pi.Write_("forwards", out_file_matrix);
  }
  out_file_matrix.close();
}


int main(int argc, char *argv[])
{
  // Defaults
  std::string file_label      = "transitions.txt";
  int THREADS                 = std::thread::hardware_concurrency()-1;
  std::string file_name       = "forwards.txt";
  int number_of_tests         = 100;
  int parameters_greedy_begin = 0;
  int RUNS                    = 200;
  std::string ic_type         = "AIC";
  int verbose                 = 0;
  int seed                    = 0;

  int opt = 0;
  static const char *optString = "f:t:L:T:G:r:I:vS:h?";
  opt = getopt(argc, argv, optString);
  while(opt!=-1) {
    switch(opt) {
    case 'f':
      file_label              = optarg;
      break;
    case 'L':
      file_name               = optarg;
      break;
    case 'T':
      number_of_tests         = atoi(optarg);
      break;
    case 't':
      THREADS                 = atoi(optarg);
      break;
    case 'G':
      parameters_greedy_begin = atoi(optarg);
      break;
    case 'r':
      RUNS                    = atoi(optarg);
      break;
    case 'I':
      ic_type                 = optarg;
      break;
    case 'v':
      verbose                 = 1;
      break;
    case 'S':
      seed                    = atoi(optarg);
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
  
    
  std::vector <std::pair<double, std::string>> lls;
  std::vector <std::vector <int>> d;
  std::string file       = file_label;
  LoadData(file, d);
  int N                  = (int)((double)d.size()/2);
  int L                  = d[0].size();
  
  State start(std::vector <int> (L, 0), 0);
  State end(std::vector <int> (L, 1), 0);
  
  // Set-up data
  std::vector <State> data_start(N, State(std::vector<int> (L,0),0));
  for(int i=0; i<N; i++) {
    if(verbose)
      std::cout << "Data item end: " << (int)((double)i/2) << std::endl;
    for(unsigned int j=0; j<d[2*i].size(); j++) {
      data_start[i].position_[j] = d[2*i][j];
      if(verbose)
	std::cout << data_start[i].position_[j] << " ";
    }
    if(verbose)
      std::cout << std::endl;
  }

  std::vector <State> data_end(N, State(std::vector<int> (L,0),0));
  for(int i=0; i<N; i++) {
    if(verbose)
      std::cout << "Data item end: " << (int)((double)i/2) << std::endl;
    for(unsigned int j=0; j<d[2*i+1].size(); j++) {
      data_end[i].position_[j] = d[2*i+1][j];
      if(verbose)
	std::cout << data_end[i].position_[j] << " ";
    }
    if(verbose)
      std::cout << std::endl;
  }

  /* Main likelihood calculation requirements */
  // 1. Load pi from posterior
  // 2. Calculate max likelihood, AIC/BIC for full FO parameters
  // 3a. Loop through matrix setting the element that causes the smallest change in likelihood to 0
  // 3b. Calculate the IC for the new parameterisation.
  //    3ba. If less than previous IC, accept new matrix and set element to 0.
  //         Otherwise, stop and return the matrix.

  std::vector <Trajectory> ts;
  std::vector <Updater> us;
  std::vector <Transition> pis;
  for(int i=0; i<THREADS; i++) {
    Trajectory t(start, end);
    Transition pi(L, 0.1);
    Updater u(pi);
    us.push_back(u);
    ts.push_back(t);
    pis.push_back(pi);
  }

  int TESTS = 1;
  Likelihood l(N, RUNS, THREADS, seed);

  Transition pi_current_aic(L, 0.1);
 
  MaxLikelihoodFromMC(file_name, lls);

  std::string out_file_name;
  std::string out_file_matrix_name;
  out_file_name        = file_name.substr(0,file_name.size()-4) + "_information-criterion-" + std::to_string(parameters_greedy_begin) + ".csv";
  out_file_matrix_name = file_name.substr(0,file_name.size()-4) + "_information-criterion-matrix-" + std::to_string(parameters_greedy_begin) + ".txt";
  
  std::cout << out_file_name << std::endl;
  std::cout << out_file_matrix_name << std::endl;

  std::ofstream out_file       (out_file_name);
  std::ofstream out_file_matrix(out_file_matrix_name);

  out_file << "parameters"                  << ","
	   << "data_points"                 << ","
	   << "log_likelihood_aic"          << ","
	   << "log_likelihood_bic"          << ","
	   << "aic_minimum"                 << ","
	   << "bic_minimum"                 << ","
	   << "label" << std::endl;
  
  for(int n=0; n<number_of_tests; n++) {
    std::cout << n << std::endl;
    pis[0].Load_(file_name,
		 "forwards",
		 lls[n].second);
    
    for(int j=0; j<THREADS; j++)
      pis[j] = pis[0];

    int parameters = 0;
    for(unsigned int i=0; i<pis[0].F_.size(); i++) {
      for(unsigned int j=0; j<pis[0].F_[i].size(); j++) {
	if(!(DoubleEqual(pis[0].F_[i][j], 0.0)))
	  parameters++;
      }
    }
    
    int parameters_current_aic = parameters;    
    pi_current_aic = pis[0];
    
    double ll_current_aic      = MaxLL(l, TESTS, pis, us, data_end, data_start, ts);
    double ll_test             = 0.0;
    double ll_max_aic          = ll_current_aic;
    double aic_current         = IC(N, parameters_current_aic, ll_current_aic, ic_type);
    double aic_test            = 0.0;
    double aic_min             = aic_current;
    
    double temp_store          = 0.0;
    int row_aic_min            = 0;
    int col_aic_min            = 0;
    std::string label = std::to_string(n) + "." + std::to_string(parameters);
    out_file << parameters                    << ","
	     << N                             << ","
	     << ll_max_aic                    << ","
	     << "NaN"                         << ","
	     << aic_min                       << ","
	     << "NaN"                         << ","
	     << label                         << std::endl;

    
    out_file_matrix << std::to_string(n) + "." + std::to_string(parameters) << "\t";
    pis[0].Write_("forwards", out_file_matrix);
    
    int parameters_start = parameters;
    while(parameters > 0) {
      if(verbose) {
	std::cout << "Parameters testing: "      << parameters  << "\t";
	std::cout << "Min ic: " << aic_current << "\t";
	std::cout << "Min ic ll: " << ll_current_aic << "\t";
	std::cout << "Min ic params: " << parameters_current_aic << std::endl;
	std::cout << "ic type: " << ic_type << std::endl;
      }
      aic_min = std::numeric_limits<double>::max();
      parameters--;
      if(parameters_start - parameters > parameters_greedy_begin) {
	for(unsigned int row=0; row<pis[0].F_.size(); row++) {
	  for(unsigned int col=0; col<pis[0].F_[row].size(); col++) {
	    temp_store = pis[0].F_[row][col];
	    if(DoubleEqual(temp_store, 0.0))
	      continue;
	    pis[0].F_[row][col] = 0.0;
	    ll_test = MaxLL(l, TESTS, pis, us, data_end, data_start, ts);
	    aic_test = IC(N, parameters, ll_test, ic_type);
	    if(aic_test < aic_min) {
	      aic_min     = aic_test;
	      ll_max_aic  = ll_test;
	      row_aic_min = row;
	      col_aic_min = col;
	    }
	    
	    pis[0].F_[row][col] = temp_store;
	  }
	}
      }
      else {
	row_aic_min = (int)(UniformRealRange(0,1)*L);
	col_aic_min = (int)(UniformRealRange(0,1)*L);
	pis[0].F_[row_aic_min][col_aic_min] = 0.0;
	ll_max_aic = MaxLL(l, TESTS, pis, us, data_end, data_start, ts);
	aic_min = IC(N, parameters, ll_max_aic, ic_type);
      }
      pis[0].F_[row_aic_min][col_aic_min] = 0.0;
      
      // Parameters,DPs,log_lik,AIC,BIC to a file for each model type
      label = std::to_string(n) + "." + std::to_string(parameters);
      out_file << parameters << ","
	       << N          << ","
	       << ll_max_aic << ","
	       << "NaN"      << ","
	       << aic_min    << ","
	       << "NaN"      << ","
	       << label      << std::endl;
      
      out_file_matrix << label << "\t";
      pis[0].Write_("forwards", out_file_matrix);
      
      if(aic_min < aic_current) {
	aic_current = aic_min;
	ll_current_aic = ll_max_aic;
	parameters_current_aic = parameters;
	pi_current_aic = pis[0];
      }
    }
  }
  out_file.close();
  out_file_matrix.close();

  std::vector <std::pair<double, std::string>> ics;
  MakeBestICFromOutputs(out_file_name, out_file_matrix_name, ics, L, number_of_tests, file_name);

  return 0;
}
