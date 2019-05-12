#include "mcmc_sampler.hpp"
#include <unistd.h>

void PrintUsageMessage()
{
  std::string s(100, '-');
  s = s + "\n";
  std::cout << "\n";
  std::cout << s
	    << "\nHyperTraPS Discrete Time - MCMC Sampler\n\n"
	    << s
	    << "\nUsage options:\n";
  std::cout << "  -f <input data file = \"transitions.txt\">\n";
  std::cout << "  -M <model type = \"second-order\">\n";
  std::cout << "     Can take: \"first-order\" or \"second-order\"\n";
  
  std::cout << "  -N <number of posterior samples  = 500>\n";
  
  std::cout << "  -r <HyperTraPS trajectories per transition = 200>\n";
  
  std::cout << "  -n <MCMC start matrix = \"zero\">\n";
  std::cout << "     Can take: \"zero\"\n"
	    << "               \"uniform\"\n";


  std::cout << "  -p <prior range in log space = 20, (-20,20)>\n";
  
  std::cout << "  -k <proposal kernel = \"mcmc\">\n";
  std::cout << "     Can take: \"mcmc\" or \"mcmc-apm\" for auxiliary pseudo marginal mcmc\n";
  
  std::cout << "  -s <sigma for gaussian proposal kernel = 0.025>\n";
  
  std::cout << "  -b <burn in = 50000>\n";
  
  std::cout << "  -i <iid steps between posterior sample draws = 100>\n";
  
  std::cout << "  -S <seed = 0>\n";
  
  std::cout << "  -t <threads or processors = hardware - 1>\n";
  
  std::cout << "  -q <auxiliary pseudo marginal proposal type = 0>\n";
  std::cout << "     Switches between proposing all new random walkers or a single\n"
               "     new random walker out of the -t threads being used\n";
  
  std::cout << "  -L <special load of previously run matrix = \"\">\n";
  std::cout << "  -v <verbosity = 0>\n";
  std::cout << "  -h or -? <print this message>\n";
  std::cout << s;
}

int main(int argc, char *argv[])
{
  // Default flags:
  std::string out_file = "out.txt";
  std::ofstream out(out_file, std::ofstream::out);
  std::string file_label = "transitions.txt";
  std::string model      = "second-order";
  int N_samples          = 5000;
  int runs               = 200;
  std::string start_type = "zero";
  int prior_range        = 20;
  std::string kernel     = "mcmc";
  double sigma           = 0.025;
  int burn_in            = 50000;
  int iid_step           = 100;
  uint32_t seed          = 0;
  int N_threads          = std::thread::hardware_concurrency()-1;
  int sigma_apm          = 0;
  std::string special_load = "";
  int verbose            = 0;

  int opt = 0;
  static const char *optString = "f:M:N:r:n:p:k:s:b:i:S:t:q:L:vh?";
  opt = getopt(argc, argv, optString);
  while(opt!=-1) {
    switch(opt) {
    case 'f':
      file_label = optarg;
      break;
    case 'M':
      model      = optarg;
      break;
    case 'N':
      N_samples  = atoi(optarg);
      break;
    case 'r':
      runs       = atoi(optarg);
      break;
    case 'n':
      start_type = optarg;
      break;
    case 'p':
      prior_range = atoi(optarg);
      break;
    case 'k':
      kernel = optarg;
      break;
    case 's':
      sigma = atof(optarg);
      break;
    case 'b':
      burn_in = atoi(optarg);
      break;
    case 'i':
      iid_step = atoi(optarg);
      break;
    case 'S':
      seed = atoi(optarg);
      break;
    case 't':
      N_threads = atoi(optarg);
      break;
    case 'q':
      sigma_apm = atoi(optarg);
      break;
    case 'L':
      special_load = optarg;
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

  // Apply program
  RandomInitialise(seed);
  int N_mcmc_steps   = N_samples*iid_step + burn_in;
  std::vector <std::vector <int>> d;
  std::string file       = file_label;
  LoadData(file, d);
  int N                  = (int)((double)d.size()/2);
  int L                  = d[0].size();
  std::string label = "";
  
  if(verbose or true) {
    // Print main variable types
    out << "In file: "     + file_label << std::endl
	<< "Data points: " + std::to_string(N) << std::endl
	<< "Model: "       + model << std::endl
	<< "N_samples: "   + std::to_string(N_samples) << std::endl
	<< "Runs: "        + std::to_string(runs) << std::endl
	<< "Start type: "  + start_type << std::endl
	<< "Prior range: " + std::to_string(prior_range) << std::endl
	<< "Kernel: "      + kernel << std::endl
	<< "Sigma: "       + std::to_string(sigma) << std::endl
	<< "Burn in: "     + std::to_string(burn_in) << std::endl
	<< "iid step: "    + std::to_string(iid_step) << std::endl
	<< "Seed: "        + std::to_string(seed) << std::endl
	<< "Threads: "     + std::to_string(N_threads) << std::endl
	<< "Sigma APM: "   + std::to_string(sigma_apm) << std::endl
	<< "Special load: " + special_load << std::endl;
  }
  
  State start(std::vector <int> (L, 0), 0);
  State end(std::vector <int> (L, 1), 0);

  // Load data from d into State vector
  std::vector <State> data_start(N, State(std::vector<int> (L,0),0));
  for(int i=0; i<N; i++) {
    if(verbose)
      out << "Data item end: " << (int)((double)i/2) << std::endl;
    for(unsigned int j=0; j<d[2*i].size(); j++) {
      data_start[i].position_[j] = d[2*i][j];
      if(verbose)
	out << data_start[i].position_[j] << " ";
    }
    if(verbose)
      out << std::endl;
  }
  
  std::vector <State> data_end(N, State(std::vector<int> (L,0),0));
  for(int i=0; i<N; i++) {
    if(verbose)
      out << "Data item end: " << (int)((double)i/2) << std::endl;
    for(unsigned int j=0; j<d[2*i+1].size(); j++) {
      data_end[i].position_[j] = d[2*i+1][j];
      if(verbose)
	out << data_end[i].position_[j] << " ";
    }
    if(verbose)
      out << std::endl;
  }
      
  /* Main likelihood calculation requirements */
  std::vector <Trajectory> ts;
  std::vector <Updater> us;
  std::vector <Transition> pis;
  for(int i=0; i<N_threads; i++) {
    Transition pi(L, sigma, prior_range);
    if(special_load != "")
      pi.Load_(special_load, "forwards");
    Updater u(pi);
    Trajectory t(start, end);
    pis.push_back(pi);
    us.push_back(u);
    ts.push_back(t);
  }

  // Likelihood setup
  Likelihood l(N, runs, N_threads, seed, 0);
  out.close();

  // Set-up mcmc object
  MCMC mcmc(N_mcmc_steps, burn_in, label, verbose, model, kernel, iid_step, start_type, N_threads, seed, out_file, sigma_apm);
  
  // Run mcmc
  mcmc.RunMCMC_(l, data_end, ts, pis, us, data_start);

  return 0;
}
