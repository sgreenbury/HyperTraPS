#ifndef __MCMC__HPP
#define __MCMC__HPP

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <random>
#include <chrono>
#include <ctime>

#include "utilities.hpp"
#include "transition.hpp"
#include "trajectory.hpp"
#include "updater_fast.hpp"
#include "state.hpp"
#include "likelihood_sampler_thread.hpp"

class MCMC
{
public:
  std::vector <double> log_likelihoods_;

  std::chrono::time_point<std::chrono::system_clock> start_t_, current_t_;
  std::chrono::duration<double> elapsed_seconds_t_;
  
  std::mt19937 gen_;
  std::vector <std::mt19937> gens_;
  std::vector <std::mt19937> save_gens_;
  std::vector <std::mt19937> proposal_gens_;
  void ProposalGens_(std::vector <std::mt19937>& gens);
  void CopyGens_(std::vector <std::mt19937>& from_gens,
		 std::vector <std::mt19937>& to_gens);

  void SetStartType_(std::string start_type);
  void SetStartPis_(std::vector <Transition>& pis);
  
  std::map <int, Transition> transitions_;
  int N_mcmc_steps_;
  int mcmc_step_;
  int verbose_;
  int vverbose_;
  int burn_in_steps_;
  int iid_step_;
  int seed_;

  int N_threads_;
  
  std::string model_;
  std::string mcmc_type_;
  std::string start_type_;
  std::string time_label_;
  
  // Option to record ordering histogram along with sampling
  // Currently switched off
  int record_ordering_;
  
  // How many steps do we go before performing a sample over the posterior
  // for ordering of traits
  // Currently switched off
  int ordering_sample_gap_;
  
  // Keeps track of the acceptance figures for the chain over whole run
  int accepts0_;
  int accepts0_lower_;
  int rejects0_;
  
  int accepts1_;
  int accepts1_lower_;
  int rejects1_;

  double log_uni_;
  double log_likelihood_difference_;

  int sigma_apm_;
  
  std::string label_;
  std::string out_mcmc_stats_label_;
  std::string out_ordering_label_;

  std::string out_forwards_label_;
  std::string out_backwards_label_;
  std::string out_other_states_label_;
  
  
  // Ouput file for ongoing write
  std::ofstream out_;
  std::ofstream outfile;
  std::ofstream outparams_f;
  std::ofstream outparams_b;
  std::ofstream outparams_o;
  
  std::vector <std::vector <double>> ordering_;
  std::vector <std::vector <double>> ordering_normalised_;

  MCMC(const int N_mcmc_steps, const int burn_in, const std::string label, const int verbose, const std::string model, const std::string mcmc_type, const int iid_step, const std::string prior_type, const int N_threads, const int seed, const std::string out, const double sigma_apm = 0);
  
  void Init_(const int N_mcmc_steps, const int burn_in, const std::string label, const int verbose, const std::string model, const std::string mcmc_type, const int iid_step, const std::string prior_type, const int N_threads, const int seed, const std::string out, const double sigma_apm = 0);

  void PerturbGens_();
  
  void OpenOutfile_();
  void OpenOutParams_(Transition& pi);
  void WriteToOutfile_();
  void WriteToFile_(Transition& pi, int iteration);
  void CloseOutfile_();
  void CloseOutParams_(Transition& pi);
    
  void PerformMCMCStep_(Likelihood& l, std::vector <State>& data, std::vector <Trajectory>& ts, std::vector <Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts);
  
  void PerformMCMCStepAPM_(Likelihood& l, std::vector <State>& data, std::vector <Trajectory>& ts, std::vector <Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts);
  
  void NormaliseOrdering_(int step);

  void RunMCMC_(Likelihood& l, std::vector <State>& data, std::vector <Trajectory>& ts, std::vector <Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts);
  
  friend std::ostream& operator<<(std::ostream& out, const MCMC& m);  
};

std::ostream& operator<<(std::ostream& out, const MCMC& m);

#endif
