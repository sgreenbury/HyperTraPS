#ifndef __LIKELIHOOD_SAMPLER_THREAD__HPP
#define __LIKELIHOOD_SAMPLER_THREAD__HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <thread>
#include <random>

#include "utilities.hpp"
#include "transition.hpp"
#include "state.hpp"
#include "trajectory.hpp"
#include "updater_fast.hpp"

class Likelihood
{
public:
  std::vector <std::vector<double>> alphas_;
  std::vector <std::mt19937> gens_;
  std::vector <std::uniform_real_distribution <double>> diss_;
  
  double likelihood_;
  double total_likelihood_;
  double total_log_likelihood_;
  std::vector <double> log_likelihoods_;
  std::vector <std::thread> threads_;
  
  int N_samples_;
  int N_runs_;
  int verbose_;
  int N_threads_;
  int N_thread_;
  int seperate_pis_;

  Likelihood(const int N_samples, const int N_runs, const int N_threads, const uint32_t SEED, const int seperate_pis=0);
    
  void Init_(const int N_samples, const int N_runs, const int N_threads, const uint32_t SEED, const int seperate_pis=0);

  void Reset_();

  void LogLikelihoodThreads_(std::vector <State>& data, std::vector <Trajectory>& ts, std::vector <Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts, const int thread, const int start_j, const int end_j);

  void EndLikelihoods_(std::vector <State>& data, std::vector <Trajectory>& ts, std::vector<Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts);

  void CalculateTotalLikelihood_();
  void CalculateTotalLogLikelihood_();

  void LoadGens_(std::vector <std::mt19937>& gens);
  void ResetDiss_();
  
  void RunLikelihoodCalculation(std::vector <State>& data, std::vector <Trajectory>& ts, std::vector <Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts);

  
  void RunLikelihoodCalculation(std::vector <State>& data, Trajectory& t, Transition& pi, Updater& u, State& start);
  
  friend std::ostream& operator<<(std::ostream& out, const Likelihood& l);  
};

std::ostream& operator<<(std::ostream& out, const Likelihood& l);

#endif
