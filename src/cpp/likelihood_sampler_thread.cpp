#include "likelihood_sampler_thread.hpp"

Likelihood::Likelihood(const int N_samples, const int N_runs, const int N_threads, const uint32_t SEED, const int seperate_pis)
{
  Init_(N_samples, N_runs, N_threads, SEED, seperate_pis);
}

void Likelihood::Init_(const int N_samples, const int N_runs, const int N_threads, const uint32_t SEED, const int seperate_pis)
{
  N_samples_    = N_samples;
  N_runs_       = N_runs;
  N_threads_    = N_threads;
  N_thread_     = (double)N_samples_/(N_threads_ - 1);
  seperate_pis_ = seperate_pis;
  
  for(int i=0; i<N_threads_; i++) {
    std::mt19937 generator;
    std::uniform_real_distribution <double> dis(0,1);
    generator.seed(SEED+1);
    gens_.push_back(generator);
    diss_.push_back(dis);
  }
  
  verbose_   = 0;
  for(int i=0; i<N_threads_; i++)
    alphas_.push_back(std::vector <double> (N_samples_, 0.0));

  threads_.reserve(N_threads_);
  log_likelihoods_.reserve(N_samples_);

  Reset_();
}

void Likelihood::Reset_()
{
  log_likelihoods_.clear();
  total_likelihood_     = 1.0;
  total_log_likelihood_ = 0.0;
  likelihood_           = 0.0;
}

void Likelihood::LoadGens_(std::vector <std::mt19937>& gens)
{
  for(unsigned int i=0; i<gens.size(); i++) {
    gens_[i] = gens[i];
  }
  ResetDiss_();
}

void Likelihood::ResetDiss_()
{
  for(unsigned int i=0; i<gens_.size(); i++) {
    diss_[i].reset();
  } 
}

void Likelihood::LogLikelihoodThreads_(std::vector <State>& data,
				       std::vector <Trajectory>& ts,
				       std::vector <Transition>& pis,
				       std::vector <Updater>& us,
				       std::vector <State>& starts,
				       const int thread,
				       const int start_j,
				       const int end_j)
{
  for(int j=start_j; j<end_j; j++) {
    alphas_[thread].resize(N_runs_);
    ts[thread].Reset_();
    ts[thread].ResetEnd_(data[j]);

    // Run likelihoods for each thread
    for(int i=0; i<N_runs_; i++) {
      ts[thread].ResetStart_(starts[j]);
      ts[thread].Reset_();
      us[thread].Clear_();
      int finished = 0;
      while(finished == 0) {
	if(seperate_pis_)
	  finished = us[thread].UpdateTrajectory_(ts[thread], pis[thread], gens_[thread], diss_[thread]);
	else
	  finished = us[thread].UpdateTrajectory_(ts[thread], pis[0], gens_[thread], diss_[thread]);
      }
      alphas_[thread][i] = ts[thread].cumulative_alpha_;
    }
    // Calculate average log likelihood from all N_runs_
    double mean_alpha = Mean(alphas_[thread]);
    if(std::isinf(log(mean_alpha)))
      log_likelihoods_[j] = log(std::numeric_limits<double>::min());
    else
      log_likelihoods_[j] = log(mean_alpha);
  }
}

// TODO: Add alternative code for when threads are equal to 1 to avoid calling threads
void Likelihood::EndLikelihoods_(std::vector <State>& data, std::vector <Trajectory>& ts, std::vector<Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts)
{
  Reset_();
  log_likelihoods_.resize(N_samples_, 0.0);
  threads_.clear();
  for(int i=0; i<N_threads_; i++) {
    int start_j = i*N_thread_;
    int end_j   = (i == N_threads_ - 1 ? data.size() : (i+1)*N_thread_);
    threads_.push_back(std::thread(&Likelihood::LogLikelihoodThreads_,
				   this,
				   std::ref(data),
				   std::ref(ts),
				   std::ref(pis),
				   std::ref(us),
				   std::ref(starts),
				   i,
				   start_j,
				   end_j));
				   
  }
  for(unsigned int i=0; i<threads_.size(); i++)
    threads_[i].join();
}

void Likelihood::CalculateTotalLikelihood_()
{
  total_likelihood_ = 1.0;
  for(unsigned int i=0; i<log_likelihoods_.size(); i++) {
    total_likelihood_ *= exp(log_likelihoods_[i]);
  }

  // TODO: Catch exception
  if(std::isinf(total_likelihood_))
    total_likelihood_ = 0.0;
}

void Likelihood::CalculateTotalLogLikelihood_()
{
  total_log_likelihood_ = 0.0;
  for(unsigned int i=0; i<log_likelihoods_.size(); i++) {
    total_log_likelihood_ += log_likelihoods_[i];
  }

  // TODO: Consider change to -INFINITY
  if(std::isinf(total_log_likelihood_))
    total_log_likelihood_ = log(std::numeric_limits<double>::min());
}

void Likelihood::RunLikelihoodCalculation(std::vector <State>& data, std::vector <Trajectory>& ts, std::vector <Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts)
{
  Reset_();
  EndLikelihoods_(data, ts, pis, us, starts);
  CalculateTotalLikelihood_();
  CalculateTotalLogLikelihood_();
}

std::ostream& operator<<(std::ostream& out, const Likelihood& l)
{
  out << "Individual sample likelihoods:\n";
  for(unsigned int i=0; i<l.log_likelihoods_.size(); i++)
    out << i << "\t" << l.log_likelihoods_[i] << "\t" << exp(l.log_likelihoods_[i]) << std::endl;
  out << "Total likelihood: " << l.total_likelihood_ << std::endl;
  out << "Total log likelihood: " << l.total_log_likelihood_ << std::endl;
  return out;
}
