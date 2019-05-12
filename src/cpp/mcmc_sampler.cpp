#include "mcmc_sampler.hpp"

MCMC::MCMC(const int N_mcmc_steps, const int burn_in, const std::string label, const int verbose, const std::string model, const std::string mcmc_type, const int iid_step, const std::string start_type, const int N_threads, const int seed, const std::string out, const double sigma_apm)
{
  Init_(N_mcmc_steps, burn_in, label, verbose, model, mcmc_type, iid_step, start_type, N_threads, seed, out, sigma_apm);
}

void MCMC::Init_(const int N_mcmc_steps, const int burn_in, const std::string label, const int verbose, const std::string model, const std::string mcmc_type, const int iid_step, const std::string start_type, const int N_threads, const int seed, const std::string out, const double sigma_apm)
{
  out_.open(out, std::ofstream::app);
  N_mcmc_steps_         = N_mcmc_steps;
  mcmc_step_            = 0;
  burn_in_steps_        = burn_in;
  verbose_              = verbose;
  vverbose_             = 0;
  
  record_ordering_      = 0;
  iid_step_             = iid_step;

  model_                = model;
  mcmc_type_            = mcmc_type;

  sigma_apm_            = sigma_apm;
  
  seed_                 = seed;
  N_threads_            = N_threads;

  gen_.seed(seed_);
  for(int i=0; i<N_threads; i++) {
    std::mt19937 gen;
    gen.seed(seed_+i+1);
    gens_.push_back(gen);
    save_gens_.push_back(gen);
    proposal_gens_.push_back(gen);
  }
  
  // Stats about the MCMC run
  accepts0_         = 0;
  accepts0_lower_   = 0;
  rejects0_         = 0;
  accepts1_         = 0;
  accepts1_lower_   = 0;
  rejects1_         = 0;

  log_uni_         = 0.0;
  log_likelihood_difference_ = 0.0;
  
  log_likelihoods_.reserve(N_mcmc_steps_);

  start_type_       = start_type;

  time_t rawtime;
  struct tm* tm;
  char buffer[80];
  time (&rawtime);
  tm = localtime(&rawtime);
  std::strftime(buffer, sizeof(buffer), "%Y%m%d", tm);

  out_ << "Label: " << label << std::endl;
  if(label != "") {
    time_label_ = (std::string(buffer)).substr(2);

    out_mcmc_stats_label_ = time_label_ + "_out_mcmc_sampler_stats"    + label + ".csv";
    out_ordering_label_   = time_label_ + "_out_mcmc_sampler_ordering" + label + ".csv";
    label_                = label.substr(1,label.size());

    out_forwards_label_       = time_label_ + "_out_mcmc_sampler_forwards" + label + ".csv";
    out_backwards_label_      = time_label_ + "_out_mcmc_sampler_backwards" + label + ".csv";
    out_other_states_label_   = time_label_ + "_out_mcmc_sampler_other-states" + label + ".csv";
  }
  else {
    out_mcmc_stats_label_     = "stats.csv";
    out_ordering_label_       = "ordering.csv";
    out_forwards_label_       = "forwards.txt";
    out_backwards_label_      = "backwards.txt";
    out_other_states_label_   = "other-states.txt";
  }
}

void MCMC::SetStartType_(std::string start_type)
{
  start_type_ = start_type;
}

void MCMC::SetStartPis_(std::vector <Transition>& pis)
{
  // Start types: zero, uniform, uniform_basal, uniform_basal_noise_interaction, noise_basal, biased
  if(model_ == "second-order" or model_ == "second-order-lu" or model_ == "second-order-u") {
    if(start_type_ == "zero")
      ;
    else
      if(start_type_ == "uniform")
	for(unsigned int i=0; i<pis.size(); i++)
	  pis[i].MakeUniform_(pis[i].F_);
      else
	if(start_type_ == "uniform-basal")
	  for(unsigned int i=0; i<pis.size(); i++)
	    pis[i].MakeUniformBasal_(pis[i].F_);
	else
	  if(start_type_ == "uniform-basal-noise-interaction")
	    for(unsigned int i=0; i<pis.size(); i++)
	      pis[i].MakeUniformBasalNoiseInteraction_(pis[i].F_);
	  else
	    if(start_type_ == "noise-basal")
	      for(unsigned int i=0; i<pis.size(); i++)
		pis[i].MakeNoiseBasal_(pis[i].F_);
	    else
	      if(start_type_ == "biased")
		for(unsigned int i=0; i<pis.size(); i++)
		  pis[i].MakeBiased_(pis[i].F_);
	      else
		if(start_type_ == "noise")
		  for(unsigned int i=0; i<pis.size(); i++)
		    pis[i].MakeNoise_(pis[i].F_);
    
  }
  else {
    if(model_ == "first-order") {
      if(start_type_ == "zero")
	;
      else
	if(start_type_ == "uniform-basal" or start_type_ == "uniform")
	  for(unsigned int i=0; i<pis.size(); i++)
	    pis[i].MakeUniformBasal_(pis[i].F_);
	else
	  if(start_type_ == "noise-basal" or start_type_ == "noise")
	    for(unsigned int i=0; i<pis.size(); i++)
	      pis[i].MakeNoiseBasal_(pis[i].F_);
    }
  }  
}

void MCMC::WriteToFile_(Transition& pi, int iteration)
{
  if(pi.F_on_) {
    outparams_f << iteration << "\t";
    pi.Write_("forwards", outparams_f);
  }
  if(pi.B_on_) {
    outparams_b << iteration << "\t";
    pi.Write_("backwards", outparams_b);
  }
  if(pi.B_on_) {
    outparams_o << iteration << "\t";
    pi.Write_("other_states", outparams_o);
  }
}

void MCMC::PerformMCMCStep_(Likelihood& l, std::vector <State>& data, std::vector <Trajectory>& ts, std::vector <Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts)
{
  double log_proposal_quotient = 0.0;
  if(mcmc_step_ == 0) {
    l.RunLikelihoodCalculation(data, ts, pis, us, starts);
    log_likelihoods_.push_back(l.total_log_likelihood_);
    accepts1_++;
  }
  else {
    if(model_ == "second-order") {
      pis[0].Perturb_();
    }
    else
      if(model_ == "first-order") {
	pis[0].PerturbDiagonal_();
      }
      else
	if(model_ == "second-order-lu"){
	  log_proposal_quotient = pis[0].PerturbLU_();
	}
	else
	  if(model_ == "second-order-u"){
	    pis[0].PerturbU_();
	  }
        
    // Proposed
    if(verbose_) {
      out_ << "Proposed: " << std::endl;
      for(unsigned int i=0; i<pis.size(); i++)
	out_ << pis[i];
    }
    
    // Previous log likelihood
    double ll_prev       = log_likelihoods_.back();

    // Calculate proposed log likelihood
    l.RunLikelihoodCalculation(data, ts, pis, us, starts);
    double ll_new        = l.total_log_likelihood_;
    
    log_uni_ = log(Uniform(gen_));
    log_likelihood_difference_ = ll_new - ll_prev + log_proposal_quotient;
    
    if(verbose_)
      out_ << ll_new << "\t"
		<< ll_prev << "\t"
		<< log_likelihood_difference_ << "\t"
		<< log_uni_ << "\t" << std::endl;
    
    if(log_likelihood_difference_ > log_uni_) {
      log_likelihoods_.push_back(ll_new);
      accepts1_++;
      if(log_likelihood_difference_ < 0)
	accepts1_lower_++;
    }
    else {
      rejects1_++;
      log_likelihoods_.push_back(ll_prev);
      pis[0].RestorePreviousTransitions_();
    }
  }
}

void MCMC::ProposalGens_(std::vector <std::mt19937>& gens)
{
  if(sigma_apm_ != -1) {
    int rand_i = Uniform()*gens.size();
    for(unsigned int i=0; i<gens.size(); i++) {
      if((sigma_apm_ == 1 and (int)i == rand_i) or sigma_apm_ == 0) {
	gen_();
	gens_[i] = gen_;
      }
    }
  }
}

void MCMC::CopyGens_(std::vector <std::mt19937>& from_gens,
		     std::vector <std::mt19937>& to_gens)
{
  for(unsigned int i=0; i<from_gens.size(); i++) {
    to_gens[i] = from_gens[i];
  }
}

void MCMC::PerformMCMCStepAPM_(Likelihood& l, std::vector <State>& data, std::vector <Trajectory>& ts, std::vector <Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts)
{
  for(int step=0; step<2; step++) {
    if(mcmc_step_ == 0 and step == 0) {
      CopyGens_(gens_, save_gens_);
      CopyGens_(gens_, proposal_gens_);
      l.LoadGens_(gens_);
      l.RunLikelihoodCalculation(data, ts, pis, us, starts);
      log_likelihoods_.push_back(l.total_log_likelihood_);
      accepts0_++;
    }
    else {
      if(step == 0) {
	CopyGens_(gens_, save_gens_);
	ProposalGens_(gens_);
	CopyGens_(gens_, proposal_gens_);
	l.LoadGens_(gens_);
      }
      else {
	if(step == 1) {
	  l.LoadGens_(gens_);
	  if(model_ == "second-order")
	    pis[0].Perturb_();
	  else
	    if(model_ == "first-order")
	      pis[0].PerturbDiagonal_();
	}
      }

      double ll_prev = log_likelihoods_.back();
      l.RunLikelihoodCalculation(data, ts, pis, us, starts);
      double ll_new  = l.total_log_likelihood_;
      
      log_uni_ = log(Uniform());
      log_likelihood_difference_ = ll_new - ll_prev;

      // Added for sigma_apm_ = 1 as proposal now asymmetric
      if(sigma_apm_ == 1 and step == 0) {
	log_likelihood_difference_ += (gens_.size() > 1 ? log((double)1/gens_.size()) : 0);
      }
      
      if(vverbose_) {
	out_ << mcmc_step_ << "\t"
	     << step << "\t"
	     << ll_new << "\t"
	     << ll_prev << "\t"
	     << log_likelihood_difference_ << "\t"
	     << log_uni_ << "\t" << std::endl;
      }
      
      if(log_likelihood_difference_ > log_uni_) {
	if(step == 0) {
	  log_likelihoods_.push_back(ll_new);
	  accepts0_++;
	  if(log_likelihood_difference_ < 0)
	    accepts0_lower_++;
	  CopyGens_(proposal_gens_, gens_);
	}
	if(step == 1) {
	  log_likelihoods_[log_likelihoods_.size()-1] = ll_new;
	  accepts1_++;
	  if(log_likelihood_difference_ < 0)
	    accepts1_lower_++;
	}
      }
      else {
	if(step == 0) {
	  log_likelihoods_.push_back(ll_prev);
	  rejects0_++;
	  CopyGens_(save_gens_, gens_);
	}
	if(step == 1) {
	  rejects1_++;
	  pis[0].RestorePreviousTransitions_();
	}
      }
    }
  }
}

void MCMC::NormaliseOrdering_(int step)
{
  int number_of_samples = (int)((double)(step - burn_in_steps_)/iid_step_) + 1;
  for(unsigned int i=0; i<ordering_.size(); i++)
    for(unsigned int j=0; j<ordering_[i].size(); j++)
      ordering_normalised_[i][j] = ordering_[i][j]/number_of_samples;
}

void MCMC::OpenOutfile_()
{
  outfile.open(out_mcmc_stats_label_, std::ofstream::out);
  start_t_ = std::chrono::system_clock::now();
  current_t_ = start_t_;
  outfile << "mcmc_step,"
	  << "log_likelihood,"
	  << "accepts,"
	  << "accepts_lower,"
	  << "rejects,"
	  << "accepts0,"
	  << "accepts0_lower,"
	  << "rejects0,"
	  << "time"
	  << std::endl;
  outfile.precision(7);
}

void MCMC::OpenOutParams_(Transition& pi)
{
  if(pi.F_on_) {
    outparams_f.open(out_forwards_label_);
    outparams_f.precision(4);
  }
  if(pi.B_on_) {
    outparams_b.open(out_backwards_label_);
    outparams_b.precision(4);
  }
  if(pi.O_on_) {
    outparams_o.open(out_other_states_label_);
    outparams_o.precision(4);
  }
}

void MCMC::CloseOutParams_(Transition& pi)
{
  if(pi.F_on_)
    outparams_f.close();
  if(pi.B_on_)
    outparams_b.close();
  if(pi.O_on_)
    outparams_o.close();
}

void MCMC::WriteToOutfile_()
{
  current_t_ = std::chrono::system_clock::now();
  elapsed_seconds_t_ = current_t_ - start_t_;
  int align = 1;
  outfile << mcmc_step_                                 << ","
	  << log_likelihoods_.back()                    << ","
	  << accepts1_ / (double)(mcmc_step_ + align)       << ","
	  << accepts1_lower_ / (double)(mcmc_step_ + align) << ","
	  << rejects1_ / (double)(mcmc_step_ + align)       << ","
	  << accepts0_ / (double)(mcmc_step_ + align)       << ","
	  << accepts0_lower_ / (double)(mcmc_step_ + align) << ","
	  << rejects0_ / (double)(mcmc_step_ + align)       << ","
    	  << elapsed_seconds_t_.count()/(60*60*24)      << " days"
	  << std::endl;
}

void MCMC::CloseOutfile_()
{
  out_.close();
  outfile.close();
}

void MCMC::RunMCMC_(Likelihood& l, std::vector <State>& data, std::vector <Trajectory>& ts, std::vector <Transition>& pis, std::vector <Updater>& us, std::vector <State>& starts)
{
  SetStartPis_(pis);

  // Set-up recording of ordering.
  if(record_ordering_ == 1) {
    ordering_.resize(pis[0].L_, std::vector <double> (pis[0].L_, 0.0));
    ordering_normalised_.resize(pis[0].L_, std::vector <double> (pis[0].L_, 0.0));
  }

  State start(std::vector <int> (pis[0].L_, 0), 0);
  State end(std::vector <int> (pis[0].L_, 1), 0);

  OpenOutfile_();
  OpenOutParams_(pis[0]);
  for(mcmc_step_=0; mcmc_step_ < N_mcmc_steps_; mcmc_step_++) {
    if(mcmc_type_ == "mcmc")
      PerformMCMCStep_(l, data, ts, pis, us, starts);
    else
      if(mcmc_type_ == "mcmc-apm")
	PerformMCMCStepAPM_(l, data, ts, pis, us, starts);
    
    // Write data to outfile
    if(verbose_) {
      out_ << "Step: " << mcmc_step_
	   << "\t" << (double)accepts1_/mcmc_step_
	   << "\t" << (double)rejects1_/mcmc_step_
	   << "\t" << (double)accepts1_lower_/mcmc_step_
	   << "NaN"
	   << std::endl;
      out_ << "Step: " << mcmc_step_ << std::endl;
      out_ << l;
    }
    
    // Write stats every iteration and transition matrix to an out file
    WriteToOutfile_();
    if(mcmc_step_ % iid_step_ == 0) {
      WriteToFile_(pis[0], mcmc_step_);
      if(verbose_)
	out_ << ts[0];
    }
    // TODO: Remove as currently non-functional
    if(mcmc_step_ >= burn_in_steps_ and record_ordering_ == 1 and mcmc_step_ % iid_step_ == 0) {
      us[0].PerformWalksFromZeroToOne(ts[0], pis[0], ordering_, start, end);
      NormaliseOrdering_(mcmc_step_+1);
      WriteMatrixCSV(ordering_normalised_, out_ordering_label_);
    }
  }
  CloseOutfile_();
  CloseOutParams_(pis[0]);
}

std::ostream& operator<<(std::ostream& out, const MCMC& m)
{
  for(unsigned int i=0; i<m.log_likelihoods_.size(); i++)
    out << i << "," << m.log_likelihoods_[i] << "\n";
  return out;
}
