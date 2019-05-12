#include "trajectory.hpp"

Trajectory::Trajectory()
{
  Init_();
}

Trajectory::Trajectory(const State& start, const State& end)
{
  Init_(start, end);
}


Trajectory::Trajectory(const Trajectory& t)
{
  Init_(t);
}

Trajectory::~Trajectory()
{
  Clear_();
}

void Trajectory::Init_()
{
  // Arbirtrary limits for number of states to visit
  cutoff_   = 1000;
  finished_ = 0;

  compatibility_count_ = 0;
  cumulative_alpha_    = 0.0;

  // For likelihood with discrete time HyperTraPS
  alpha_    = 1.0;
  // TODO: for adaptation to CT-HTPS with rates
  beta_     = 0.0;

  // TODO: Additional off hypercube moves currently not in use
  to_health_  = 0;
  to_stopped_ = 0;
}


void Trajectory::Init_(const State& start, const State& end)
{
  start_    = start;
  end_      = end;
  current_  = start_;

  Init_();

  visited_.resize(cutoff_, State (std::vector <int> (start_.position_.size(), 0), 0));

  n_visited_ = 0;
}

void Trajectory::Init_(const Trajectory& t)
{
  start_   = t.start_;
  end_     = t.end_;
  current_ = t.current_;

  cutoff_ = t.cutoff_;
  finished_ = t.finished_;
  
  alpha_ = t.alpha_;
  compatibility_count_ = t.compatibility_count_;
  cumulative_alpha_    = t.cumulative_alpha_;
  beta_  = t.beta_;

  n_visited_ = t.n_visited_;

  visited_.resize(t.visited_.size());

  for(unsigned int i=0; i<t.visited_.size(); i++)
    visited_[i] = t.visited_[i];  
}


Trajectory& Trajectory::operator=(const Trajectory &t)
{
  if( &t != this ) {
    Clear_();
    Init_( t );
  }
  return *this;
}

void Trajectory::TestUpdateCurrent()
{
  current_.position_[0] = 1;
  visited_[n_visited_] = current_;
  n_visited_++;
}

void Trajectory::Reset_()
{
  current_             = start_;
  finished_            = 0;
  alpha_               = 1.0;
  compatibility_count_ = 0;
  cumulative_alpha_    = 0.0;
  n_visited_ = 0;
}

int Trajectory::ResetStart_(State& new_start)
{
  start_ = new_start;
  return 0;
}

int Trajectory::ResetStart_(std::vector <int>& new_start)
{
  start_.position_ = new_start;
  return 0;
}

int Trajectory::ResetEnd_(State& new_end)
{
  end_ = new_end;
  return 0;
}

int Trajectory::ResetFull_(State& new_start, State& new_end)
{
  ResetStart_(new_start);
  ResetEnd_(new_end);
  Reset_();
  return 0;
}

int Trajectory::ResetEnd_(std::vector <int>& new_end)
{
  end_.position_ = new_end;
  return 0;
}

void Trajectory::Clear_()
{
  return;
}


std::ostream& operator<< (std::ostream& out, Trajectory const& t)
{
  out << "Trajectory:\n";
  out << "Start:\n";
  out << t.start_;
  out << "Current:\n";
  out << t.current_;
  out << "End:\n";
  out << t.end_;
  out << "Alpha: " << t.alpha_ << "\n";
  out << "Compatibility count: " << t.compatibility_count_ << "\n";
  out << "Cumulative alpha: " << t.cumulative_alpha_ << "\n";
  out << "Visited:\n";
  for(int i=0; i<t.n_visited_; i++)
    out << t.visited_[i];
  out << std::endl;
  return out;
}
