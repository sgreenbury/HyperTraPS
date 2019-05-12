#ifndef __TRAJECTORY__HPP
#define __TRAJECTORY__HPP

#include <iostream>
#include <fstream>
#include <vector>
#include "utilities.hpp"
#include "state.hpp"

class Trajectory
{
public:
  State start_;
  State end_;
  State current_;

  std::vector <State> visited_;

  int cutoff_;
  int finished_;
  int n_visited_;
  
  double alpha_;

  int compatibility_count_;
  double cumulative_alpha_;
  
  // TODO: included for adaptation to CT-HTPS
  double beta_;

  int to_health_;
  int to_stopped_;

  Trajectory();
  ~Trajectory();
  Trajectory(const State& start, const State& end);

  Trajectory(const Trajectory& t);
  Trajectory& operator=(const Trajectory& t);

  void Init_();
  void Init_(const State& start, const State& end);
  void Init_(const Trajectory& t);

  void TestUpdateCurrent();

  int ResetFull_(State& new_start, State& new_end);
  void Reset_();
  int ResetStart_(State& new_start);
  int ResetStart_(std::vector <int>& new_start);
  int ResetEnd_(State& new_end);
  int ResetEnd_(std::vector <int>& new_end);
  
  void Clear_();
  
  friend std::ostream& operator<<(std::ostream& out, const Trajectory& t);  
};

std::ostream& operator<<(std::ostream& out, const Trajectory& t);

#endif
