#ifndef __UPDATER_FAST__HPP
#define __UPDATER_FAST__HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <limits>

#include "utilities.hpp"
#include "transition.hpp"
#include "trajectory.hpp"

class Updater
{
public:

  int CalculateForwardsWeights_(Trajectory& t, Transition& pi, int compatible);
  int UpdateForwardsWeights_(Trajectory& t, Transition& pi, int compatible);

  // TODO: not currently functional until trajectory and transition
  //       classes have backward and other well defined
  int CalculateBackwardsWeights_(Trajectory& t, Transition& pi, int compatible);
  int CalculateOtherStatesWeights_(Trajectory& t, Transition& pi, int compatible);
  int CalculateTotalWeight_(Trajectory& t, Transition& pi, int compatible);
  int UpdateTotalWeight_(Trajectory& t, Transition& pi, int compatible);
  
  int UpdateAlpha_(Trajectory& t);
  int UpdateCProbabilityVectors_(Transition& pi);

  int UpdateCompatibilityCountAndCumulativeAlpha_(Trajectory& t, const int incomplete_variable);
  
  // Perform move based on cumulative probability vectors
  int MovePossible_();
  int MoveForwards_(Trajectory& t, int pos);
  int MoveBackwards_(Trajectory& t, int pos);
  int MoveOtherStates_(Trajectory& t, int pos);
 
  int UpdateTrajectory_(Trajectory& t, Transition& pi);
  int UpdateTrajectory_(Trajectory& t, Transition& pi, std::mt19937& gen);
  int UpdateTrajectory_(Trajectory& t, Transition& pi, std::mt19937& gen, std::uniform_real_distribution <double>& dis);

  // TODO: non-functional method for now as related to incomplete data problem
  int UpdateTrajectoryToCompatible_(Trajectory& t,
  				    Transition& pi,
  				    int must_move);

  void PerformWalkUntilCompatible_(Trajectory& t,
				   Transition& pi,
				   int must_move);


  int PerformWalksFromZeroToOne(Trajectory& t,
				Transition& pi,
				std::vector <std::vector <double>>& ordering,
				State& start,
				State& end);

  int PerformWalksFromZeroToOne(Trajectory& t,
				Transition& pi,
				std::vector <std::vector <double>>& ordering,
				State& start,
				State& end,
				std::ofstream& f);
  

  
  
  void PerformWalk(Trajectory& t,
		   Transition& pi,
		   State& start,
		   State& end,
		   int start_count,
		   int end_count);
  
  State PerformWalk(Trajectory& t,
		    Transition& pi,
		    std::vector <std::vector <double>>& ordering,
		    State& start,
		    State& end,
		    std::ofstream& f,
		    int start_count,
		    int end_count,
		    const std::string record);

  State PerformWalk(Trajectory& t,
		    Transition& pi,
		    std::vector <std::vector <double>>& ordering,
		    State& start,
		    State& end,
		    std::ofstream& f,
		    int start_count,
		    int end_count,
		    std::map <std::pair<long long, long long>, int>& transition_counter,
		    const std::string record);
  
  State PerformWalk(Trajectory& t,
		    Transition& pi,
		    std::vector <std::vector <double>>& ordering,
		    State& start,
		    State& end,
		    std::ofstream& f,
		    int start_count,
		    int end_count,
		    std::map <std::pair<long long, long long>, int>& transition_counter,
		    std::map <long long, int>& state_counter,
		    const std::string record);
    
    
  Updater(Transition& pi);
  void Init_(Transition& pi);

  Updater(const Updater& u);
  Updater& operator=(const Updater& u);

  void Init_(const Updater& u);

  void Clear_();

  double count_;
  double small_positive_;
  double max_count_;
  
  int    TESTS_;
  
  double forwards_weight_;
  double backwards_weight_;
  double other_states_weight_;
  
  double total_weight_;
  double total_compatible_weight_;

  int previous_feature_change_;
  int previous_direction_;
  
  std::vector <double> forwards_probabilities_;
  std::vector <double> forwards_compatible_probabilities_;
  std::vector <double> cforwards_probabilities_;


  // TODO: not currently functional until trajectory and transition
  //       classes have backward and other well defined
  std::vector <double> backwards_probabilities_;
  std::vector <double> backwards_compatible_probabilities_;
  std::vector <double> cbackwards_probabilities_;
  
  std::vector <double> other_states_probabilities_;
  std::vector <double> cother_states_probabilities_;


  friend std::ostream& operator<<(std::ostream& os, const Updater& u);  
};

std::ostream& operator<<(std::ostream& out, const Updater& u);

#endif
