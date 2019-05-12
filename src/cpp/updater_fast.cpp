#include "updater_fast.hpp"

Updater::Updater(Transition& pi)
{
  Init_(pi);
}

Updater::Updater(const Updater& u)
{
  Init_(u);
}

void Updater::Init_(Transition& pi)
{
  small_positive_     = std::numeric_limits <double>::min();

  // Maximum value that an exponential weight can take
  // This is hard-wired as 1e-4 smaller than max.
  max_count_          = 0.0001*std::numeric_limits <double>::max();
  
  count_              = 0.0;
  TESTS_              = 100;
  
  forwards_weight_    = 0.0;
  backwards_weight_   = 0.0;
  other_states_weight_= 0.0;
  
  total_weight_            = 0.0;
  total_compatible_weight_ = 0.0;

  previous_feature_change_ = -2;
  previous_direction_ = -2;
  
  forwards_probabilities_.resize                (pi.F_.size(), 0.0);
  forwards_compatible_probabilities_.resize     (pi.F_.size(), 0.0);
  cforwards_probabilities_.resize               (pi.F_.size(), 0.0);

  backwards_probabilities_.resize               (pi.B_.size(), 0.0);
  backwards_compatible_probabilities_.resize    (pi.B_.size(), 0.0);
  cbackwards_probabilities_.resize              (pi.B_.size(), 0.0);

  other_states_probabilities_.resize            (pi.O_.size(), 0.0);
  cother_states_probabilities_.resize           (pi.O_.size(), 0.0);
}

void Updater::Init_(const Updater& u)
{
  count_                                  = u.count_;
  small_positive_                         = u.small_positive_;
  max_count_                              = u.max_count_;
  TESTS_                                  = u.TESTS_;

  forwards_weight_                        = u.forwards_weight_;
  backwards_weight_                       = u.backwards_weight_;
  other_states_weight_                    = u.other_states_weight_;
  
  total_weight_                           = u.total_weight_;
  total_compatible_weight_                = u.total_compatible_weight_;

  previous_feature_change_                = u.previous_feature_change_;
  previous_direction_                     = u.previous_direction_;

  forwards_probabilities_                 = u.forwards_probabilities_;
  forwards_compatible_probabilities_      = u.forwards_compatible_probabilities_;
  cforwards_probabilities_                = u.cforwards_probabilities_;

  backwards_probabilities_                = u.backwards_probabilities_;
  backwards_compatible_probabilities_     = u.backwards_compatible_probabilities_;
  cbackwards_probabilities_               = u.cbackwards_probabilities_;
  other_states_probabilities_             = u.other_states_probabilities_;
  cother_states_probabilities_            = u.cother_states_probabilities_;
}


Updater& Updater::operator=(const Updater &u)
{
  if( &u != this ) {
    Clear_();
    Init_(u);
  }
  return *this;
}

int Updater::CalculateForwardsWeights_(Trajectory& t, Transition& pi, int compatible)
{  
  forwards_weight_ = 0.0;
  if(pi.F_on_) {
    for(unsigned int j=0; j<pi.F_features_.size(); j++) {
      count_ = 0.0;
      if(pi.F_features_[j]) {
	if((t.current_.position_[j] == 0 and
	    t.end_.position_[j] != 0 and
	    compatible == 1)
	   or
	   (t.current_.position_[j] == 0 and
	    compatible == 0)
	   ) {
	  count_ = pi.F_[j][j];
	  
	  for(unsigned int i=0; i<pi.F_.size(); i++)
	    count_ += t.current_.position_[i] * pi.F_[j][i];
	  count_   = std::min(std::max(exp(count_),2*small_positive_), max_count_);
	}
	if(compatible == 0)
	  forwards_probabilities_[j] = count_;
	if(compatible == 1)
	  forwards_compatible_probabilities_[j] = count_;
      }
      forwards_weight_ += count_;
    }
  }
  return 0;
}

int Updater::UpdateForwardsWeights_(Trajectory& t, Transition& pi, int compatible)
{
  forwards_weight_ = 0.0;
  if(pi.F_on_) {
    // Given the element that has changed, update forwards weights
    // If previous step was forwards:
    //   - For each j:
    //     - If el that has changed, set element that has changed to 0.0
    //     - multiply forwards_probabilities by: exp(pi.F_[j][previous_element_change])
    //     - Add to forwards weight
    // If previous step was backwards (currently does not occur):
    //   - For each j:
    //     - If el that has changed, calculate element fully
    //     - otherwise divide by exp(pi.F_[j][previous_element_change])
    // Check size of new weight
    // Add to forwards weight.
    for(unsigned int j=0; j<pi.F_features_.size(); j++) {
      if(compatible == 0)
	count_ = forwards_probabilities_[j];
      if(compatible == 1)
	count_ = forwards_compatible_probabilities_[j];

      if(pi.F_features_[j]) {
	if((t.current_.position_[j] == 0 and
	    t.end_.position_[j] != 0 and
	    compatible == 1)
	   or
	   (t.current_.position_[j] == 0 and
	    compatible == 0)
	   ) {

	  // Previous step forwards
	  if(previous_direction_ == 1) {
	    if((int)j == previous_feature_change_)
	      count_ = 0.0;
	    else
	      count_ *= exp(pi.F_[j][previous_feature_change_]);
	  }
	  // Previous step backwards
	  if(previous_direction_ == -1) {
	    if((int)j == previous_feature_change_) {
	      count_ = 0.0;
	      for(unsigned int i=0; i<pi.F_.size(); i++)
		count_ += t.current_.position_[i] * pi.F_[j][i];
	    }
	    else
	      count_ /= exp(pi.F_[j][previous_feature_change_]);
	  }
	  count_   = std::min(std::max(count_,2*small_positive_), max_count_);
	}
	else
	  count_ = 0.0;
	
	forwards_weight_ += count_;
      }

      // Update the probabilities based on whether compatible or not.
      if(compatible == 0) {
	forwards_probabilities_[j] = count_;
      }
      if(compatible == 1) {
	forwards_compatible_probabilities_[j] = count_;
      }
    }
  }
  return 0;
}


// TODO: not currently functional until trajectory and transition
//       classes have backward and other well defined
int Updater::CalculateBackwardsWeights_(Trajectory& t, Transition& pi, int compatible)
{
  backwards_weight_ = 0.0;
  if(pi.B_on_) {
    for(unsigned int j=0; j<pi.B_features_.size(); j++) {
      count_ = 0.0;
      if(pi.B_features_[j]) {
	if((t.current_.position_[j] == 1 and
	    t.end_.position_[j] != 1 and
	    compatible == 1)
	   or
	   (t.current_.position_[j] == 1 and
	    compatible == 0)
	   ) {
	  count_ = pi.B_[j][j];
	  for(unsigned int i=0; i<pi.B_.size(); i++) {
	    count_ += t.current_.position_[i] * pi.B_[j][i];
	  }
	  count_   = std::min(std::max(exp(count_),2*small_positive_), max_count_);
	}
	backwards_weight_ += count_;
      }
      if(compatible == 0)
	backwards_probabilities_[j] = count_;

      if(compatible == 1)
	backwards_compatible_probabilities_[j] = count_;
    }
  }
  return 0;
}

// TODO: not currently functional until trajectory and transition
//       classes have backward and other well defined
int Updater::CalculateOtherStatesWeights_(Trajectory& t, Transition& pi, int compatible)
{
  other_states_weight_ = 0.0;
  if(pi.O_on_) {
    for(unsigned int j=0; j<pi.O_features_.size(); j++) {
      count_ = 0.0;
      if(pi.O_features_[j] and j==0 and t.to_health_ == 1) {
	count_ = pi.O_[j][j];
	for(unsigned int i=0; i<pi.O_.size(); i++) {
	  count_ += t.current_.position_[i] * pi.O_[j][i];
	}
	count_   = std::min(std::max(exp(count_),2*small_positive_), max_count_);
      }
      // Handling transitions to stopped/death state
      if(pi.O_features_[j] and j==1 and t.to_stopped_ == 1) {
	count_ = pi.O_[j][j];
	for(unsigned int i=0; i<pi.O_.size(); i++) {
	  count_ += t.current_.position_[i] * pi.O_[i][j];
	}
	count_   = std::min(std::max(exp(count_),2*small_positive_), max_count_);
      }
      other_states_weight_ += count_;
      other_states_probabilities_[j] = count_; 
    }
  }
  return 0;
}

int Updater::CalculateTotalWeight_(Trajectory& t, Transition& pi, int compatible)
{
  // Sum over each possible type of move to contribute to total_weight
  CalculateForwardsWeights_(t, pi, compatible);
  CalculateBackwardsWeights_(t, pi, compatible);
  CalculateOtherStatesWeights_(t, pi, compatible);
  if(compatible)
    total_compatible_weight_ = forwards_weight_ + backwards_weight_ + other_states_weight_;
  else
    total_weight_ = forwards_weight_ + backwards_weight_ + other_states_weight_;
  return 0;
}

int Updater::UpdateTotalWeight_(Trajectory& t, Transition& pi, int compatible)
{
  // Sum over each possible type of move to contribute to total_weight
  UpdateForwardsWeights_(t, pi, compatible);
  CalculateBackwardsWeights_(t, pi, compatible);
  CalculateOtherStatesWeights_(t, pi, compatible);
  if(compatible)
    total_compatible_weight_ = forwards_weight_ + backwards_weight_ + other_states_weight_;
  else
    total_weight_ = forwards_weight_ + backwards_weight_ + other_states_weight_;
  return 0;
}

int Updater::UpdateCProbabilityVectors_(Transition& pi)
{
  // Calculate the probability of making forwards or backwards moves
  if(pi.F_on_) {
    for(unsigned int i=0; i<forwards_probabilities_.size(); i++)
      cforwards_probabilities_[i] = (i == 0 ? 
				     0.0 :
				     cforwards_probabilities_[i-1])
	+
	forwards_compatible_probabilities_[i]
	/total_compatible_weight_;
  }
  
  if(pi.B_on_) {
    for(unsigned int i=0; i<backwards_probabilities_.size(); i++)
      cbackwards_probabilities_[i] = (i == 0 ?
				     0.0 :
				     cbackwards_probabilities_[i-1])
	+
	backwards_probabilities_[i]/total_compatible_weight_;
  }

  if(pi.O_on_) {
    for(unsigned int i=0; i<other_states_probabilities_.size(); i++)
      cbackwards_probabilities_[i] = (i == 0 ?
				     0.0 :
				     cother_states_probabilities_[i-1])
	+
	other_states_probabilities_[i]/total_compatible_weight_;
  }
  return 0;
}

int Updater::MovePossible_()
{
  // Return:
  // 0 if compatible move possible,
  // 1 if only non-compatible move possible,
  //-1 if no move possible.
  
  if(total_weight_ > small_positive_ and total_compatible_weight_ < small_positive_)
    return 1;
  if(total_weight_ > small_positive_ and total_compatible_weight_ > small_positive_)
    return 0;

  return -1;
}

int Updater::MoveForwards_(Trajectory& t, int pos)
{
  t.visited_[t.n_visited_] = t.current_;
  t.n_visited_++;
  t.current_.position_[pos] = 1;
  return 0;
}

// TODO: not currently functional until trajectory and transition
//       classes have backward and other well defined
int Updater::MoveBackwards_(Trajectory& t, int pos)
{
  t.visited_[t.n_visited_] = t.current_;
  t.n_visited_++;
  t.current_.position_[pos] = 0;
  return 0;
}

// TODO: not currently functional until trajectory and transition
//       classes have backward and other well defined
int Updater::MoveOtherStates_(Trajectory& t, int pos)
{
  t.visited_[t.n_visited_] = t.current_;
  t.n_visited_++;
  if(pos == 0)
    // Reset to 0^L state
    for(unsigned int i=0; i<t.current_.position_.size(); i++)
      t.current_.position_[i] = 0;
  if(pos == 1)
    // Move to stopped state
    t.current_.stopped_ = 1;

  return 0;
}

int Updater::UpdateAlpha_(Trajectory& t)
{
  // Sum over each possible type of move to contribute to total_weight
  t.alpha_ *= total_compatible_weight_/total_weight_;
  return 0;
}

int Updater::UpdateCompatibilityCountAndCumulativeAlpha_(Trajectory& t,
							 const int incomplete_variable)
{
  // Sum over each possible type of move to contribute to total_weight
  if(t.current_.IsCompatible(t.end_, incomplete_variable) == 1) {
    t.compatibility_count_ += 1;
    t.cumulative_alpha_ += t.alpha_;
  }
  return 0;
}

int Updater::UpdateTrajectory_(Trajectory& t, Transition& pi, std::mt19937& gen)
{
  // Incomplete data currently represented as a 2.
  int incomplete_variable = 2;

  // 1. Calculations for compatibility and weight of edges for possible moves.
  if(previous_feature_change_ == -2) {
    CalculateTotalWeight_(t, pi, 0);
    CalculateTotalWeight_(t, pi, 1);
  }
  else {
   UpdateTotalWeight_(t, pi, 0);
   UpdateTotalWeight_(t, pi, 1);
  } 
  UpdateCProbabilityVectors_(pi);

  // 2. With known transition possibilities, count compatibility of current state.
  UpdateCompatibilityCountAndCumulativeAlpha_(t, incomplete_variable);

  // 3. Can any move be made.
  if(t.current_.IsSame(t.end_)
     or
     t.current_.stopped_ == 1
     or
     MovePossible_() != 0) {
    t.visited_[t.n_visited_] = t.current_;
    t.n_visited_++;
    return 1;
  }

  // 4. Now update alpha in preparation for the next state that is moved to
  UpdateAlpha_(t);
  
  // 5. Draw a Uniform and update the trajectory based upon the cumulative probs.
  double rnd = Uniform(gen);
  if(pi.F_on_) {
    for(unsigned int i=0; i<cforwards_probabilities_.size(); i++)
      if(rnd < cforwards_probabilities_[i]) {
	MoveForwards_(t, i);
	previous_feature_change_ = i;
	previous_direction_ = 1;
	return 0;
      }
  }

  if(pi.B_on_) {
    for(unsigned int i=0; i<cbackwards_probabilities_.size(); i++)
      if(rnd < cbackwards_probabilities_[i]) {
	MoveBackwards_(t, i);
	previous_feature_change_ = i;
	previous_direction_ = -1;
	return 0;
      }
  }

  if(pi.O_on_) {
    for(unsigned int i=0; i<cother_states_probabilities_.size(); i++)
      if(rnd < cother_states_probabilities_[i]) {
	MoveOtherStates_(t, i);
	return 0;
      }
  }
  // If reaches here without returning, then problem and return 1 not 0
  return 1;
}

int Updater::UpdateTrajectory_(Trajectory& t, Transition& pi, std::mt19937& gen, std::uniform_real_distribution<double>& dis)
{
  // Incomplete data currently represented as a 2.
  int incomplete_variable = 2;

  // 1. Calculations for compatibility and weight of edges for possible moves.
  if(previous_feature_change_ == -2) {
    CalculateTotalWeight_(t, pi, 0);
    CalculateTotalWeight_(t, pi, 1);
  }
  else {
   UpdateTotalWeight_(t, pi, 0);
   UpdateTotalWeight_(t, pi, 1);
  } 
  UpdateCProbabilityVectors_(pi);

  // 2. With known transition possibilities, count compatibility of current state.
  UpdateCompatibilityCountAndCumulativeAlpha_(t, incomplete_variable);

  // 3. Can any move be made.
  if(t.current_.IsSame(t.end_)
     or
     t.current_.stopped_ == 1
     or
     MovePossible_() != 0 ) {
    t.visited_[t.n_visited_] = t.current_;
    t.n_visited_++;
    return 1;
  }

  // 4. Now update alpha in preparation for the next state that is moved to
  UpdateAlpha_(t);
  
  // 5. Draw a Uniform and update the trajectory based upon the cumulative probs.
  double rnd = dis(gen);
  if(pi.F_on_) {
    for(unsigned int i=0; i<cforwards_probabilities_.size(); i++)
      if(rnd < cforwards_probabilities_[i]) {
	MoveForwards_(t, i);
	previous_feature_change_ = i;
	previous_direction_ = 1;
	return 0;
      }
  }

  if(pi.B_on_) {
    for(unsigned int i=0; i<cbackwards_probabilities_.size(); i++)
      if(rnd < cbackwards_probabilities_[i]) {
	MoveBackwards_(t, i);
	previous_feature_change_ = i;
	previous_direction_ = -1;
	return 0;
      }
  }

  if(pi.O_on_) {
    for(unsigned int i=0; i<cother_states_probabilities_.size(); i++)
      if(rnd < cother_states_probabilities_[i]) {
	MoveOtherStates_(t, i);
	return 0;
      }
  }

  // If reaches here without returning, then problem and return 1 not 0
  return 1;
}


int Updater::UpdateTrajectory_(Trajectory& t, Transition& pi)
{
  // Incomplete data currently represented as a 2.
  int incomplete_variable = 2;

  // 1. Calculations for compatibility and weight of edges for possible moves.
  if(previous_feature_change_ == -2) {
    CalculateTotalWeight_(t, pi, 0);
    CalculateTotalWeight_(t, pi, 1);
  }
  else {
   UpdateTotalWeight_(t, pi, 0);
   UpdateTotalWeight_(t, pi, 1);
  } 
  UpdateCProbabilityVectors_(pi);

  // 2. With known transition possibilities, count compatibility of current state.
  UpdateCompatibilityCountAndCumulativeAlpha_(t, incomplete_variable);

  // 3. Can any move be made.
  if(t.current_.IsSame(t.end_)
     or
     t.current_.stopped_ == 1
     or
     MovePossible_() != 0 ) {
    t.visited_[t.n_visited_] = t.current_;
    t.n_visited_++;
    return 1;
  }

  // 4. Now update alpha in preparation for the next state that is moved to
  UpdateAlpha_(t);
  
  // 5. Draw a Uniform and update the trajectory based upon the cumulative probs.
  double rnd = Uniform();
  if(pi.F_on_) {
    for(unsigned int i=0; i<cforwards_probabilities_.size(); i++)
      if(rnd < cforwards_probabilities_[i]) {
	MoveForwards_(t, i);
	previous_feature_change_ = i;
	previous_direction_ = 1;
	return 0;
      }
  }

  if(pi.B_on_) {
    for(unsigned int i=0; i<cbackwards_probabilities_.size(); i++)
      if(rnd < cbackwards_probabilities_[i]) {
	MoveBackwards_(t, i);
	previous_feature_change_ = i;
	previous_direction_ = -1;
	return 0;
      }
  }

  if(pi.O_on_) {
    for(unsigned int i=0; i<cother_states_probabilities_.size(); i++)
      if(rnd < cother_states_probabilities_[i]) {
	MoveOtherStates_(t, i);
	return 0;
      }
  }

  // If reaches here without returning, then problem and return 1 not 0
  return 1;
}

// TODO: non-functional method for now as related to incomplete data problem
int Updater::UpdateTrajectoryToCompatible_(Trajectory& t, Transition& pi, int must_move)
{
  // Incomplete data currently represented as a 2.
  int incomplete_variable = 2;

  // 1. Calculations for compatibility and weight of edges for possible moves.
  if(previous_feature_change_ == -2) {
    CalculateTotalWeight_(t, pi, 0);
    CalculateTotalWeight_(t, pi, 1);
  }
  else {
   UpdateTotalWeight_(t, pi, 0);
   UpdateTotalWeight_(t, pi, 1);
  } 
  UpdateCProbabilityVectors_(pi);

  // 2. With known transition possibilities, count compatibility of current state.
  // UpdateCompatibilityCountAndCumulativeAlpha_(t, incomplete_variable);


  // 3. Is current state compatible with target and not already visited
  // A. Compatible
  // B. Must move this round?
  // A & B - return
  int A = t.current_.IsCompatible(t.end_, incomplete_variable) == 1;
  int B = (must_move == 1 ? 0 : 1);
  if((A and B)) {
    // -- NEW
    //t.visited_[t.n_visited_] = t.current_;
    //t.n_visited_++;
    t.cumulative_alpha_ = t.alpha_;
    //std::cout << "1-Alpha: " << t.alpha_ << std::endl;
    return 1;
  }
  
  // 4. Now update alpha in preparation for the next state that is moved to
  UpdateAlpha_(t);
  
  // 5. Draw a Uniform and update the trajectory based upon the cumulative probs.
  double rnd = Uniform();
  if(pi.F_on_) {
    for(unsigned int i=0; i<cforwards_probabilities_.size(); i++)
      if(rnd < cforwards_probabilities_[i]) {
	MoveForwards_(t, i);
	previous_feature_change_ = i;
	previous_direction_ = 1;
	return 0;
      }
  }

  if(pi.B_on_) {
    for(unsigned int i=0; i<cbackwards_probabilities_.size(); i++)
      if(rnd < cbackwards_probabilities_[i]) {
	MoveBackwards_(t, i);
	previous_feature_change_ = i;
	previous_direction_ = -1;
	return 0;
      }
  }

  if(pi.O_on_) {
    for(unsigned int i=0; i<cother_states_probabilities_.size(); i++)
      if(rnd < cother_states_probabilities_[i]) {
	MoveOtherStates_(t, i);
	return 0;
      }
  }

  // If reaches here without returning, then problem and return 1 not 0
  return 1;
}


void Updater::PerformWalkUntilCompatible_(Trajectory& t,
					  Transition& pi,
					  int must_move)
{
  Clear_();
  int stop = 0;
  must_move = 1;
  while(stop == 0) {
    stop = UpdateTrajectoryToCompatible_(t, pi, must_move);
    must_move = 0;
  }
}



int Updater::PerformWalksFromZeroToOne(Trajectory& t,
				       Transition& pi,
				       std::vector <std::vector <double>>& ordering,
				       State& start,
				       State& end)
{
  for(int i=0; i<TESTS_; i++) {
    t.ResetFull_(start, end);
    Clear_();
    int finished = 0;
    while(finished == 0)
      finished = UpdateTrajectory_(t, pi);
      

      
    if(t.n_visited_ != pi.L_+1) {
      std::cout << "n_visited_: " << t.n_visited_ << std::endl;
      std::cout << pi;
      std::cout << t;
      throw std::invalid_argument("Walk not fully progressed.");
    }
    if(t.n_visited_ < 0)
      throw std::invalid_argument("n_visited < 0.");
    for(int j=0; j<t.n_visited_-1; j++) {
      int el = t.visited_[j+1].Difference(t.visited_[j]);
      if(el != -1)
	ordering[j][el] += (double)1/TESTS_;
      else
	throw std::invalid_argument("Not different by one position.");
    }
  }
  return 0;
}

int Updater::PerformWalksFromZeroToOne(Trajectory& t,
				       Transition& pi,
				       std::vector <std::vector <double>>& ordering,
				       State& start,
				       State& end,
				       std::ofstream& f)
{
  for(int i=0; i<TESTS_; i++) {
    t.ResetFull_(start, end);
    Clear_();
    int finished = 0;
    int step_count = 0;
    while(finished == 0) {
      finished = UpdateTrajectory_(t, pi);
      if(!finished)
	f << step_count << "," << previous_feature_change_ << "\n";
      step_count++;
    }

      
    if(t.n_visited_ != pi.L_+1) {
      std::cout << "n_visited_: " << t.n_visited_ << std::endl;
      std::cout << pi;
      std::cout << t;
      throw std::invalid_argument("Walk not fully progressed.");
    }
    if(t.n_visited_ < 0)
      throw std::invalid_argument("n_visited < 0.");
    for(int j=0; j<t.n_visited_-1; j++) {
      int el = t.visited_[j+1].Difference(t.visited_[j]);
      if(el != -1)
	ordering[j][el] += (double)1/TESTS_;
      else
	throw std::invalid_argument("Not different by one position.");
    }
  }
  return 0;
}

void Updater::PerformWalk(Trajectory& t,
			  Transition& pi,
			  State& start,
			  State& end,
			  int start_count,
			  int end_count)			  
{
  t.ResetFull_(start, end);
  Clear_();
  int step_count = start_count;
  int finished = 0;
  while(step_count < end_count) {
    finished = UpdateTrajectory_(t, pi);
    step_count++;
  }
  return;
}

State Updater::PerformWalk(Trajectory& t,
			   Transition& pi,
			   std::vector <std::vector <double>>& ordering,
			   State& start,
			   State& end,
			   std::ofstream& f,
			   int start_count,
			   int end_count,
			   const std::string record)
{
  t.ResetFull_(start, end);
  Clear_();
  int step_count = start_count;
  int finished = 0;
  while(step_count < end_count) {
    finished = UpdateTrajectory_(t, pi);
    if(!finished and record == "yes") {
      f << step_count << "," << previous_feature_change_ << "\n";
      ordering[step_count][previous_feature_change_] += 1.0;
    }
    step_count++;
  }
  return t.current_;
}

// Methods for performing walks over graph structure representation of data
State Updater::PerformWalk(Trajectory& t,
			   Transition& pi,
			   std::vector <std::vector <double>>& ordering,
			   State& start,
			   State& end,
			   std::ofstream& f,
			   int start_count,
			   int end_count,
			   std::map <std::pair<long long, long long>, int>& transition_counter,
			   const std::string record)
{
  std::pair <long long, long long> p = std::make_pair(0, 0);
  t.ResetFull_(start, end);
  Clear_();
  int step_count = start_count;
  int finished = 0;
  long long prev = start.ConvertToLong_();
  while(step_count < end_count) {
    finished = UpdateTrajectory_(t, pi);
    if(!finished and record == "yes") {
      f << step_count << "," << previous_feature_change_ << "\n";
      ordering[step_count][previous_feature_change_] += 1.0;
      p.first = prev;
      p.second = t.current_.ConvertToLong_();
      if(transition_counter.find(p) != transition_counter.end())
	transition_counter[p]++;
      else
	transition_counter[p] = 1;
      prev = p.second;
    }
    step_count++;
  }
  return t.current_;
}

State Updater::PerformWalk(Trajectory& t,
			   Transition& pi,
			   std::vector <std::vector <double>>& ordering,
			   State& start,
			   State& end,
			   std::ofstream& f,
			   int start_count,
			   int end_count,
			   std::map <std::pair<long long, long long>, int>& transition_counter,
			   std::map <long long, int>& state_counter,
			   const std::string record)
{
  std::pair <long long, long long> p = std::make_pair(0, 0);
  t.ResetFull_(start, end);
  Clear_();
  int step_count = start_count;
  int finished = 0;
  long long prev = start.ConvertToLong_();

  while(step_count < end_count) {
    finished = UpdateTrajectory_(t, pi);
    if(!finished and record == "yes") {
      f << step_count << "," << previous_feature_change_ << "\n";
      ordering[step_count][previous_feature_change_] += 1.0;

      p.first = prev;
      p.second = t.current_.ConvertToLong_();


      if(step_count == start_count) {
	if(state_counter.find(p.first) != state_counter.end())
	  state_counter[p.first]++;
	else
	  state_counter[p.first] = 1;
      }
      if(state_counter.find(p.second) != state_counter.end())
	state_counter[p.second]++;
      else
	state_counter[p.second] = 1;
      
      if(transition_counter.find(p) != transition_counter.end())
	transition_counter[p]++;
      else
	transition_counter[p] = 1;
      prev = p.second;
    }
    step_count++;
  }
  return t.current_;
}

void Updater::Clear_()
{
  // Make everything 0.0
  forwards_probabilities_.assign(forwards_probabilities_.size(), 0.0);
  forwards_compatible_probabilities_.assign(forwards_probabilities_.size(), 0.0);
  cforwards_probabilities_.assign(cforwards_probabilities_.size(), 0.0);
  backwards_probabilities_.assign(backwards_probabilities_.size(), 0.0);
  cbackwards_probabilities_.assign(cbackwards_probabilities_.size(), 0.0);
  other_states_probabilities_.assign(other_states_probabilities_.size(), 0.0);
  cother_states_probabilities_.assign(cother_states_probabilities_.size(), 0.0);

  previous_feature_change_ = -2;
  previous_direction_ = -2;

  return;
}

std::ostream& operator<<(std::ostream& out, const Updater& u)
{
  out << "Forwards:" << std::endl;
  for(unsigned int i=0; i<u.cforwards_probabilities_.size(); i++)
    out << i << "\t" << u.cforwards_probabilities_[i] << std::endl;

  out << "Backwards:" << std::endl;
  for(unsigned int i=0; i<u.cbackwards_probabilities_.size(); i++)
    out << i << "\t" << u.cbackwards_probabilities_[i] << std::endl;

  out << "Other states:" << std::endl;
  for(unsigned int i=0; i<u.cother_states_probabilities_.size(); i++)
    out << i << "\t" << u.cother_states_probabilities_[i] << std::endl;
  
  return out;
}
