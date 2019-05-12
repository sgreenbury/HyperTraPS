#include "transition.hpp"

Transition::Transition()
{
}

Transition::Transition(const int L, const double sigma)
{
  Init_(L, sigma);
}

Transition::Transition(const int L, const double sigma, const int limit)
{
  Init_(L, sigma, limit);
}

Transition::Transition(const int L, const double sigma, const std::string direction, const std::string other_states, const int n_other_states)
{
  Init_(L, sigma, direction, other_states, n_other_states);
}

Transition::Transition(const Transition& t)
{
  Init_(t);
}

Transition::~Transition()
{
  Clear_();
}

void Transition::Init_(const int L, const double sigma)
{
  L_     = L;
  sigma_ = sigma;
  n_other_states_ = 0;
  Init_();
}

void Transition::Init_(const int L, const double sigma, const int limit)
{
  L_     = L;
  sigma_ = sigma;
  n_other_states_ = 0;
  Init_();
  limit_ = limit;
  limit_u_ = limit_;
  limit_l_ = -limit_;
}

void Transition::Init_(const int L, const double sigma, const std::string direction, const std::string other_states, const int n_other_states)
{
  L_              = L;
  sigma_          = sigma;
  n_other_states_ = n_other_states;

  Init_();
  UpdateFHD_(direction, other_states);
}

void Transition::Init_()
{
  F_on_ = 1;
  B_on_ = 0;
  O_on_ = 0;

  limit_ = 10;
  limit_u_ = limit_;
  limit_l_ = -limit_;

  p_     = 1.;

  F_weight_ = 0.0;
  B_weight_ = 0.0;
  O_weight_ = 0.0;
  
  F_.resize(L_, std::vector <double> (L_, 0.0));
  F_prev_.resize(L_, std::vector <double> (L_, 0.0));
  F_features_.resize(L_, 1);
}

void Transition::UpdateFHD_(const std::string direction, const std::string other_states)
{
  if(direction == "forwards") {
    B_on_ = 0;
  }

  else
    if(direction == "both") {
      B_on_ = 1;
      B_.resize(L_, std::vector <double> (L_, 0.0));
      B_prev_.resize(L_, std::vector <double> (L_, 0.0));
      B_features_.resize(L_, 1);
    }
    else
      throw std::invalid_argument("Received an invalid direction for trajectories.");

  if(other_states == "yes") {
    O_on_ = 1;
    O_.resize(L_+1, std::vector <double> (n_other_states_, 0.0));
    O_prev_.resize(L_+1, std::vector <double> (n_other_states_, 0.0));
    O_features_.resize(n_other_states_, 1);
  }
  else
    if(other_states == "no") {
      O_on_ = 0;
    }
    else
      throw std::invalid_argument("Received an invalid direction for trajectories.");
}

void Transition::Init_(const Transition& t)
{
  F_on_ = t.F_on_;
  B_on_ = t.B_on_;
  O_on_ = t.O_on_;

  // Copy all elements across
  L_              = t.L_;
  sigma_          = t.sigma_;
  n_other_states_ = t.n_other_states_;
  limit_          = t.limit_;
  limit_u_        = t.limit_u_;
  limit_l_        = t.limit_l_;
  p_              = t.p_;
  
  // Copy across forwards moves
  if(t.F_on_) {
    F_on_ = t.F_on_;
    if(F_.size() != t.F_.size())
      F_.resize(t.F_.size(),std::vector <double> (L_));

    for(unsigned int i=0; i<t.F_.size(); i++)
      F_[i] = t.F_[i];

    if(F_prev_.size() != t.F_prev_.size())
      F_prev_.resize(t.F_prev_.size(), std::vector <double> (L_));

    for(unsigned int i=0; i<t.F_prev_.size(); i++)
      F_prev_[i] = t.F_prev_[i];
    
    F_weight_ = t.F_weight_;

    F_features_ = t.F_features_;
  }

  // TODO: Reverse moves to be implemented
  // Copy across backwards moves
  if(t.B_on_) {
    B_on_ = t.B_on_;

    if(B_.size() != t.B_.size())
      B_.resize(t.B_.size(), std::vector <double> (L_));
      
    for(unsigned int i=0; i<t.B_.size(); i++)
      B_[i] = t.B_[i];

    if(B_prev_.size() != t.B_prev_.size())
      B_prev_.resize(t.B_prev_.size(), std::vector <double> (L_));
      
    for(unsigned int i=0; i<t.B_prev_.size(); i++)
      B_prev_[i] = t.B_prev_[i];
    
    B_weight_ = t.B_weight_;

    B_features_ = t.B_features_;
  }

  // TODO: Off hypercube moves to be implemented
  // Copy across other states dependent on the presence or absence of traits
  if(t.O_on_) {
    O_on_ = t.O_on_;

    if(O_.size() != t.O_.size())
      O_.resize(t.O_.size(), std::vector <double> (n_other_states_));
      
    for(unsigned int i=0; i<t.O_.size(); i++) {
      O_[i] = t.O_[i];
    }

    if(O_prev_.size() != t.O_prev_.size())
      O_prev_.resize(t.O_prev_.size(), std::vector <double> (n_other_states_));
    
    for(unsigned int i=0; i<t.O_prev_.size(); i++)
      O_prev_[i] = t.O_prev_[i];

    O_weight_ = t.O_weight_;

    O_features_ = t.O_features_;
  }
}


Transition& Transition::operator=(const Transition &t)
{
  if( &t != this ) {
    Clear_();
    Init_(t);
  }
  return *this;
}

void Transition::Load_(const std::string file_name, const std::string load_type)
{
  // Loads each of the different possible types of matrix
  std::ifstream in_file(file_name);
  std::string line;
  int row = 0;
  int col = 0;
  double d;
  while(std::getline(in_file, line)) {
    col = 0;
    std::istringstream line_ss(line);
    while(line_ss >> d) {
      if((col >= L_ and load_type != "other_states") or (col > 2 and load_type == "other_states"))
	throw std::invalid_argument("Too many columns in data file.");
      if(load_type == "forwards")
	F_[row][col] = d;
      if(load_type == "backwards")
	B_[row][col] = d;
      if(load_type == "other_states")
	O_[row][col] = d;
      col++;
    }
    row++;
  }
  in_file.close();
}

int Transition::FindL_(const std::string file_name, const std::string load_type)
{
  std::ifstream in_file(file_name);
  std::string line;
  std::string line_label;
  std::string row_string;
  int row = 0;
  int col = 0;
  double d;
  while(std::getline(in_file, line)) {
    std::stringstream linestream(line);
    std::getline(linestream, line_label, '\t');
    for(row = 0; row < 1; row++) {
      std::getline(linestream, row_string, '\t');
      std::stringstream rowstream(row_string);
      col = 0;
      while(rowstream >> d) {
	col++;
      }
      L_ = col;
    }
    in_file.close();
    return L_;
  }
  return -1;
}


int Transition::Load_(const std::string file_name, const std::string load_type, const std::string label)
{
  // Loads each of the different possible types of matrix
  std::ifstream in_file(file_name);
  std::string line;
  std::string line_label;
  std::string row_string;
  int row = 0;
  int col = 0;
  double d;
  while(std::getline(in_file, line)) {
    std::stringstream linestream(line);
    std::getline(linestream, line_label, '\t');
    if(line_label != label)
      continue;
    else {
      for(row = 0; row < L_; row++) {
	std::getline(linestream, row_string, '\t');
	std::stringstream rowstream(row_string);
	col = 0;
	while(rowstream >> d) {
	  if((col >= L_ and load_type != "other_states") or (col > 2 and load_type == "other_states"))
	    throw std::invalid_argument("Too many columns in data file.");
	  if(load_type == "forwards")
	    F_[row][col] = d;
	  if(load_type == "backwards")
	    B_[row][col] = d;
	  if(load_type == "other_states")
	    O_[row][col] = d;
	  col++;
	}
      }
      in_file.close();
      return 0;
    }
  }
  in_file.close();
  return 1;
}


void Transition::Write_(const std::string type, std::ofstream& out_file)
{
  // Loads each of the different possible types of matrix
  std::vector <std::vector <double>>* mptr = 0;
  if(type == "forwards")
    mptr = &F_;
  else
    if(type == "backwards")
      mptr = &B_;
    else
      if(type == "other_states")
	mptr = &O_;
    
    
  for(unsigned int i=0; i<(*mptr).size(); i++) {
    for(unsigned int j=0; j<(*mptr).size(); j++) {
      out_file << (*mptr)[i][j];
      if(j < (*mptr).size() - 1)
	out_file << " ";
      else
	if(i < (*mptr).size() - 1)
	  out_file << "\t";
	else
	  out_file << std::endl;
    }
  }
}

int Transition::Write_(const std::string file_name_index)
{
  // Loads each of the different possible types of matrix
  std::ofstream out_file;
  if(F_on_ == 1) {
    out_file.open(file_name_index + "_forwards.txt");
    for(unsigned int i=0; i<F_.size(); i++) {
      for(unsigned int j=0; j<F_[i].size(); j++)
	out_file << F_[i][j] << " ";
      out_file << std::endl;
    }
    out_file.close();
  }
  
  if(B_on_ == 1) {
    out_file.open(file_name_index + "_backwards.txt");
    for(unsigned int i=0; i<B_.size(); i++) {
      for(unsigned int j=0; j<B_[i].size(); j++)
	out_file << B_[i][j] << " ";
      out_file << std::endl;
    }
    out_file.close();
  }
  
  if(O_on_ == 1) {
    out_file.open(file_name_index + "_other-states.txt");
    for(unsigned int i=0; i<O_.size(); i++) {
      for(unsigned int j=0; j<O_[i].size(); j++)
	out_file << O_[i][j] << " ";
      out_file << std::endl;
    }
    out_file.close();
  }
  
  return 0;
}


void Transition::ChangeSigma_(double sigma)
{
  sigma_ = sigma;
}

void Transition::RestorePreviousTransitions_()
{
  if(F_on_)
    for(unsigned int i=0; i<F_.size(); i++)
      for(unsigned int j=0; j<F_[i].size(); j++)
	F_[i][j] = F_prev_[i][j];

  if(B_on_)
    for(unsigned int i=0; i<B_.size(); i++)
      for(unsigned int j=0; j<B_[i].size(); j++)
	B_[i][j] = B_prev_[i][j];

  if(O_on_)
    for(unsigned int i=0; i<O_.size(); i++)
      for(unsigned int j=0; j<O_[i].size(); j++)
	O_[i][j] = O_prev_[i][j];
}

int Transition::MakeBiasedMatrix_(std::vector <std::vector <double> >& m)
{
  // maximum large log number
  double large_log = log(std::numeric_limits<double>::max());
  for(unsigned int i=0; i<m.size(); i++) {
    for(unsigned int j=0; j<m.size(); j++) {
      if(i == 0 and j==0)
	m[i][j] = large_log;
      else
	if(j+1 == i)
	  m[i][j] = large_log;
	else
	  m[i][j] = 0.0;
    }
  }

  return 0;
}

int Transition::MakeTwoWayMatrix_(std::vector <std::vector <double> >& m)
{
  // maximum large log number
  //double large_log = log(std::numeric_limits<double>::max());
  double large_log = 10;
  for(unsigned int i=0; i<m.size(); i++) {
    for(unsigned int j=0; j<m.size(); j++) {
      if((i == 0 and j==0) or (i == m.size()-1 and j == m.size()-1))
	m[i][j] = 0;
      else
	if(i==j)
	  m[i][j] = -large_log;
	else
	  if((j+1 == i) or (j==i+1))
	    m[i][j] = large_log;
	  else
	    m[i][j] = -large_log;
    }
  }

  return 0;
}

int Transition::MakeZeroOrderMatrix_(std::vector <std::vector <double> >& m)
{
  // maximum large log number
  double small_log = log(std::numeric_limits<double>::min());
  for(unsigned int i=0; i<m.size(); i++) {
    for(unsigned int j=0; j<m.size(); j++) {
      if(i == j)
	m[i][j] = 0.0;
      else
	m[i][j] = small_log;
    }
  }

  return 0;
}


int Transition::MakeRandomMatrix_(std::vector <std::vector <double> >& m)
{
  for(unsigned int i=0; i<m.size(); i++) {
    for(unsigned int j=0; j<m.size(); j++) {
      // Arbitrarily chosen width of gaussian
      m[i][j] = Gaussian(0.0, 1.0);
    }
  }

  return 0;
}

int Transition::MakeUniformMatrix_(std::vector <std::vector <double> >& m)
{
  for(unsigned int i=0; i<m.size(); i++) {
    for(unsigned int j=0; j<m.size(); j++) {
      m[i][j] = 0.0;
    }
  }
  return 0;
}

// Redundant: Now called MakeUniform
int Transition::MakeRandomLimitMatrix_(std::vector <std::vector <double> >& m)
{
  for(unsigned int i=0; i<m.size(); i++) {
    for(unsigned int j=0; j<m.size(); j++) {
      m[i][j] = UniformRealRange((double)limit_l_, (double)limit_u_);
    }
  }
  return 0;
}

int Transition::MakeUniform_(std::vector <std::vector <double> >& m)
{
  for(unsigned int i=0; i<m.size(); i++) {
    for(unsigned int j=0; j<m.size(); j++) {
      m[i][j] = UniformRealRange((double)limit_l_, (double)limit_u_);
    }
  }
  return 0;
}

int Transition::MakeUniformBasal_(std::vector <std::vector <double> >& m)
{
  for(unsigned int i=0; i<m.size(); i++) {
    m[i][i] = UniformRealRange((double)limit_l_, (double)limit_u_);
  }
  return 0;
}

int Transition::MakeUniformBasalNoiseInteraction_(std::vector <std::vector <double> >& m)
{
  for(unsigned int i=0; i<m.size(); i++) {
    for(unsigned int j=0; j<m.size(); j++) {
      if(i==j)
	m[i][j] = UniformRealRange((double)limit_l_, (double)limit_u_);
      else
	m[i][j] = Gaussian(0.0, limit_/(double)10);
    }
  }
  return 0;
}

int Transition::MakeNoiseBasal_(std::vector <std::vector <double> >& m)
{
  for(unsigned int i=0; i<m.size(); i++) {
    m[i][i] = Gaussian(0.0, limit_/(double)10);
  }
  return 0;
}

int Transition::MakeBiased_(std::vector <std::vector <double> >& m)
{
  MakeBiasedMatrix_(m);
  return 0;
}

int Transition::MakeNoise_(std::vector <std::vector <double> >& m)
{
  for(unsigned int i=0; i<m.size(); i++) {
    for(unsigned int j=0; j<m.size(); j++) {
      m[i][j] = Gaussian(0.0, limit_/(double)10);
    }
  }
  return 0;
}

int Transition::Perturb_()
{
  // Assume equal types of normal perturbation across all elements that are swithched on
  if(F_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<F_.size(); i++)
      for(unsigned int j=0; j<F_[i].size(); j++)
	if(F_features_[j] and uniform < p_) {
	// ** if(F_features_[j]) {
	  F_prev_[i][j] = F_[i][j];
	  F_[i][j] += Gaussian(0.0, sigma_);
	  // **Check the effect of keeping to with -10,10
	  if(F_[i][j] > limit_u_)  F_[i][j] =  limit_u_;
	  if(F_[i][j] < limit_l_)  F_[i][j] =  limit_l_;
	}
  }
  if(B_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<B_.size(); i++)
      for(unsigned int j=0; j<B_[i].size(); j++)
	if(B_features_[j] and uniform < p_) {
	// ** if(B_features_[j]) {
	  B_prev_[i][j] = B_[i][j];
	  B_[i][j] += Gaussian(0.0, sigma_);
	  if(B_[i][j] > limit_u_)  B_[i][j] =  limit_u_;
	  if(B_[i][j] < limit_l_)  B_[i][j] =  limit_l_;
	}
  }
  if(O_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<O_.size(); i++)
      for(unsigned int j=0; j<O_[i].size(); j++) {
	if(O_features_[j] and uniform < p_) {
	//if(O_features_[j]) {
	  O_prev_[i][j] = O_[i][j];
	  O_[i][j] += Gaussian(0.0, sigma_);
	  if(O_[i][j] > limit_u_)  O_[i][j] =  limit_u_;
	  if(O_[i][j] < limit_l_)  O_[i][j] =  limit_l_;
	}
      }
  }
  return 0;
}

double Transition::PerturbLU_()
{
  // Returns the proposal quotient to be used in MH update
  double proposal_counter = 0.0;
  if(F_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<F_.size(); i++)
      for(unsigned int j=0; j<F_[i].size(); j++)
	if(F_features_[j] and uniform < p_) {
	// ** if(F_features_[j]) {
	  F_prev_[i][j] = F_[i][j];
	  double draw = LogUniform((double)1/sigma_, sigma_);
	  proposal_counter += -2*draw;
	  F_[i][j] += draw;
	  // **Check the effect of keeping to with -10,10
	  if(F_[i][j] > limit_u_)  F_[i][j] =  limit_u_;
	  if(F_[i][j] < limit_l_)  F_[i][j] =  limit_l_;
	}
  }
  return proposal_counter;
  /* Non-functional currently */
  /*
  if(B_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<B_.size(); i++)
      for(unsigned int j=0; j<B_[i].size(); j++)
	if(B_features_[j] and uniform < p_) {
	// ** if(B_features_[j]) {
	  B_prev_[i][j] = B_[i][j];
	  B_[i][j] += LogUniform((double)1/sigma_, sigma_);
	  if(B_[i][j] > limit_)  B_[i][j] =  limit_;
	  if(B_[i][j] < -limit_) B_[i][j] = -limit_;
	}
  }
  if(O_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<O_.size(); i++)
      for(unsigned int j=0; j<O_[i].size(); j++) {
	if(O_features_[j] and uniform < p_) {
	//if(O_features_[j]) {
	  O_prev_[i][j] = O_[i][j];
	  O_[i][j] += LogUniform((double)1/sigma_, sigma_);
	  if(O_[i][j] > limit_)  O_[i][j] =  limit_;
	  if(O_[i][j] < -limit_) O_[i][j] = -limit_;
	}
      }
  }
  return 0;
  */
}

double Transition::PerturbU_()
{
  // Returns the proposal quotient to be used in MH update
  double proposal_counter = 0.0;
  if(F_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<F_.size(); i++)
      for(unsigned int j=0; j<F_[i].size(); j++)
	if(F_features_[j] and uniform < p_) {
	// ** if(F_features_[j]) {
	  F_prev_[i][j] = F_[i][j];
	  double draw = UniformRealRange(-sigma_, sigma_);
	  //proposal_counter += -2*draw;
	  F_[i][j] += draw;
	  // **Check the effect of keeping to with -10,10
	  if(F_[i][j] > limit_u_)  F_[i][j] =  limit_u_;
	  if(F_[i][j] < limit_l_)  F_[i][j] =  limit_l_;
	}
  }
  return proposal_counter;
  /* Non-functional currently */
  /*
  if(B_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<B_.size(); i++)
      for(unsigned int j=0; j<B_[i].size(); j++)
	if(B_features_[j] and uniform < p_) {
	// ** if(B_features_[j]) {
	  B_prev_[i][j] = B_[i][j];
	  B_[i][j] += LogUniform((double)1/sigma_, sigma_);
	  if(B_[i][j] > limit_)  B_[i][j] =  limit_;
	  if(B_[i][j] < -limit_) B_[i][j] = -limit_;
	}
  }
  if(O_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<O_.size(); i++)
      for(unsigned int j=0; j<O_[i].size(); j++) {
	if(O_features_[j] and uniform < p_) {
	//if(O_features_[j]) {
	  O_prev_[i][j] = O_[i][j];
	  O_[i][j] += LogUniform((double)1/sigma_, sigma_);
	  if(O_[i][j] > limit_)  O_[i][j] =  limit_;
	  if(O_[i][j] < -limit_) O_[i][j] = -limit_;
	}
      }
  }
  return 0;
  */
}


int Transition::PerturbDiagonal_()
{
  // Assume equal types of normal perturbation across all elements that are swithched on
  if(F_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<F_.size(); i++)
      if(F_features_[i] and uniform < p_) {
	// ** if(F_features_[j]) {
	F_prev_[i][i] = F_[i][i];
	F_[i][i] += Gaussian(0.0, sigma_);
	// **Check the effect of keeping to with -10,10
	if(F_[i][i] > limit_u_)  F_[i][i] =  limit_u_;
	if(F_[i][i] < limit_l_)  F_[i][i] =  limit_l_;
      }
  }
  if(B_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<B_.size(); i++)
      if(B_features_[i] and uniform < p_) {
	// ** if(B_features_[j]) {
	B_prev_[i][i] = B_[i][i];
	B_[i][i] += Gaussian(0.0, sigma_);
	if(B_[i][i] > limit_u_)  B_[i][i] =  limit_u_;
	if(B_[i][i] < limit_l_)  B_[i][i] =  limit_l_;
      }
  }
  if(O_on_) {
    // ** temp added
    double uniform = Uniform();
    for(unsigned int i=0; i<O_.size(); i++)
      if(O_features_[i] and uniform < p_) {
	//if(O_features_[j]) {
	O_prev_[i][i] = O_[i][i];
	O_[i][i] += Gaussian(0.0, sigma_);
	if(O_[i][i] > limit_u_)  O_[i][i] =  limit_u_;
	if(O_[i][i] < limit_l_)  O_[i][i] =  limit_l_;
      }
  }
  return 0;
}

std::ostream& operator<< (std::ostream& out, Transition const& t)
{
  out << "Forward transition matrix:\n";
  if(t.F_on_) {
    for(unsigned int i=0; i<t.F_.size(); i++) {
      for(unsigned int j=0; j<t.F_[i].size(); j++) {
	out << t.F_[i][j] << "\t";
      }
      out << std::endl;
    }
    out << "Limit upper: " << t.limit_u_ << "\n";
    out << "Limit lower: " << t.limit_l_ << std::endl;
  }

  if(t.B_on_) {
    out << "\nBackward transition matrix:\n";
    for(unsigned int i=0; i<t.B_.size(); i++) {
      for(unsigned int j=0; j<t.B_[i].size(); j++) {
	out << t.B_[i][j] << "\t";
      }
      out << std::endl;
    }
  }

  if(t.O_on_) {
    out << "\nOther transition matrix:\n";
    for(unsigned int i=0; i<t.O_.size(); i++) {
      for(unsigned int j=0; j<t.O_[i].size(); j++) {
	out << t.O_[i][j] << "\t";
      }
      out << std::endl;
    }
  }
  out << std::endl;
  return out;
}

void Transition::Clear_()
{
  F_.clear();
  F_prev_.clear();
  B_features_.clear();
  
  B_.clear();
  B_prev_.clear();
  B_features_.clear();

  O_.clear();
  O_prev_.clear();
  O_features_.clear();
}
