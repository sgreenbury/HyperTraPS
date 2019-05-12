#ifndef __TRANSITION__HPP
#define __TRANSITION__HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <limits>
#include "utilities.hpp"


class Transition
{
public:
  std::vector <std::vector <double> > F_;
  std::vector <std::vector <double> > B_;
  std::vector <std::vector <double> > O_;

  std::vector <std::vector <double> > F_prev_;
  std::vector <std::vector <double> > B_prev_;
  std::vector <std::vector <double> > O_prev_;
  
  std::vector <int> F_features_;
  std::vector <int> B_features_;
  std::vector <int> O_features_;

  int F_on_;
  int B_on_;
  int O_on_;
  
  double F_weight_;
  double B_weight_;
  double O_weight_;

  int L_;
  double sigma_;
  int n_other_states_;

  int limit_;
  double limit_u_;
  double limit_l_;
  double p_;
  
  Transition();
  Transition(const int L, const double sigma);
  Transition(const int L, const double sigma, const int limit);
  Transition(const int L, const double sigma, const std::string direction, const std::string other_states, const int n_other_states);
  ~Transition();
  Transition(const Transition& t);
  Transition& operator=(const Transition& t);

  //private:                                                                                       
  void Init_();
  void Init_(const int L, const double sigma);
  void Init_(const int L, const double sigma, const int limit);
  void Init_(const int L, const double sigma, const std::string direction, const std::string other_states, const int n_other_states);
  void Init_(const Transition& t);

  int FindL_(const std::string file_name, const std::string load_type);


  void UpdateFHD_(const std::string direction, const std::string other_states);

  void RestorePreviousTransitions_();
  
  void ChangeSigma_(double sigma);
  void Clear_();
  
  void Load_(const std::string file_name, const std::string load_type);
  int  Load_(const std::string file_name, const std::string load_type, const std::string label);
  
  //  int LoadInline_(const std::string file_name, const std::string load_type);
  void Write_(const std::string type, std::ofstream& out_file);
  int Write_(const std::string file_name);

  int MakeBiasedMatrix_(std::vector <std::vector <double> >& m);
  int MakeTwoWayMatrix_(std::vector <std::vector <double> >& m);
  int MakeRandomMatrix_(std::vector <std::vector <double> >& m);
  int MakeRandomLimitMatrix_(std::vector <std::vector <double> >& m);
  int MakeUniformMatrix_(std::vector <std::vector <double> >& m);
  int MakeZeroOrderMatrix_(std::vector <std::vector <double> >& m);

  int MakeUniform_(std::vector <std::vector <double> >& m);
  int MakeUniformBasal_(std::vector <std::vector <double> >& m);
  int MakeUniformBasalNoiseInteraction_(std::vector <std::vector <double> >& m);
  int MakeNoiseBasal_(std::vector <std::vector <double> >& m);
  int MakeBiased_(std::vector <std::vector <double> >& m);
  int MakeNoise_(std::vector <std::vector <double> >& m);
  
  int Perturb_();
  int PerturbDiagonal_();
  
  double PerturbU_();
  double PerturbLU_();
  
  friend std::ostream& operator<<(std::ostream& out, const Transition& t);  
};

// For showing current state using in streams
std::ostream& operator<<(std::ostream& out, const Transition& t);


#endif
