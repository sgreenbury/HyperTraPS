#ifndef __STATE__HPP
#define __STATE__HPP

#include <iostream>
#include <fstream>
#include <vector>
#include "utilities.hpp"

class State
{
public:
  std::vector <int> position_;
  int stopped_;
   
  State();
  ~State();
  State(const std::vector <int>& position, const int stopped);
  
  State(const State& s);
  State& operator=(const State& s);

  int IsSame(const State& s);
  int IsCompatible(const State& s, const int incomplete_variable);
  int IsCompatible(const std::vector <int>& s, const int incomplete_variable);
  int Difference(const State& s);
  int MovePossible(const State& s);
  int MovePossible(const std::vector <int>& s);
  
  void Init_();
  void Init_(const State& s);
  void Init_(const std::vector <int>& position, const int stopped);

  void Clear_();

  long long ConvertToLong_();
  int CountOnes_();
  
  friend std::ostream& operator<<(std::ostream& out, const State& s);  
};

std::ostream& operator<<(std::ostream& out, const State& s);

#endif
