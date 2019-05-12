#include "state.hpp"

State::State()
{
}

State::State(const State& s)
{
  Init_(s);
}

State::State(const std::vector <int>& position, const int stopped)
{
  Init_(position, stopped);
}

State::~State()
{
  Clear_();
}


void State::Init_(const std::vector <int>& position, const int stopped)
{
  position_ = position;
  stopped_  = stopped;
}

void State::Init_(const State& s)
{
  position_= s.position_;
  stopped_ = s.stopped_;
}

void State::Clear_()
{
  return;
}

State& State::operator=(const State &s)
{
  if( &s != this ) {
    Clear_();
    Init_(s);
  }
  return *this;
}

int State::IsSame(const State& s)
{
  if(position_.size() != s.position_.size())
    return 0;

  for(unsigned int i=0; i<position_.size(); i++)
    if(position_[i] != s.position_[i])
      return 0;

  if(stopped_ != s.stopped_)
    return 0;

  return 1;
}

int State::IsCompatible(const State& s, const int incomplete_variable)
{
  if(position_.size() != s.position_.size())
    return 0;

  for(unsigned int i=0; i<position_.size(); i++)
    if(position_[i] != s.position_[i]
       and
       position_[i] != incomplete_variable
       and
       s.position_[i] != incomplete_variable)
      return 0;

  if(stopped_ != s.stopped_)
    return 0;

  return 1;
}

int State::IsCompatible(const std::vector <int>& s, const int incomplete_variable)
{
  if(position_.size() != s.size())
    return 0;

  for(unsigned int i=0; i<position_.size(); i++)
    if(position_[i] != s[i]
       and
       position_[i] != incomplete_variable
       and
       s[i] != incomplete_variable)
      return 0;

  return 1;
}


int State::MovePossible(const State& s)
{
  if(position_.size() != s.position_.size())
    return 0;

  for(unsigned int i=0; i<position_.size(); i++)
    if(position_[i] == 0
       and
       s.position_[i] != 0)
      return 1;

  if(stopped_ != 0)
    return 0;

  return 0;
}

int State::MovePossible(const std::vector <int>& s)
{
  if(position_.size() != s.size())
    return 0;

  for(unsigned int i=0; i<position_.size(); i++)
    if(position_[i] == 0
       and
       s[i] != 0)
      return 1;
  
  if(stopped_ != 0)
    return 0;

  return 0;
}



int State::Difference(const State& s)
{
  if(position_.size() != s.position_.size()
     or
     stopped_ != s.stopped_
     )
    return -1;

  int el=-1;
  int count = 0;
  for(unsigned int i=0; i<position_.size(); i++) {
    if(position_[i] != s.position_[i]) {
      count++;
      el = i;
    }
  }
  if(count == 1)
    return el;

  return -1;
}

long long State::ConvertToLong_()
{
  long long acc = 0;
  for(unsigned int i=0; i<position_.size(); i++)
    acc += (position_[i] == 1 ? (long)pow(2,i) : 0);
  return acc;
}

int State::CountOnes_()
{
  int acc = 0;
  for(unsigned int i=0; i<position_.size(); i++)
    acc += (position_[i] == 1 ? 1 : 0);
  return acc;
}

std::ostream& operator<< (std::ostream& out, State const& s)
{
  for(unsigned int i=0; i<s.position_.size(); i++) {
    out << s.position_[i];
  }
  out << "\t";
  out << s.stopped_;
  out << std::endl;
  return out;
}
