#ifndef __UTILITIES__HPP
#define __UTILITIES__HPP

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iomanip>
#include <random>
#include <fstream>
#include <sstream>
#include <limits>
#include <cassert>
#include <algorithm>

void RandomInitialise(uint32_t SEED);
void OutputGeneratorState(std::mt19937& g);
void OutputGeneratorState2(std::mt19937& g);

void LoadGeneratorState(std::mt19937& g);
void CopyGeneratorState(std::mt19937& from_gen, std::mt19937& to_gen);

double Uniform();
long double Uniforml();
double Uniform(std::mt19937& g);
double Uniform2();

double UniformRealRange(double lower, double upper);
double Gaussian(double mu, double sigma);
double LogUniform(double a, double b);

template <typename T>
inline bool DoubleEqual(T a, T b, T error_factor=1.0)
{
  return a==b || 
    std::abs(a-b)<std::abs(std::min(a,b))*std::numeric_limits<T>::epsilon()*error_factor;
}

template <typename T>
std::string ToStringWithPrecision(const T a_value, const int n = 3)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

template <typename T>
long double Mean(std::vector <T>& v)
{
  return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

template <typename T>
long double MeanExp(std::vector <T>& v)
{
  long double total = 0.0;
  for(unsigned int i=0; i<v.size(); i++) {
    total += expl(v[i]);
  }
  return total/v.size();
}

template <typename T>
long double SDunbiased(std::vector <T>& v)
{
  long double mean = Mean(v);
  long double n = v.size();
  long double count = 0.0;
  for(unsigned int i=0; i<v.size(); i++)
    count += (v[i]-mean)*(v[i]-mean);
  return sqrtl(1/(n-1) * count);
}


template <typename T>
int WriteMatrixCSV(std::vector <std::vector<T>>& matrix, std::string file_name)
{
  std::ofstream f;
  f.open(file_name);
  for(unsigned int i=0; i<matrix.size(); i++) {
    for(unsigned int j=0; j<matrix[i].size(); j++) {
      if(j < matrix[i].size()-1)
	f << matrix[i][j] << ",";
      else
	f << matrix[i][j] << "\n";
    }
  }
  f.close();
  return 0;
}

template <typename T>
int LoadData(const std::string file_name, std::vector <std::vector <T> >& v)
{
  std::ifstream in_file(file_name);
  std::string line;
  int row = 0;
  int col = 0;
  double d;
  while(std::getline(in_file, line)) {
    std::vector <T> line_d;
    std::istringstream line_ss(line);
    while(line_ss >> d) {
      line_d.push_back((T)d);
      col++;
    }
    v.push_back(line_d);
    row++;
  }
  in_file.close();
  return 0;
}

#endif
