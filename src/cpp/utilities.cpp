#include "utilities.hpp"

std::mt19937 generator;
std::mt19937 generator2;

std::uniform_real_distribution <double> ud(0,1.0);
std::uniform_real_distribution <long double> uld(0,1.0);

void RandomInitialise(uint32_t SEED)
{
  generator .seed(SEED);
  generator2.seed(SEED);
}

void OutputGeneratorState(std::mt19937& g)
{
  g = generator;
}

void OutputGeneratorState2(std::mt19937& g)
{
  g = generator2;
}

void LoadGeneratorState(std::mt19937& g)
{
  generator = g;
}

void CopyGeneratorState(std::mt19937& from_gen, std::mt19937& to_gen)
{
  to_gen = from_gen;
}

double Uniform()
{
  return ud(generator);
}

long double Uniforml()
{
  return uld(generator);
}

double Uniform(std::mt19937& g)
{
  return ud(g);
}

double Uniform2()
{
  return ud(generator2);
}

double UniformRealRange(double lower, double upper)
{
  std::uniform_real_distribution<double> d(lower, upper);
  return d(generator);
}

double Gaussian(double mu, double sigma)
{
  std::normal_distribution<double>d(mu, sigma);
  return d(generator);
}

double LogUniform(double a, double b)
{
  std::uniform_real_distribution<double>d(a, b);
  return log(d(generator));
}

