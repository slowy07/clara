#ifndef STAT_H_
#define STAT_H_

#include <random>

namespace clara {
namespace stat {
extern std::random_device _rd;
extern std::mt19937 _rng;

class NormalDistribution {
std::normal_distribution<> d;

public:
  NormalDistribution(double mean = 0, double sigma = 1) {
    d = std::normal_distribution<>(mean, sigma);
  }
  double sample() {
    return d(_rng);
  }
};

class UniformRealDistribution {
std::uniform_real_distribution<> d;
public:
  UniformRealDistribution(double a = 0, double b = 1) {
    d = std::uniform_real_distribution<>(a, b);
  }
  double sample() {
    return d(_rng);
  }
};
}
}

#endif
