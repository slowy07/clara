#ifndef CLASSFUNCTION_RANDOM_DEVICES_H_
#define CLASSFUNCTION_RANDOM_DEVICES_H_

#include <random>

#include "../internal/classFunction/singleton.h"

namespace clara {
/**
 * @class clara::RandomDevices
 * @brief singleton class that manages the source randomness in the library
 * consist of a wrapper around an std::mt19937 mersenne twister random
 * number geberator engine and std::random_device engine
 */
class RandomDevices final : public internal::Singleton<RandomDevices> {
  friend class internal::Singleton<RandomDevices>;
  std::random_device rd_;

 public:
  std::mt19937 rng_;

 private:
  // initialize and seed the random number generator
  RandomDevices() : rd_{}, rng_{rd_()} {}
  // default constructor
  ~RandomDevices() = default;
};

}  // namespace clara

#endif  // !CLASSFUNCTION_RANDOM_DEVICES_H_
