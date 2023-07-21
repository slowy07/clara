#ifndef CLASSFUNCTION_RANDOM_DEVICES_H_
#define CLASSFUNCTION_RANDOM_DEVICES_H_

#include <istream>
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
  std::mt19937 prng_;

 public:
  /**
   * @brief returns a reference to the internal PRNG object
   * @return the reference to the internal PRNG object
   */
  std::mt19937& get_prng() { return prng_; }

  /**
   * @brief load the state of the PRNG from an input stream
   * @param input stream
   * @return input stream
   */
  std::istream& load(std::istream& is) { return is >> prng_; }

  /**
   * @brief save the state of the PRNG to an output stream
   * @param os output stream
   * @return the output stream
   */
  std::ostream& save(std::ostream& os) const { return os << prng_; }

 private:
  /*
   * @brief intialize and seed the random number generator
   */
  RandomDevices() : rd_{}, prng_{rd_()} {}

  // default constructor
  ~RandomDevices() = default;
};

}  // namespace clara

#endif  // !CLASSFUNCTION_RANDOM_DEVICES_H_
