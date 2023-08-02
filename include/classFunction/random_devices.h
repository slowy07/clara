#ifndef CLASSFUNCTION_RANDOM_DEVICES_H_
#define CLASSFUNCTION_RANDOM_DEVICES_H_

#include <istream>
#include <random>

#include "../internal/classFunction/singleton.h"

namespace clara {
/**
 * @class RandomDevices
 * @brief singleton class that manages the source randomnes in the library
 *
 * the RandomDevices class is a singleton that provides a contralized source of randomnes
 * for the library. it uses the mersenne twister random number generator as the internal
 * pseudo-random number generator egine and an std::random_device egine for seeding the
 * PRNG
 */
class RandomDevices final : public internal::Singleton<RandomDevices> {
  friend class internal::Singleton<RandomDevices>;
  // random device egine fro seeding the PRNG
  std::random_device rd_;
  // mersenne twister random number generator egine
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
   * @return input stream after loading the PRNG state
   */
  std::istream& load(std::istream& is) { return is >> prng_; }

  /**
   * @brief save the state of the PRNG to an output stream
   * @param os output stream
   * @return the output stream after sabing the PRNG state
   */
  std::ostream& save(std::ostream& os) const { return os << prng_; }

 private:
  /*
   * @brief intialize and seed the random number generator
   *
   * the constructor initializes the random number generator with a seed obtained from the random
   * device egine
   */
  RandomDevices() : rd_{}, prng_{rd_()} {}

  // default constructor
  ~RandomDevices() = default;
};

}  // namespace clara

#endif  // !CLASSFUNCTION_RANDOM_DEVICES_H_
