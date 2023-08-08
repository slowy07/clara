#include <complex>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>

#include "../../../include/clara.h"
#include "../../../include/classFunction/random_devices.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_random_device_test, GetPrng) {
  std::mt19937& prng = RandomDevices::get_instance().get_prng();
  EXPECT_NE(prng.state_size, std::mt19937::default_seed);

  int random_number = prng();
  EXPECT_NE(random_number, 0);
}

TEST(clara_random_device_test, loadGetPrng) {
  std::stringstream ss;
  clara::rdevs.save(ss);
  cmat A1 = rand<cmat>(4, 4);
  ss.seekg(0);
  clara::rdevs.load(ss);
  cmat A2 = rand<cmat>(4, 4);
  EXPECT_EQ(0, norm(A1 - A2));
}

TEST(clara_random_device_test, loadGetPrngTest) {
  std::stringstream ss;
  clara::rdevs.save(ss);
  bigint b1 = rand(static_cast<bigint>(-100), 100);
  ss.seekg(0);
  clara::rdevs.load(ss);
  bigint b2 = rand(static_cast<bigint>(-100), 100);
  EXPECT_EQ(b1, b2);
}
