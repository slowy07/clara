#include <chrono>
#include <thread>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_timer, AllTests) {
  using namespace std::chrono;

  Timer<> t;
  std::this_thread::sleep_for(milliseconds(100));
  t.tic();
  t.toc();
  EXPECT_NEAR(t.tics(), 0, 0.05);
}
