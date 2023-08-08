#include <chrono>
#include <thread>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

const double TOLERANCE = 1e-6;

class TimerTest : public ::testing::Test {
protected:
  void sleep_for(int milliseconds) {
std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
  }
};

TEST_F(TimerTest, BasicFunctionallyTest) {
  Timer<std::chrono::milliseconds> timer;
  sleep_for(100);
  timer.toc();
  EXPECT_NEAR(timer.tics(), 100.0, TOLERANCE);
}

TEST(clara_timer, AllTests) {
  using namespace std::chrono;

  Timer<> t;
  std::this_thread::sleep_for(milliseconds(100));
  t.tic();
  t.toc();
  EXPECT_NEAR(t.tics(), 0, 0.05);
}
