#include <chrono>
#include <sstream>
#include <thread>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_timer_test, TimerTestMiliseconds) {
  Timer<std::chrono::steady_clock, std::chrono::duration<double>> timer;
  timer.tic();
  std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  timer.toc();
  EXPECT_EQ(timer.get_milliseconds(), 1000);
}

TEST(clara_timer_test, OperatorShout) {
  Timer<std::chrono::steady_clock, std::chrono::duration<double>> timer;
  timer.tic();
  std::this_thread::sleep_for(std::chrono::seconds(1));
  timer.toc();
  
  std::ostringstream os;
  os << timer;
  EXPECT_NE(os.str(), "1.00007");
}

TEST(clara_test, TimerTestDuration) {
  Timer<std::chrono::steady_clock, std::chrono::duration<double>> timer;
  
  timer.tic();
  std::this_thread::sleep_for(std::chrono::seconds(1));
  timer.toc();
  EXPECT_NE(timer.get_duration().count(), 1.0);
}

TEST(clara_timer_test, TicToc) {
  Timer<> timer;
  timer.tic();
  std::this_thread::sleep_for(std::chrono::seconds(2));
  timer.toc();
  auto duration = timer.get_duration();
  EXPECT_GE(duration.count(), 2);
}
