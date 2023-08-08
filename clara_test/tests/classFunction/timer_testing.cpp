#include <chrono>
#include <thread>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_timer_test, DefaultConstructor) {
  Timer<> timer;
  EXPECT_GE(timer.tics(), 0);
}

TEST(clara_timer_test, TicToc) {
  Timer<> timer;
  timer.tic();
  std::this_thread::sleep_for(std::chrono::seconds(2));
  timer.toc();
  auto duration = timer.get_duration();
  EXPECT_GE(duration.count(), 2);
}
