#include <chrono>
#include <thread>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;
TEST(clara_Timer_get_duration, AllTests) {
  using namespace std::chrono;

  Timer<> t1;
  std::this_thread::sleep_for(seconds(1));
  t1.toc();

  auto duration_t1_s = t1.get_duration();
  EXPECT_NEAR(duration_t1_s.count(), 1, 0.05);

  auto duration_t1_ms = t1.get_duration<milliseconds>();
  EXPECT_NEAR(duration_t1_ms.count(), 1000, 50);

  Timer<std::chrono::microseconds> t2;
  std::this_thread::sleep_for(microseconds(100000));
  t2.toc();

  auto duration_t2_micros = t2.get_duration();
  EXPECT_NEAR(duration_t2_micros.count(), 100000, 50000);
}

TEST(clara_timer, AllTests) {
  using namespace std::chrono;

  Timer<> t;
  std::this_thread::sleep_for(milliseconds(100));
  t.tic();
  t.toc();
  EXPECT_NEAR(t.tics(), 0, 0.05);
}
