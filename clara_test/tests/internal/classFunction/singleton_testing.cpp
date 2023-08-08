#include "../../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

class ClaraSingletonTesting : public ::testing::Test {
 protected:
  class TestingSingleton : public clara::internal::Singleton<TestingSingleton> {
   public:
    int getValue() const { return value_; }

    void setValue(int value) { value_ = value; }

   private:
    int value_ = 0;
  };
};

TEST_F(ClaraSingletonTesting, SingletonInstance) {
  TestingSingleton& instance1 = TestingSingleton::get_instance();
  instance1.setValue(42);
  TestingSingleton& instance2 = TestingSingleton::get_instance();
  EXPECT_EQ(instance1.getValue(), instance2.getValue());
}

#ifndef NO_THREAD_LOCAL
TEST_F(ClaraSingletonTesting, ThreadLocalSingletonInstance) {
  TestingSingleton& instance1 = TestingSingleton::get_thread_local_instance();
  instance1.setValue(50);
  TestingSingleton& instance2 = TestingSingleton::get_thread_local_instance();
  EXPECT_EQ(instance1.getValue(), instance2.getValue());
}
#endif  // !NO_THREAD_LOCAL
