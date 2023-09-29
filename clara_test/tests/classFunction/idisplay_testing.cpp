#include <ostream>
#include <sstream>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

class MockDisplay : public IDisplay {
 public:
  virtual std::ostream& display(std::ostream& os) const override {
    os << "MockTestDisplay";
    return os;
  }
};

TEST(clara_display_test, DisplayMockObject) {
  MockDisplay mock;
  std::ostringstream oss;
  oss << mock;
  EXPECT_EQ(oss.str(), "MockTestDisplay");
}

class OtherTestMockDisplay : public IDisplay {
 public:
  virtual std::ostream& display(std::ostream& os) const override {
    os << "OtherTestMockDisplay";
    return os;
  }
};

TEST(clara_display_test, OtherDisplayMockObject) {
  OtherTestMockDisplay anotherMock;
  std::ostringstream oss;
  oss << anotherMock;

  EXPECT_EQ(oss.str(), "OtherTestMockDisplay");
}
