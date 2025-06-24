#include <ostream>
#include <sstream>

#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

class MockDisplay : public InterfaceDisplay {
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
