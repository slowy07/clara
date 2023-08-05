#include <cmath>
#include "../../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(CodeTest, CodeWordFiveQubit) {
  ket codeword = Codes::get_instance().codeword(Codes::Type::FIVE_QUBIT, 0);
  EXPECT_NEAR(abs(codeword.norm()), 1.0, 1e-12);
}

TEST(CodeTest, Codeword7Qubit) {
  ket codeword = Codes::get_instance().codeword(Codes::Type::SEVEN_QUBIT_STEANE, 0);
  EXPECT_NEAR(abs(codeword.norm()), 1.0, 1e-12);
}

TEST(CodeTest, Codeword9Qubit) {
  ket codeword = Codes::get_instance().codeword(Codes::Type::NINE_QUBIT_SHOR, 0);
  EXPECT_NEAR(abs(codeword.norm()), 1.0, 1e-12);
}
