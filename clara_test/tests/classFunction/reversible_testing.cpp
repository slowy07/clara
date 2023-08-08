#include "../../../include/classFunction/reversible.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_reversibe_test, X) {
  Bit_circuit circuit(4);
  circuit.X(0);
  EXPECT_EQ(circuit.get(0), true);
}

TEST(clara_reversibe_test, CNOT) {
  Bit_circuit circuit(4);
  circuit.X(0);
  circuit.CNOT({0, 1});
  EXPECT_EQ(circuit.get(1), true);
}

TEST(clara_reversibe_test, TOF) {
  Bit_circuit circuit(4);
  circuit.X(0);
  circuit.X(1);
  circuit.TOF({0, 1, 2});
  EXPECT_EQ(circuit.get(2), true);
}
