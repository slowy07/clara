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

TEST(clara_reversibe_test, SetAndGet) {
  Dynamic_bitset bitset(8);
  bitset.set(3);
  EXPECT_TRUE(bitset.get(3));
  EXPECT_FALSE(bitset.get(2));
}

TEST(clara_reversibe_test, InitializationCircuit) {
  Dynamic_bitset bitset(8);
  Bit_circuit circuit(bitset);
  EXPECT_EQ(circuit.size(), 8);
}

TEST(clara_reversibe_test, InitializationCircuiTenBitset) {
  Dynamic_bitset bitset(10);
  EXPECT_EQ(bitset.size(), 10);
  EXPECT_EQ(bitset.storage_size(), 1);
}

TEST(clara_reversibe_test, XGate) {
  Dynamic_bitset bitset(8);
  Bit_circuit circuit(bitset);
  circuit.X(3);
  EXPECT_TRUE(circuit.get(3));
}
