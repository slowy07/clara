#ifndef CLASSFUNCTION_REVESIBLE_H_
#define CLASSFUNCTION_REVESIBLE_H_

#include <strings.h>

#include <bitset>
#include <cassert>
#include <climits>
#include <locale>
#include <ostream>
#include <random>
#include <string>
#include <vector>

#include "../types.h"
#include "idisplay.h"

namespace clara {

/**
 * @brief a dynamic bitset calss that implements a flexible-size arrays of bits
 * @details this class represent a dynamic bitset, which is an array of bits
 * that can dynamically grow as needed. it is used for efficient storage and manipulation
 * of binary data
 */
class Dynamic_bitset : public IDisplay {
 public:
  // type used to store individual bits
  using value_type = unsigned int;
  // type for the storage container
  using storage_type = std::vector<value_type>;

 protected:
  // number of storage elements needed
  idx storage_size_;
  // total number of bits
  idx N_;
  // storage vector for the bitset
  std::vector<value_type> v_;

  /**
   * @brief calculate the index of the storage element containing a bit at the given position
   * @param pos the position of the bit
   * @return the index of the storage element
   */
  idx index_(idx pos) const { return pos / (sizeof(value_type) * CHAR_BIT); }
  /**
   * @brief calculate the offset within the storage element fora bit at the given position
   * @param pos the position of the bit
   * @return the offset within the storage element
   */
  idx offset_(idx pos) const { return pos % (sizeof(value_type) * CHAR_BIT); }

 public:
  /**
   * @brief constructor to create a dynamic bitset of a given size
   * @param N the total number of bits in the bitset
   */
  Dynamic_bitset(idx N)
      : storage_size_{N / (sizeof(value_type) * CHAR_BIT) + 1}, N_{N}, v_(storage_size_) {}

  /**
   * @brief get a reference to the internal stroage internal container
   * @return a const reference to the internal storage vector
   */
  const storage_type& data() const { return v_; }
  /**
   * @brief get a reference to the internal storage container
   * @return const reference to the internal storage vector
   */
  idx size() const noexcept { return N_; }
  /**
   * @brief get the number of storage elements used in the bitset
   * @return the number of storage elements
   */
  idx storage_size() const noexcept { return storage_size_; }

  /**
   * @brief count the number of set bits (ones) in the bitset
   * @return the count of set bits in the bitset
   */
  idx count() const noexcept {
    idx result = 0;
    idx bitset_size = this->size();
    for (idx i = 0; i < bitset_size; ++i) {
      if (this->get(i))
        ++result;
    }
    return result;
  }

  /**
   * @brief get the value of the bit at the specified position
   * @param pos position of the bit to retrieve
   * @return the value of the bit at the given position
   */
  bool get(idx pos) const noexcept { return 1 & (v_[index_(pos)] >> offset_(pos)); }

  /**
   * @brief check if not bits are set in the bitset
   * @return true if not bits are set, false otherwise
   */
  bool none() const noexcept {
    bool result = true;
    idx bitset_storage_size = this->storage_size();
    for (idx i = 0; i < bitset_storage_size; ++i) {
      if (v_[i]) {
        return false;
      }
    }
    return result;
  }

  /**
   * @brief check if all bits are set in the bitset
   * @return true if all bits are set, false otherwise
   */
  bool all() const noexcept {
    bool result = true;
    idx bitset_storage_size = this->storage_size();
    for (idx i = 0; i < bitset_storage_size; ++i) {
      if (~v_[i]) {
        return false;
      }
    }
    return result;
  }

  /**
   * @brief check if any bits are set in the bitset
   * @return true if any bit is set, false if no bits are set
   */
  bool any() const noexcept { return !(this->none()); }

  /**
   * @brief set the value of the bit at the specified position
   * @param pos the position of the bit to set
   * @param value the value to set the bit
   * @return reference to the modified Dynamic_bitset
   */
  Dynamic_bitset& set(idx pos, bool value = true) {
    value ? v_[index_(pos)] |= (1 << offset_(pos)) : v_[index_(pos)] &= ~(1 << offset_(pos));
    return *this;
  }

  /**
   * @brief set all bits in the bitset to the specified value
   * @param value the value to set all bits (default is true)
   * @return reference to the modified Dynamic_bitset
   */
  Dynamic_bitset& set() noexcept {
    idx bitset_storage_size = this->storage_size();
    for (idx i = 0; i < bitset_storage_size; ++i) {
      v_[i] = ~0;
    }
    return *this;
  }

  /**
   * @brief randomly set the value of the bit at the specified position
   * @param pos the position of the bit to set randomly
   * @param p the probability of setting the bit of true (default is 0.5)
   * @return reference to the modified Dynamic_bitset
   */
  Dynamic_bitset& rand(idx pos, double p = 0.5) {
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::bernoulli_distribution d{p};

    this->set(pos, d(gen));
    return *this;
  }

  /**
   * @brief randomly set the value of all bits in the bitset
   * @param p the probability of setting each bit to true (default is 0.5)
   * @return reference to the modified Dynamic_bitset
   */
  Dynamic_bitset& rand(double p = 0.5) {
    idx bitset_size = this->size();
    for (idx i = 0; i < bitset_size; ++i) {
      this->rand(i, p);
    }
    return *this;
  }

  /**
   * @brief reset the value of the bit at the specified position to false
   * @param pos the position of the bit to reset
   * @return A reference to the modified Dynamic_bitset
   */
  Dynamic_bitset& reset(idx pos) {
    v_[index_(pos)] &= ~(1 << offset_(pos));
    return *this;
  }

  /**
   * @brief reset the value of all bits in the bitset to false
   * @return reference to the modified Dynamic_bitset
   */
  Dynamic_bitset& reset() noexcept {
    idx bitset_storage_size = this->storage_size();
    for (idx i = 0; i < bitset_storage_size; ++i) {
      v_[i] = 0;
    }
    return *this;
  }

  /**
   * @brief flip the value of the bit at the specified position
   * @param pos the position of the bit to flip
   * @return reference to the modified Dynamic_bitset
   */
  Dynamic_bitset& flip(idx pos) {
    v_[index_(pos)] ^= 1 << (offset_(pos));
    return *this;
  }

  /**
   * @brief flip the value of all bits in the bitset
   * @return a reference to the modified Dynamic_bitset
   */
  Dynamic_bitset& flip() noexcept {
    idx bitset_storage_size = this->storage_size();
    for (idx i = 0; i < bitset_storage_size; ++i) {
      v_[i] = ~v_[i];
    }
    return *this;
  }

  /**
   * @brief compare two Dynamic_bitset object for equality
   * @param rhs the right-hand side Dynamic_bitset to compare with.
   * @return true if the bitset are equal, false otherwise
   */
  bool operator==(const Dynamic_bitset& rhs) const noexcept {
    assert(this->size() == rhs.size());
    bool result = true;
    idx n = std::min(this->storage_size(), rhs.storage_size());
    for (idx i = 0; i < n; ++i) {
      if (v_[i] != rhs.v_[i]) {
        return false;
      }
    }
    return result;
  }

  /**
   * @brief compare two Dynamic_bitset object for inequality
   * @param ths the right-hand side Dynamic_bitset to compare with
   * @return true if the bitset are not equal, false if the are equal.
   */
  bool operator!=(const Dynamic_bitset& rhs) const noexcept { return !(*this == rhs); }

  /**
   * @brief calculate the hamming distance between two Dynamic_bitset objects
   * @param rhs right-hand side Dynamic_bitset to calculate the distance to
   * @return the hamming distance between the bitset
   */
  idx operator-(const Dynamic_bitset& rhs) const noexcept {
    idx result = 0;
    idx bitset_size = this->size();
    for (idx i = 0; i < bitset_size; ++i) {
      if (this->get(i) != rhs.get(i))
        ++result;
    }
    return result;
  }

  /**
   * @brief convert the Dynamic_bitset to a string representation
   * @tparam CharT are character type for the string (default is 'char')
   * @tparam Traits the character traits type (default is 'std::char_traits<CharT>')
   * @tparam Allocator the allocator type for memory management (default is 'std::allocator<CharT>')
   * @param zero the character representing a bit value of 0
   * @param one the character representing a bit value of 1
   * @return A string representation of the dynamic bitset
   *
   * NOTE: this function converts the Dynamic_bitset to a string representation using the specified
   * characters for bit values 0 and 1. the resulting string is constructed by iterating through the
   * bits in the Dynamic_bitset and appending the correspoinding character to the string
   */
  template <class CharT = char, class Traits = std::char_traits<CharT>,
            class Allocator = std::allocator<CharT>>
  std::basic_string<CharT, Traits, Allocator> to_string(CharT zero = CharT('0'),
                                                        CharT one = CharT('1')) const {
    std::basic_string<CharT, Traits, Allocator> result;
    idx bitset_size = this->size();
    result.resize(bitset_size);

    // iterate through the bits in reverse order and append characters to the string
    for (idx i = bitset_size; i-- > 0;) {
      if (!this->get(i)) {
        result[bitset_size - i - 1] = zero;
      } else {
        result[bitset_size - i - 1] = one;
      }
    }
    return result;
  }

 private:
  /**
   * @brief display the binary representation of Dynamic_bitset to an output stream
   * @param os the output stream to write the binary repersentation to
   * @return the update output stream after writing the binary representation
   *
   * @override IDisplay::display
   *
   * NOTE: this private member function is used to display the binary repersentation of the
   * Dynamic_bitset to the specified output stream. it iterates through the bits in reverse
   * order and writes each bit's value (0 or 1) to the output stream. this function is intended
   * to be used as an override for the IDisplay::display function
   */
  std::ostream& display(std::ostream& os) const override {
    idx bitset_size = this->size();
    for (idx i = bitset_size; i-- > 0;) {
      os << this->get(i);
    }
    return os;
  }
};

/**
 * @brief A class repersenting a circuit of bit operations, derived from Dynamic_bitset
 *
 * NOTE: this class extend the functionally of Dynamic_bitset to repersent a circuit of bit.
 * it includes a struct Gate_count to keep track of different gate counts within the circuit
 */
class Bit_circuit : public Dynamic_bitset {
 public:
  /**
   * @brief a structure to keep track of gate counts within the circuit
   */
  struct Gate_count {
    // count of NOT Gates
    idx NOT = 0;
    // alias for the count of X gates (same as NOT gates)
    idx& X = NOT;

    // count CNOT gates
    idx CNOT = 0;
    // count SWAP gates
    idx SWAP = 0;
    // count FREDKIN gates
    idx FRED = 0;
    // count of toffoli gates
    idx TOF = 0;
  } gate_count{};

  /**
   * @brief constructor that forwards arguments to the base Dynamic_bitset constructor
   * @tparam Args parameter for the Dynamic_bitset constructor
   * @param args erguments to forward to the Dynamic_bitset constructor
   */
  using Dynamic_bitset::Dynamic_bitset;
  Bit_circuit(const Dynamic_bitset& dynamic_bitset) : Dynamic_bitset{dynamic_bitset} {}
  Bit_circuit& X(idx pos) {
    this->flip(pos);
    ++gate_count.X;
    return *this;
  }
  Bit_circuit& NOT(idx pos) {
    this->flip(pos);
    ++gate_count.NOT;
    return *this;
  }

  /**
   * @brief apply CNOT gate to the circuit
   * @param pos A vector two position: control bit and target bit
   * @return a reference to the modified Bit_circuit
   */
  Bit_circuit& CNOT(const std::vector<idx>& pos) {
    v_[index_(pos[1])] ^= (1 & (v_[index_(pos[0])] >> offset_(pos[0]))) << offset_(pos[1]);
    ++gate_count.CNOT;
    return *this;
  }

  /**
   * @brief apply toffoli gate to the circuit
   * @param pos A vector of three positions: control bit 1, control bit 2, and target bit
   * @return reference to the modified Bit_circuit
   */
  Bit_circuit& TOF(const std::vector<idx>& pos) {
    v_[index_(pos[2])] ^= ((1 & (v_[index_(pos[1])] >> offset_(pos[1]))) &
                           (1 & (v_[index_(pos[0])] >> offset_(pos[0]))))
                          << offset_(pos[2]);
    ++gate_count.TOF;
    return *this;
  }

  /**
  * @brief apply a SWAP gate to the circuit
  * @param pos A vector of two position to SWAP
  * @return A reference to the modified Bit_circuit
  *
  * NOTE: function swaps the value of the two specified position in the circuit
  */
  Bit_circuit& SWAP(const std::vector<idx>& pos) {
    if (this->get(pos[0]) != this->get(pos[1])) {
      this->X(pos[0]);
      this->X(pos[0]);
      this->X(pos[1]);
    }
    ++gate_count.SWAP;
    return *this;
  }

  /**
  * @brief apply FREDKIN gate to the circuit
  * @param pos A vector of three position: control bit, target bit 1, target bit 2
  * @return reference to the modified Bit_circuit
  *
  * NOTE: if the control bit is set, this function swaps the values of target bits 1 and 2 in the
  * circuit
  */
  Bit_circuit& FRED(const std::vector<idx>& pos) {
    if (this->get(pos[0])) {
      this->SWAP({pos[1], pos[2]});
    }
    ++gate_count.FRED;
    return *this;
  }

  /**
  * @brief reset the circuit and gat counts to their initial state
  * @return reference to the modified Bit_circuit
  *
  * NOTE: this function reset all bits in the circuit to zero and reset the gate count
  * to zero
  */
  Bit_circuit& reset() noexcept {
    gate_count.NOT = gate_count.X = 0;
    gate_count.CNOT = gate_count.SWAP = 0;
    gate_count.FRED = gate_count.TOF = 0;
    Dynamic_bitset::reset();

    return *this;
  }
};
}  // namespace clara

#endif  // !CLASSFUNCTION_REVESIBLE_H_
