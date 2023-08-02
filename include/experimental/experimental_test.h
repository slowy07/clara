#ifndef EXPERIMENTAL_EXPERIMENTAL_TEST_H_
#define EXPERIMENTAL_EXPERIMENTAL_TEST_H_

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstddef>
#include <iterator>
#include <ostream>
#include <random>
#include <string>
#include <vector>

using idx = std::size_t;

namespace clara {
namespace experimental {

/**
 * @class Dynamic_bitset
 * @brief Dynamic bitest class, allows the specification of the number of number bits at runtime
 *
 * the Dynamic_bitset class is a dynamic bitset that allows the specification of the number of bits
 * at runtime. it provides functionalities to manipulate individual bits, set random bits, and
 * perform bit operations.
 */
class Dynamic_bitset {
 public:
  using value_type = unsigned int;
  using storage_type = std::vector<value_type>;

 protected:
  // storage size
  idx storage_size_;
  // number of bits
  idx N_;
  // storage space
  std::vector<value_type> v_;

  /**
   * @brief index of the pos bit in the storage space
   * @return index of the pos bit in thestorage space
   */
  idx index_(idx pos) const { return pos / (sizeof(value_type) * CHAR_BIT); }

  /**
   * @brief constructor, intialize all bits to false
   * @return offset of the pos bit in the storage space relative to its index
   */
  idx offset_(idx pos) const { return pos % (sizeof(value_type) * CHAR_BIT); }

  /**
   * @brief offset of the pos bit in the storage space relativea
   *       to its index
   * @return offset of the pos bit in the storage space relative its
   *         index
   */
 public:
  Dynamic_bitset(idx N)
      : storage_size_{N / (sizeof(value_type) * CHAR_BIT) + 1}, N_{N}, v_(storage_size_) {}

  /**
   * @brief raw storage space of the bitset
   * @return const reference to the underlying
   */
  const storage_type& data() const { return v_; }

  /**
   * @brief number of bits stored in the bitset
   * @return number of bits
   */
  idx size() const { return N_; }

  /**
   * @brief size of the underlying stroage space
   * @return size of the underlying storage space
   */
  idx storage_size() const { return storage_size_; }

  idx count() const noexcept {
    std::size_t result = 0;
    for (idx i = 0; i < size(); ++i) {
      if (this->get(i))
        ++result;
    }
    return result;
  }

  /**
   * @brief get the value of the bit at the given position
   * @param pos position of the bit to check
   * @return true if the bit set (1), false otherwise (0)
   */
  bool get(idx pos) const { return 1 & (v_[index_(pos)] >> offset_(pos)); }

  /**
   * @dbrief check if none of the bits are set in the bitset
   * @return true if none of the bits are set, false otherwise
   */
  bool none() const noexcept {
    bool result = true;
    for (idx i = 0; i < storage_size(); ++i) {
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
    for (idx i = 0; i < storage_size(); ++i) {
      if (~v_[i]) {
        return false;
      }
    }
    return false;
  }

  /**
   * @brief check if any of the bits are set in the bitset
   * @return true if any of the bits are set, false otherwise
   */
  bool any() const noexcept { return !(this->none()); }

  /**
   * @brief set the bit at the given position to the specified value
   * @param pos position of the bitset
   * @param value value to set the bit
   * @return reference to the modified bitset
   */
  Dynamic_bitset& set(idx pos, bool value = true) {
    value ? v_[index_(pos)] |= (1 << offset_(pos)) : v_[index_(pos)] &= ~(1 << offset_(pos));
    return *this;
  }

  /*
   * @brief all bits in the bitset to the specified value (true = 1, false = 0)
   * @return reference to the modified bitset
   */
  Dynamic_bitset& set() noexcept {
    for (idx i = 0; i < storage_size(); ++i) {
      v_[i] = ~0;
    }
    return *this;
  }

  /**
   * @brief set the bit at the given position randomly according to a Bernoulli distribution
   * @param pos position of the bit to set
   * @param p probability of setting the bit to 1
   * @return reference to the modified bitset
   */
  Dynamic_bitset& rand(idx pos, double p = 0.5) {
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::bernoulli_distribution d{p};

    this->set(pos, d(gen));
    return *this;
  }

  /**
   * @brief set all bits in the bitset randomly according to a Bernoulli distribution
   * @param p probability of setting each bit to 1
   * @return reference to the modified bitset
   */
  Dynamic_bitset& rand(double p = 0.5) {
    for (idx i = 0; i < size(); ++i) {
      this->rand(i, p);
    }
    return *this;
  }

  /**
   * @brief reset the bit at the given position (set 0)
   * @param pos position of the bit to reset
   * @return reference to the modified bitset
   */
  Dynamic_bitset& reset(idx pos) {
    v_[index_(pos)] &= ~(1 << offset_(pos));
    return *this;
  }

  /**
   * @brief reset all bits in the bitset
   * @return reference to the modified bitset
   */
  Dynamic_bitset& reset() noexcept {
    for (idx i = 0; i < storage_size(); ++i) {
      v_[i] = 0;
    }
    return *this;
  }

  /**
   * @brief flip the bit at the given position
   * @param position orf the bit to flip
   * @return reference to the modified bitset
   */
  Dynamic_bitset& flip(idx pos) {
    v_[index_(pos)] ^= 1 << (offset_(pos));
    return *this;
  }

  /**
   * @brief equality comparsion operator for bitsets
   * @param rhs bitset to compare with rhs
   * @return true if the bitsets are equal, false otherwise
   */
  Dynamic_bitset& flip() noexcept {
    for (idx i = 0; i < storage_size(); ++i) {
      v_[i] = ~v_[i];
    }
    return *this;
  }

  /**
   * @brief equality comprasion operator for bitsets
   * @param rhs bitset to compare with
   * @return true if the bitsets are equal, false otherwise
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
   * @brief inequality comprasion oeprator for bitset
   * @param rhs bitset to compare with
   * @return true if the bitsets are not equal, false otherwise
   */
  bool operator!=(const Dynamic_bitset& rhs) const noexcept { return !(*this == rhs); }

  /**
   * @brief output stream operator to print the bitset in binary representation
   * @param os output stream
   * @param rhs bitset to print
   * @return reference to the output stream
   */
  friend std::ostream& operator<<(std::ostream& os, const Dynamic_bitset& rhs) {
    for (idx i = rhs.size(); i-- > 0;) {
      os << rhs.get(i);
    }
    return os;
  }

  /**
   * @brief convert the bitset to a string representation
   * @param zero character to represent 0 bit
   * @param one character to represent 1 bit
   * @return string representation of the bitset
   */
  template <class CharT = char, class Traits = std::char_traits<CharT>,
            class Allocator = std::allocator<CharT>>
  std::basic_string<CharT, Traits, Allocator> to_string(CharT zero = CharT('0'),
                                                        CharT one = CharT('1')) const {
    std::basic_string<CharT, Traits, Allocator> result;
    idx bitset_size = this->size();
    result.resize(bitset_size);

    for (idx i = bitset_size; i-- > 0;) {
      if (!this->get(i)) {
        result[bitset_size - i - 1] = zero;
      } else {
        result[bitset_size - i - 1] = one;
      }
    }
    return result;
  }
};

/**
 * @class Bit_circuit
 * @brief classical revesible circuit simulator
 *
 * the Bit_circuit class is a classical
 */
class Bit_circuit : public Dynamic_bitset {
 public:
  struct Gate_count {
    // 1 bit gates
    idx NOT = 0;
    idx& X = NOT;

    // 2 bit gates
    idx CNOT = 0;
    idx SWAP = 0;

    // 3 bit gates
    idx FRED = 0;
    idx TOF = 0;
  } gate_count{};

  // inherit the constructor
  using Dynamic_bitset::Dynamic_bitset;

  // flip the bit
  Bit_circuit& X(idx pos) {
    this->flip(pos);
    ++gate_count.X;
    return *this;
  }

  // flips the bit
  Bit_circuit& NOT(idx pos) {
    this->flip(pos);
    ++gate_count.NOT;
    return *this;
  }

  // toffoli control-control-target
  Bit_circuit& TOF(const std::vector<idx>& pos) {
    v_[index_(pos[2])] ^= ((1 & (v_[index_(pos[1])] >> offset_(pos[1]))) &
                           (1 & (v_[index_(pos[0])] >> offset_(pos[0]))))
                          << offset_(pos[2]);
    ++gate_count.TOF;
    return *this;
  }

  // SWAP 2 bits
  Bit_circuit& SWAP(const std::vector<idx>& pos) {
    if (this->get(pos[0]) != this->get(pos[1])) {
      this->X(pos[0]);
      this->X(pos[1]);
    }
    ++gate_count.SWAP;
    return *this;
  }

  Bit_circuit& reset() noexcept {
    gate_count.NOT = gate_count.X = 0;
    gate_count.CNOT = gate_count.SWAP = 0;
    gate_count.FRED = gate_count.TOF = 0;

    return *this;
  }
};

}  // namespace experimental
}  // namespace clara

#endif  // !EXPERIMENTAL_EXPERIMENTAL_TEST_H_
