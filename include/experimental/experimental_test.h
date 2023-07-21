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
 * @class clara::Dynamic_bitset
 * @brief dynamic bitset class, allow the specification of
 *         the number of bits at runtime
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

  bool get(idx pos) const { return 1 & (v_[index_(pos)] >> offset_(pos)); }
  bool none() const noexcept {
    bool result = true;
    for (idx i = 0; i < storage_size(); ++i) {
      if (v_[i]) {
        return false;
      }
    }
    return result;
  }

  bool all() const noexcept {
    bool result = true;
    for (idx i = 0; i < storage_size(); ++i) {
      if (~v_[i]) {
        return false;
      }
    }
    return false;
  }

  bool any() const noexcept { return !(this->none()); }

  Dynamic_bitset& set(idx pos, bool value = true) {
    value ? v_[index_(pos)] |= (1 << offset_(pos)) : v_[index_(pos)] &= ~(1 << offset_(pos));
    return *this;
  }

  Dynamic_bitset& set() noexcept {
    for (idx i = 0; i < storage_size(); ++i) {
      v_[i] = ~0;
    }
    return *this;
  }

  // set the bit according to a random bernouli distribution
  Dynamic_bitset& rand(idx pos, double p = 0.5) {
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::bernoulli_distribution d{p};

    this->set(pos, d(gen));
    return *this;
  }

  // set all bits according to a random bernouli distribution
  Dynamic_bitset& rand(double p = 0.5) {
    for (idx i = 0; i < size(); ++i) {
      this->rand(i, p);
    }
    return *this;
  }

  // set bit false
  Dynamic_bitset& reset(idx pos) {
    v_[index_(pos)] &= ~(1 << offset_(pos));
    return *this;
  }

  // set all bits 0
  Dynamic_bitset& reset() noexcept {
    for (idx i = 0; i < storage_size(); ++i) {
      v_[i] = 0;
    }
    return *this;
  }

  // flips the bit
  Dynamic_bitset& flip(idx pos) {
    v_[index_(pos)] ^= 1 << (offset_(pos));
    return *this;
  }

  Dynamic_bitset& flip() noexcept {
    for (idx i = 0; i < storage_size(); ++i) {
      v_[i] = ~v_[i];
    }
    return *this;
  }

  // operator
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

  bool operator!=(const Dynamic_bitset& rhs) const noexcept { return !(*this == rhs); }

  friend std::ostream& operator<<(std::ostream& os, const Dynamic_bitset& rhs) {
    for (idx i = rhs.size(); i-- > 0;) {
      os << rhs.get(i);
    }
    return os;
  }

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
 * @class clara::Bit_circuit
 * @brief classical reversible circuit simulator
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
