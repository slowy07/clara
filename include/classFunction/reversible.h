// Copyright (c) 2023 arfy slowy
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef CLASSFUNCTION_REVESIBLE_H_
#define CLASSFUNCTION_REVESIBLE_H_

#include <strings.h>

#include <bitset>
#include <cassert>
#include <climits>
#include <locale>
#include <optional>
#include <ostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "../types.h"
#include "idisplay.h"
#include "interface_json.h"

namespace clara {

/**
 * @brief a dynamic bitset calss that implements a flexible-size arrays of bits
 * @details this class represent a dynamic bitset, which is an array of bits
 * that can dynamically grow as needed. it is used for efficient storage and manipulation
 * of binary data
 */
class Dynamic_bitset : public InterfaceDisplay {
 public:
  // type used to store individual bits
  using value_type = unsigned int;
  // type for the storage container
  using storage_type = std::vector<value_type>;

 protected:
  // number of storage elements needed
  idx storage_size_;
  // total number of bits
  idx n_;
  // storage vector for the bitset
  // std::vector<value_type>
  std::vector<value_type> v_;

  /**
   * @brief calculate the index of the storage element containing a bit at the given position
   * @param pos the position of the bit
   * @return the index of the storage element
   */
  inline static idx index_(idx pos) { return pos / (sizeof(value_type) * CHAR_BIT); }
  /**
   * @brief calculate the offset within the storage element fora bit at the given position
   * @param pos the position of the bit
   * @return the offset within the storage element
   */
  inline static idx offset_(idx pos) { return pos % (sizeof(value_type) * CHAR_BIT); }

 public:
  /**
   * @brief constructor to create a dynamic bitset of a given size
   * @param N the total number of bits in the bitset
   */
  explicit Dynamic_bitset(idx n)
      : storage_size_{static_cast<idx>(n / (sizeof(value_type) * CHAR_BIT) + 1)},
        n_{n},
        v_(storage_size_) {}

  /**
   * @brief construct dynamic bitset from string of binary digits
   *
   * string is interpreted as a squance of binary characters,
   * where each character representing single bit, by default '0' is interpreted
   * as false (bit 0), and '1' as true (bit 1)
   *
   * @param str string contains binary digit character
   * @param zero cahracter used to representing 0-bit
   * @param one character used to representing a 1-bit
   */
  explicit Dynamic_bitset(std::string str, char zero = '0', [[maybe_unused]] char one = '1')
      : Dynamic_bitset{static_cast<idx>(str.size())} {
    idx n = str.size();
    // iterate over each character in the string
    for (idx i = 0; i < n; ++i) {
      // ensure each character is either zero or one
      assert(str[i] == zero || str[i] == one);
      // set the bit at position (n - 1 - i) to match string value
      // this makse sure the first character in the string treated
      // as the MSB
      this->set(n - 1 - i, str[i] != zero);
    }
  }

  ~Dynamic_bitset() override = default;

  /**
   * @brief get a reference to the internal stroage internal container
   * @return a const reference to the internal storage vector
   */
  const storage_type& data() const { return v_; }
  /**
   * @brief get a reference to the internal storage container
   * @return const reference to the internal storage vector
   */
  idx size() const noexcept { return n_; }
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
  bool get(idx pos) const noexcept {
    assert(pos < size());
    return 1 & (v_[index_(pos)] >> offset_(pos));
  }

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
  bool any() const noexcept { return !(none()); }

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
  virtual Dynamic_bitset& reset() noexcept {
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
  virtual std::string to_string(char zero = '0', char one = '1') const {
    std::string result;
    idx bitset_size = size();
    result.resize(bitset_size);

    for (idx i = bitset_size; i-- > 0;) {
      result[bitset_size - i - 1] = get(i) ? one : zero;
    }
    return result;
  }

 protected:
  /**
   * @brief display the binary representation of Dynamic_bitset to an output stream
   * @param os the output stream to write the binary repersentation to
   * @return the update output stream after writing the binary representation
   *
   * @override InterfaceDisplay::display
   *
   * NOTE: this private member function is used to display the binary repersentation of the
   * Dynamic_bitset to the specified output stream. it iterates through the bits in reverse
   * order and writes each bit's value (0 or 1) to the output stream. this function is intended
   * to be used as an override for the InterfaceDisplay::display function
   */
  std::ostream& display(std::ostream& os) const override { return os << this->to_string(); }
};

/**
 * @brief A class repersenting a circuit of bit operations, derived from Dynamic_bitset
 *
 * NOTE: this class extend the functionally of Dynamic_bitset to repersent a circuit of bit.
 * it includes a struct Gate_count to keep track of different gate counts within the circuit
 */
class Bit_circuit : public Dynamic_bitset, public InterfaceJson {
  std::unordered_map<std::string, idx> count_{};
  std::unordered_map<std::string, idx> depth_{};
  Dynamic_bitset bNOT_, bCNOT_, bSWAP_, bTOF_, bFRED_, btotal_;

 public:
  /**
   * @brief construct a bit with specified number of zero initialized
   *
   * initialized all operation counters and tracker to default value
   *
   * @param n number of bit in the circuit
   */
  explicit Bit_circuit(idx n)
      : Dynamic_bitset{n}, bNOT_{n}, bCNOT_{n}, bSWAP_{n}, bTOF_{n}, bFRED_{n}, btotal_{n} {}

  /**
   * @brief construct a Bit_circuit from a binary string
   *
   * parse a string where each character representing a bit value
   * and initialize the bitset accordingly
   */
  explicit Bit_circuit(std::string str, char zero = '0', [[maybe_unused]] char one = '1')
      : Dynamic_bitset{static_cast<idx>(str.size())},
        bNOT_{size()},
        bCNOT_{size()},
        bSWAP_{size()},
        bTOF_{size()},
        bFRED_{size()},
        btotal_{size()} {
    idx n = str.size();
    for (idx i = 0; i < n; ++i) {
      assert(str[i] == zero || str[i] == one);
      this->set(i, str[i] != zero);
    }
  }

  /**
   * @brief construct a Bit_circuit from binary string
   *
   * parsing string where each character representing bit value
   * the leftomst character correponding to the first bit
   *
   * @param str string containign binary digit char
   * @param zero character used to representing 0bit
   * @param one character used to representing a 1bit
   */
  explicit Bit_circuit(const Dynamic_bitset& dynamic_bitset)
      : Dynamic_bitset{dynamic_bitset},
        bNOT_{size()},
        bCNOT_{size()},
        bSWAP_(size()),
        bTOF_{size()},
        bFRED_{size()},
        btotal_{size()} {}

  /**
   * @brief convenience wrapper for apply not gate
   *
   * @param i index of the bit to flip
   * @return reference to this object for method chaining
   */
  Bit_circuit& X(idx i) {
    assert(i < size());
    NOT(i);
    return *this;
  }

  /**
   * @brief destructor
   *
   * virtual to allow safe inheritance and polymorphism
   */
  ~Bit_circuit() override = default;

  /**
   * @brief applies NOT gate to specified bit
   *
   * flips the bit at index `i`, updateing the NOT tracker and icrement counter
   *
   * @param i index of the bit to flip
   * @return reference to this object for method chaining
   */
  Bit_circuit& NOT(idx i) {
    assert(i < size());
    flip(i);  // flip bit at position i

    ++count_["NOT"];
    ++count_[__FILE__ "__total__"];

    if (count_["NOT"] == 1) {
      depth_["NOT"] = 1;
    }

    bNOT_.flip(i);

    // if the bit has been reset after being flip, increment global depth
    if (!bNOT_.get(i)) {
      bNOT_ = Dynamic_bitset{n_};      // reset tracker
      btotal_.set(i);                  // mark as part of new depth layer
      ++depth_[__FILE__ "__total__"];  // increment total depth counter
    }
    return *this;
  }

  Bit_circuit& CNOT(idx ctrl, idx target) {
    [[maybe_unused]] auto n = size();
    assert(ctrl < n && target < n);
    assert(ctrl != target);

    v_[index_(target)] ^= (static_cast<value_type>(1) & (v_[index_(ctrl)] >> offset_(ctrl)))
                          << offset_(target);
    ++count_["CNOT"];
    ++count_[__FILE__ "__total__"];

    if (count_["CNOT"] == 1) {
      depth_["CNOT"] = 1;
    }

    bCNOT_.flip(ctrl).flip(target);
    if (!bCNOT_.get(ctrl) || !bCNOT_.get(target)) {
      bCNOT_ = Dynamic_bitset{n_};

      bCNOT_.set(ctrl).set(target);
      ++depth_["CNOT"];
    }

    if (count_[__FILE__ "__total__"] == 1) {
      depth_[__FILE__ "__total__"] = 1;
    }

    btotal_.flip(ctrl).flip(target);

    if (!btotal_.get(ctrl) || !btotal_.get(target)) {
      btotal_ = Dynamic_bitset{n_};
      btotal_.set(ctrl).set(target);
      ++depth_[__FILE__ "__total__"];
    }
    return *this;
  }

  Bit_circuit& TOF(idx i, idx j, idx k) {
    [[maybe_unused]] auto n = size();
    assert(i < n && j < n && k < n);

    v_[index_(k)] ^= ((static_cast<value_type>(1) && (v_[index_(j)] >> offset_(j))) &
                      (static_cast<value_type>(1) & (v_[index_(i)] >> offset_(i))))
                     << offset_(k);
    ++count_["TOF"];
    ++count_[__FILE__ "__total__"];

    if (count_["TOF"] == 1) {
      depth_["TOF"] = 1;
    }

    bTOF_.flip(i).flip(j).flip(k);

    if (!bTOF_.get(i) || !bTOF_.get(j) || !bTOF_.get(k)) {
      bTOF_ = Dynamic_bitset{n_};

      bTOF_.set(i).set(j).set(k);
      ++depth_["TOF"];
    }

    if (count_[__FILE__ "__total__"] == 1) {
      depth_[__FILE__ "__total__"] = 1;
    }

    btotal_.flip(i).flip(j).flip(k);

    if (!btotal_.get(i) || !btotal_.get(j) || !btotal_.get(k)) {
      btotal_ = Dynamic_bitset{n_};
      btotal_.set(i).set(j).set(k);
      ++depth_[__FILE__ "__total__"];
    }
    return *this;
  }

  Bit_circuit& SWAP(idx i, idx j) {
    [[maybe_unused]] auto n = size();
    assert(i < n && j < n);
    assert(i != j);

    if (get(i) != get(j)) {
      X(i);
      X(j);
    }
    ++count_["SWAP"];
    ++count_[__FILE__ "__total__"];

    if (count_["SWAP"] == 1) {
      depth_["SWAP"] = 1;
    }

    btotal_.flip(i).flip(j);
    if (!btotal_.get(i) || !btotal_.get(j)) {
      btotal_ = Dynamic_bitset{n_};
      btotal_.set(i).set(j);
      ++depth_[__FILE__ "__total__"];
    }
    return *this;
  }

  Bit_circuit& FRED(idx i, idx j, idx k) {
    [[maybe_unused]] auto n = size();
    assert(i < n && j < n && k < n);
    assert(i != j && i != k && j != k);

    if (get(i)) {
      SWAP(j, k);
    }
    ++count_["FRED"];
    ++count_[__FILE__ "__total__"];

    if (count_["FRED"] == 1) {
      depth_["FRED"] = 1;
    }

    bFRED_.flip(i).flip(j).flip(k);
    if (!bFRED_.get(i) || !bFRED_.get(j) || !bFRED_.get(k)) {
      bFRED_ = Dynamic_bitset{n_};
      bFRED_.set(i).set(j).set(k);
      ++depth_["FRED"];
    }

    if (count_[__FILE__ "__total__"] == 1) {
      depth_[__FILE__ "__total__"] = 1;
    }

    btotal_.flip(i).flip(j).flip(k);
    if (!btotal_.get(i) || !btotal_.get(j) || !btotal_.get(k)) {
      btotal_ = Dynamic_bitset{n_};
      btotal_.set(i).set(j).set(k);
      ++depth_[__FILE__ "__total__"];
    }
    return *this;
  }

  Bit_circuit& reset() noexcept override { return *this; }

  idx get_gate_count(std::optional<std::string> name = std::nullopt) const {
    idx result;

    if (name.has_value()) {
      try {
        result = (name == "X") ? count_.at("NOT") : count_.at(name.value());
      } catch (...) {
        return 0;
      }
    } else {
      try {
        result = count_.at(__FILE__ "__total__");
      } catch (...) {
        return 0;
      }
    }
    return result;
  }

  idx get_gate_depth(std::optional<std::string> name = std::nullopt) const {
    idx result;

    if (name.has_value()) {
      try {
        result = (name == "X") ? depth_.at("NOT") : depth_.at(name.value());
      } catch (...) {
        return 0;
      }
    } else {
      try {
        result = depth_.at(__FILE__ "__total__");
      } catch (...) {
        return 0;
      }
    }
    return result;
  }

  std::string to_JSON(bool enclosed_in_curly_brackets = true) const override {
    std::string result;

    if (enclosed_in_curly_brackets) {
      result += "{";
    }
    result += "\"n\" : " + std::to_string(n_);
    result += ", \"total gate count\" : " + std::to_string(get_gate_count());
    result += ", \"total gate depth\" : " + std::to_string(get_gate_depth());
    result += R"(, "bit state" : ")" + this->to_string() + '\"';
    result += ", \"Hamming weight\" : " + std::to_string(count());
    if (enclosed_in_curly_brackets) {
      result += "}";
    }
    return result;
  }

  std::string to_string(char zero = '0', char one = '1') const override {
    std::string result;
    idx bitset_size = size();
    result.resize(bitset_size);

    for (idx i = 0; i < bitset_size; ++i) {
      result[i] = get(i) ? one : zero;
    }
    return result;
  }

 private:
  std::ostream& display(std::ostream os) const {
    os << "n = " << n_ << '\n';
    os << "total gate count: " << get_gate_count() << '\n';
    os << "total gate depth: " << get_gate_depth() << '\n';
    os << "bit state: " << this->to_string() << '\n';
    os << "Hamming weight: " << count();
    return os;
  }
};
}  // namespace clara

#endif  // !CLASSFUNCTION_REVESIBLE_H_
