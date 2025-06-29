#ifndef CLASSFUNCTION_LAYOUT_H_
#define CLASSFUNCTION_LAYOUT_H_

#include <vector>

#include "../functions.h"
#include "../internal/util.h"
#include "../types.h"
#include "exception"
#include "exception.h"

namespace clara {

/**
 * @class InterfaceLayout
 * @brief abstract interface for mapping between flat indices and muti-dimensional
 *        coordinate
 *
 * this class serve base class for implement different layout of multi dimension
 * such as row-major, column-major, or custom ordering
 */
class InterfaceLayout {
 public:
  /**
   * @brief convert multi-dimensional coordinate (a vector indices) to flat linear
   *
   * @param xs vector representing position in each dimension
   * @return corresponding flat index
   */
  virtual idx operator()(const std::vector<idx>& xs) const = 0;

  /**
   * @brief converrt flat linear index into its corresponding multi-dimensional coordinate
   *
   * @param i flat index to convert
   * @return vector representing the coordinate in multi-dimensional space
   */
  virtual std::vector<idx> to_coordinate(idx i) const = 0;

  /**
   * @brief retunr the dimension of the layout
   *
   * @return vector of dimension size
   */
  virtual std::vector<idx> get_dims() const = 0;

  /**
   * @brief virtual destructor
   *
   * make sure proper cleanup of derived object class
   */
  virtual ~InterfaceLayout() = default;
};

/**
 * @class Lattice
 * @brief concrete implementation of InterfaceLayout for representing multi-dimensional
 *
 * example
 *
 * given dimension = {2, 3}, the mapping would be
 * (0, 0) -> 0  (0, 1) -> 1  (0, 2) -> 2
 * (1, 0) -> 3  (1, 1) -> 4  (1, 2) -> 5
 */
class Lattice : public InterfaceLayout {
 public:
  /**
   * @brief construct lattice from a vector of dimension
   *
   * @param dims vector specifying the size of each dimension
   * @throws exception::ZeroSize if dims is empty
   */
  explicit Lattice(const std::vector<idx>& dims) : dims_{dims} {
    if (dims.empty()) {
      throw exception::ZeroSize("clara::Lattice::Lattice()", "dims");
    }
  }

  /**
   * @brief construct a lattice from a variadic list of dimension
   *
   * @tparam Ts variadic tempalte types convertible to idx
   * @param ds parameter pack of dimension sizes
   *
   * example
   * ```
   * Lattice lattice(2, 3, 4); // this will create 3d lattice size 2x3x4
   * ```
   */
  template <typename... Ts>
  explicit Lattice(Ts... ds) : clara::Lattice(std::vector<idx>{static_cast<idx>(ds)...}) {}

  /**
   * @brief convert a multi-dimensional coordinate to flat linear index
   *
   * @param xs vector representing the coordinate in each dimension.
   * @return linear index corresponding to the coordinate.
   */
  idx operator()(const std::vector<idx>& xs) const override {
    if (xs.size() != dims_.size()) {
      throw exception::SizeMismatch("clara::Lattice::operator", "xs");
    }

    for (idx i = 0; i < static_cast<idx>(dims_.size()); ++i) {
      if (xs[i] >= dims_[i]) {
        throw exception::OutOfRange("clara::Lattice::operator()", "xs");
      }
    }
    return internal::multiidx2n(xs.data(), dims_.size(), dims_.data());
  }

  /**
   * @brief variadic overload of operator() for convertion coordinate flat index
   *
   * @tparam Ts types convertible to idx
   * @param xs variadic list of indices
   * @return flat index corresponding to the coordinate
   */
  template <typename... Ts>
  idx operator()(Ts... xs) const {
    return operator()(std::vector<idx>{static_cast<idx>(xs)...});
  }

  /**
   * @brief convert flat index to its corresponding multi-dimensional coordinate
   *
   * @param i linear index to convert
   * @return vector representing the coordinate in each dimension
   */
  std::vector<idx> to_coordinate(idx i) const override {
    if (i >= prod(dims_)) {
      throw exception::OutOfRange("clara::Lattice::to_coordinate()", "i");
    }

    std::vector<idx> result(dims_.size());
    internal::n2multiidx(i, dims_.size(), dims_.data(), result.data());

    return result;
  }

  /**
   * @brief return the dimension of the lattice
   *
   * @return vector containing the size of each dimension
   */
  std::vector<idx> get_dims() const override { return dims_; }

 protected:
  std::vector<idx> dims_;
};

/**
 * @class PeriodicBoundaryLattice
 * @brief a lattice with periodic boundary conditions
 *
 * inherit from clara::Lattice and overriding its indexing behaviour
 * so that out-of-bounds coordinate are wrapped around using modulo
 * arithmetic.
 */
class PeriodicBoundaryLattice : public Lattice {
 public:
  // allow PeriodicBoundaryLattice to be construct using any valid constructor
  // from its parent class Lattice such as
  // - from vector dimenson
  // - from parameter pack of dimension size
  using clara::Lattice::Lattice;
  // inherit operator() from lattice from coordinate -> index mapping
  using clara::Lattice::operator();

  /**
   * @brief convert possibly out-of-bounds coordinate to flat linear index with wrapping
   *
   * @param xs vector representing the coordinate in each dimension
   * @return linear index corresponding to the wrapped coordinate
   *
   */
  idx operator()(const std::vector<idx>& xs) const override {
    // validating input size matches lattice rank
    if (xs.size() != dims_.size()) {
      throw exception::SizeMismatch("clara::Lattice::operator()", "xs");
    }

    // make copy of input coordinate to modify them
    std::vector<idx> xs_copy = xs;
    // apply modulo operation to each coordinate to enforce periodic boundaries
    for (idx i = 0; i < static_cast<idx>(dims_.size()); ++i) {
      xs_copy[i] = xs[i] % dims_[i];
    }
    // delegate to base class implementation for final index calculation
    return Lattice::operator()(xs_copy);
  }
};

}  // namespace clara

#endif  // !CLASSFUNCTION_LAYOUT_H_
