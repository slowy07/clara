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

#ifndef INTERNAL_CLASSFUNCTION_IOMANIP_H_
#define INTERNAL_CLASSFUNCTION_IOMANIP_H_

#include <complex>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../../classFunction/interface_display.h"
#include "../util.h"

namespace clara {
namespace internal {

/**
 * @class IOManipScalar
 * @brief display manipulator for scalar value with optional zero-chopping
 *
 * this class providing customizable formatting for scalar value when printing to
 * output stream.
 *
 * design for use in numerical lib where clean and precise output is important,
 * specially when dealing with quantum state, matrices, or tensor component
 */
template <typename Scalar>
class IOManipScalar : public InterfaceDisplay {
  Scalar scalar_;
  IOManipScalarOpts opts_;

 public:
  /**
   * @brief construct and IOManipScalar object
   *
   * @param scalar scalar value to format and display
   * @param opts display options controlling formatting behavior
   */
  IOManipScalar(Scalar scalar, IOManipScalarOpts opts) : scalar_{scalar_}, opts_{opts} {}

 private:
  /**
   * @brief display the scalar value to the given output stream
   *
   * if abs value of the scalar if below the chopping threshold,
   * it will display as `0` to suppress floating-point noise
   *
   * @param os output stream to write to
   * @return reference to the modified output stream
   */
  std::ostream& display(std::ostream& os) const override {
    if (!(std::abs(scalar_) < opts_.chop)) {
      os << scalar_;
    } else {
      os << 0;
    }
    return os;
  }
};

/**
 * @class IOManipScalar<std::complex<T>>
 * @brief Specialization of IOManipScalar for complex number
 *
 * - suppressing near-zero real or imag parts based on threshold
 * - custom delimiter around the complex number
 * - custom symbol for operator and imaginary unit suffix
 */
template <typename T>
class IOManipScalar<std::complex<T>> : public InterfaceDisplay {
  std::complex<T> z_;
  IOManipComplexOpts opts_;

 public:
  /**
   * @brief construct an IOManipScalar object for a complex number
   *
   * @param z complex number to format and display
   * @param opts display options controlling formatting behavior
   */
  IOManipScalar(std::complex<T> z, IOManipComplexOpts opts) : z_{z}, opts_{std::move(opts)} {}

 private:
  /**
   * @brief display complex number to the given output stream with custom formatting
   *
   * handling
   * - both access real and imaginary part are near zero -> print "0"
   * - only imaginary part is significant -> print "bi"
   * - only real part is significant -> print "a"
   * - both part are significant -> print "a += bi"
   *
   * @param os output stream to write to
   * @return reference to the modified output stream
   */
  std::ostream& display(std::ostream& os) const override {
    realT re = std::real(z_);
    realT im = std::imag(z_);
    realT abs_re = std::abs(re);
    realT abs_im = std::abs(im);

    os << opts_.left;

    // both real and imaginary parts are negligible
    if (abs_re < opts_.chop && abs_im < opts_.chop) {
      os << '0';
    }

    // only imaginary part is significant
    if (abs_re < opts_.chop && !(abs_im < opts_.chop)) {
      os << im << opts_.im_suffix;
    }

    // only real part significant
    if (!(abs_re < opts_.chop) && abs_im < opts_.chop) {
      os << re;
    }

    // both part are significant
    if (!(abs_re < opts_.chop) && !(abs_im < opts_.chop)) {
      os << re;
      // print operator '+' or '-'
      os << (im < 0 ? opts_.minus_op : opts_.plus_op);
      os << abs_im << opts_.im_suffix;
    }
    // end with right delimiter
    os << opts_.right;
    return os;
  }
};

/**
 * @class IOManipRange
 * @brief display manipulator for printing range of values with customizable format
 *
 * wraps pair of input iterator and provides custom formatting when printinga
 * the element in that range to an output stream
 */
template <typename InputIterator>
class IOManipRange : public InterfaceDisplay {
  InputIterator first_, last_;
  IOManipRangeOpts opts_;

 public:
  /**
   * @brief construct an IOManipRange object from range and formatting options
   *
   * @param first operator pointing the begining of the range
   * @param last iterator pointing to the end of the range
   * @param opts options controlling how the range should be formatted
   */
  explicit IOManipRange(InputIterator first, InputIterator last, IOManipRangeOpts opts)
      : first_{first}, last_{last}, opts_{std::move(opts)} {}

  IOManipRange(const IOManipRange&) = default;
  IOManipRange& operator=(const IOManipRange&) = default;

 private:
  /**
   * @brief display range to the given output stream with custom formatting
   *
   * @param os output stream to write to
   * @return reference to the modified output stream
   */
  std::ostream& display(std::ostream& os) const override {
    os << opts_.left;

    // variabel to hold the separator string
    std::string sep;
    // iterate through the range and print each element
    for (InputIterator it = first_; it != last_; ++it) {
      os << sep << abs_float_or_cplx_chop(*it, opts_.chop);
      sep = opts_.sep;
    }
    // end with the right delimiter
    os << opts_.right;
    return os;
  }
};

/**
 * @class IOManipPointer
 * @brief display manipulator for printing array-like data via raw pointer
 *
 * wrap a raw pointer and length, and providing customizable formatting when
 * printing the pointed-to data to an output stream
 */
template <typename PointerType>
class IOManipPointer : public InterfaceDisplay {
  const PointerType* p_;
  idx N_;
  IOManipPointerOpts opts_;

 public:
  /**
   * @brief construct IOManipPointer object from pointer, length, and formatting options
   *
   * @param p pointer to the start of the data to display
   * @param N number of element to displaying
   * @param opts options controlling how the range should be formatted
   */
  explicit IOManipPointer(const PointerType* p, idx N, IOManipPointerOpts opts)
      : p_{p}, N_{N}, opts_{std::move(opts)} {}

  /**
   * @brief copy assign operator
   */
  IOManipPointerOpts& operator=(const IOManipPointer&) = default;

 private:
  std::ostream& display(std::ostream& os) const override {
    os << opts_.left;

    for (idx i = 0; i < N_ - 1; ++i) {
      os << abs_float_or_cplx_chop(p_[i], opts_.chop) << opts_.sep;
    }

    if (N_ > 0) {
      os << abs_float_or_cplx_chop(p_[N_ - 1], opts_.chop);
    }
    os << opts_.right;
    return os;
  }
};

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif  // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)

class IOManipEigen : public InterfaceDisplay, private Display_Impl_ {
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPLER) && (__GNUC__ == 4) && \
    (__GNUC_MINOR__ == 8)
#pragma diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPLER) && (__GNUC__ == 4)
        // && (__GNUC_MINOR__ == 8)
  cmat A_;
  IOManipEigenOpts opts_;

 public:
  template <typename Derived>
  explicit IOManipEigen(const Eigen::MatrixBase<Derived>& A, IOManipEigenOpts opts)
      : A_{A.template cast<cplx>()}, opts_{std::move(opts)} {}

 private:
  std::ostream& display(std::ostream& os) const override { return display_impl_(A_, os, opts_); }
};

template <typename Scalar>
class IOManipDirac : public InterfaceDisplay {
  clara::dirac_t<Scalar> A_;
  IOManipDiracOpts opts_{};

 public:
  explicit IOManipDirac(const dirac_t<Scalar>& A, IOManipDiracOpts opts)
      : A_{A}, opts_{std::move(opts)} {}

 private:
  void display_ket_dits_(std::ostream& os, const std::vector<idx>& dits) const {
    os << "|";
    os << IOManipRange(dits.begin(), dits.end(), IOManipRangeOpts{}.set_left("").set_right(""));
  }

  void display_bra_dits_(std::ostream& os, const std::vector<idx>& dits) const {
    os << "<";
    os << IOManipRange(dits.begin(), dits.end(), IOManipRangeOpts{}.set_sep("").set_right(""));
    os << "|";
  }

  void display_mat_dits_(std::ostream& os, const std::vector<idx>& dits, idx n_subsys_rows) const {
    auto split_it = std::next(dits.begin(), n_subsys_rows);
    std::vector<idx> row_dits(dits.begin(), split_it);
    std::vector<idx> col_dits(split_it, dits.end());
    display_ket_dits_(os, row_dits);
    display_bra_dits_(os, col_dits);
  }

  std::ostream& display(std::ostream& os) const override {
    if (A_.states.empty()) {
      return os << "0";
    }

    idx D_rows = prod(A_.dims_rows);
    idx D_cols = prod(A_.dims_cols);

    bool is_scalar = D_rows == 1 && D_cols == 1;
    bool is_col_vec = D_rows > 1 && D_cols == 1;
    bool is_row_vec = D_rows == 1 && D_cols > 1;
    bool is_matrix = D_rows > 1 && D_cols > 1;

    idx n_subsys_rows = A_.dims_rows.size();

    auto display_coeff = [&](Scalar coeff) {
      if (std::abs(std::real(coeff)) > opts_.cplx_opts.chop &&
          std::abs(std::imag(coeff)) > opts_.cplx_opts.chop) {
        os << IOManipScalar<Scalar>{
            coeff, IOManipComplexOpts{opts_.cplx_opts}.set_left("(").set_right(")")};
      } else {
        os << IOManipScalar<Scalar>{coeff,
                                    IOManipComplexOpts{opts_.cplx_opts}.set_left("").set_right("")};
      }
    };

    auto display_dits = [&](const std::vector<idx>& dits) {
      if (is_col_vec) {
        display_ket_dits_(os, dits);
      } else if (is_row_vec) {
        display_bra_dits_(os, dits);
      } else if (is_matrix) {
        display_mat_dits_(os, dits, n_subsys_rows);
      }
    };

    if (is_scalar) {
      display_coeff(A_.states[0].first);
      return os;
    }

    bool first = true;
    for (auto&& elem : A_.states) {
      auto coeff = elem.first;
      auto abs_re = std::abs(std::real(coeff));
      auto abs_im = std::abs(std::imag(coeff));
      if (abs_re < opts_.cplx_opts.chop && abs_im < opts_.cplx_opts.chop && opts_.discard_zeros) {
        continue;
      }

      auto dits = elem.second;
      if (!first) {
        os << opts_.plus_op;
      } else {
        first = false;
      }

      if (opts_.amplitude_after) {
        display_dits(dits);
        os << opts_.mul_op;
        display_coeff(coeff);
      } else {
        display_coeff(coeff);
        os << opts_.mul_op;
        display_dits(dits);
      }
    }

    if (first) {
      os << 0;
    }
    return os;
  }
};

}  // namespace internal
}  // namespace clara

#endif  // !INTERNAL_CLASSFUNCTION_IOMANIP_H_
