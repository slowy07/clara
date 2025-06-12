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

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <string>
#include <utility>

#include "types.h"

namespace clara {

/**
 * @brief double chop used in clara::disp() for setting to zero numbers
 * that have their absolute value smaller than clara:chop
 */
constexpr realT chop = 1e-10;

/**
 * @brief configuration option for format scalar value during I/O operation
 *
 * this struct allow fine-gained control over how scalar value are displayed
 */
struct IOManipScalarOpts {
  std::string left = "";
  std::string right = "";
  realT chop = clara::chop;

  /**
   * @brief set the left delimiter string
   *
   * @param left_delimiter string to preprend before the scalar value
   * @return reference to this object for method chaining
   */
  IOManipScalarOpts& set_left(std::string left_delimiter) {
    left = std::move(left_delimiter);
    return *this;
  }

  /**
   * @brief set the right delimiter string
   *
   * @param right_delimiter string to append after sclar value
   * @return reference to this object for method chaining
   */
  IOManipScalarOpts& set_right(std::string right_delimiter) {
    right = std::move(right_delimiter);
    return *this;
  }

  /**
   * @brief set chopping threshold for near-zero values
   *
   * @param chop_at threshold below which value are displayed as zero
   * @return reference to this object for method chaining
   */
  IOManipScalarOpts& set_chop(realT chop_at) {
    chop = chop_at;
    return *this;
  }
};

/**
 * @brief configuration options for formatting complex number using I/O operation
 *
 * struct will allow full control over how complex number are displayed
 * making it useful in context like matrix printing, quantum state visualization.
 */
struct IOManipComplexOpts {
  std::string im_suffix = "i";
  std::string plus_op = " + ";
  std::string minus_op = " - ";
  std::string left = "";
  std::string right = "";
  realT chop = clara::chop;

  /**
   * @brief set imaginary unit suffix
   *
   * @param imaginary_suffix new suffix string for imaginary components
   * @return reference to this object for method chaining
   */
  IOManipComplexOpts& set_im_suffix(std::string imaginary_suffix) {
    im_suffix = std::move(imaginary_suffix);
    return *this;
  }

  /**
   * @brief set the string used for addition between real and imaginary parts
   *
   * @param plus_operator string to use for positive imaginary parts
   * @return reference to this object for method chaining
   */
  IOManipComplexOpts& set_plus_operation(std::string plus_operator) {
    plus_op = std::move(plus_operator);
    return *this;
  }

  /**
   * @brief string used for subtract between real and imaginary part
   *
   * @param minus_operator string to use for negative imaginary part
   * @return reference to object for method chaining
   */
  IOManipComplexOpts& set_minus_operation(std::string minus_operator) {
    minus_op = std::move(minus_operator);
    return *this;
  }

  /**
   * @brief set left delimiter string used in formatted output
   *
   * delimiter typically using at the beginning of matrix vector representation
   *
   * @param left_delimiter new left delimiter string
   * @return reference to this object for method chaining
   */
  IOManipComplexOpts& set_left(std::string left_delimiter) {
    left = std::move(left_delimiter);  // move to avoid unnecessary copy
    return *this;
  }

  /**
   * @brief set the right delimiter string used in formatted output
   *
   * delimiter typically used at the end of a matrix or vector representation
   *
   * @param right_delimiter new right delimiter string
   * @return reference to this object for method chaining
   */
  IOManipComplexOpts& set_right(std::string right_delimiter) {
    right = std::move(right_delimiter);  // move to avoid unnecessary copy
    return *this;
  }

  /**
   * @brief set the choping threshold for near-zero values
   *
   * value with absolute magnitude less than this threshold will be display as zero
   * this help reduce floating-point noise in printed output
   *
   * @param chop_at threshold value below wich number are considered zero
   * @return reference to this object for method chaining
   */
  IOManipComplexOpts& set_chop(realT chop_at) {
    chop = chop_at;
    return *this;
  }
};

/**
 * @brief options class for customizing eigen-like matrix output behaviour
 *
 * this class aggregate formatting optons for use when displaying matrices,
 * particularly complex-valued one, it encapsulate complex-specific formatting
 * through an instance of `IOManipComplexOpts`
 */
struct IOManipEigenOpts {
  IOManipComplexOpts cplx_opts{};

  /**
   * @brief set the full set of complex number formatting options
   *
   * allow replacing formatting settings with a new `IOManipComplexOpts`
   * object, enabling fine-grained control over how complex value are displayed
   *
   * @param complex_opts new complex formatting options to apply
   * @return reference to this object for method chaining
   */
  IOManipEigenOpts& set_complex_opts(IOManipComplexOpts complex_opts) {
    cplx_opts = std::move(complex_opts);  // move assignment for efficiency
    return *this;
  }
};

/**
 * @brief option class for customizing display format of range-based containers
 *
 * this struct provide customizable formatting options for printing sequence
 * like vector
 */
struct IOManipRangeOpts {
  std::string sep = " ";
  std::string left = "[";
  std::string right = "]";
  realT chop = clara::chop;

  /**
   * @brief set separator string used between element in formatted output
   *
   * @param separator new separator string to insert between element
   * @return reference to this object for method chaining
   */
  IOManipRangeOpts& set_sep(std::string separator) {
    sep = std::move(separator);
    return *this;
  }

  /**
   * @brief set left delimiter string used to enclose range output
   *
   * @param left_delimiter new left delimiter string
   * @return reference to this object for method chaining
   */
  IOManipRangeOpts& set_left(std::string left_delimiter) {
    left = std::move(left_delimiter);
    return *this;
  }

  /**
   * @brief set right delimiter string used to enclose range output
   *
   * @param right_delimiter new right delimiter string
   * @return reference to this object for method chaining
   */
  IOManipRangeOpts& set_right(std::string right_delimiter) {
    right = std::move(right_delimiter);
    return *this;
  }

  /**
   * @brief set the chopping threshold for near-zero value
   *
   * value whose absolute magnitude is less than this threshold will be shown
   * as `0`, this help reduce visual clutter caused by floating-point imprcision
   *
   * @param chop_at threshold value below which number are treated as zero
   * @return reference to this object for method chaining
   */
  IOManipRangeOpts& set_chop(realT chop_at) {
    chop = chop_at;
    return *this;
  }
};

/**
 * @brief options class for customizing output of general containers
 *
 * this struct provide flexible formatting option for printing container like vector
 * it can be implicity convert to `IOManipRangeOpts` when needed for compatibility
 * with range-based display function
 */
struct IOManipContainerOpts {
  std::string sep = " ";
  std::string left = "[";
  std::string right = "]";
  realT chop = clara::chop;

  /**
   * @brief set the separator string used betwen element in formatted output
   *
   * @param separator new separator string to insert between elements
   * @return refernece to this object for method chaining
   */
  IOManipContainerOpts& set_sep(std::string separator) {
    sep = std::move(separator);
    return *this;
  }

  /**
   * @brief set the left dlimiter string used to enclose the container output
   *
   * @param left_delimiter new left delimiter string
   * @return reference to this object for method chainig
   */
  IOManipContainerOpts& set_left(std::string left_delimiter) {
    left = std::move(left_delimiter);
    return *this;
  }

  /**
   * set right delimiter string used to enclose the container output
   *
   * @Param right_delimiter new right delimiter string
   * @return reference to this object for method chaining
   */
  IOManipContainerOpts& set_right(std::string right_delimiter) {
    right = std::move(right_delimiter);
    return *this;
  }

  /**
   * @brief set the choping threshold for near-zero values
   *
   * @param chop_at threshold value below which number are treated as zero
   * @return reference to this object for method chaining
   */
  IOManipContainerOpts& set_chop(realT chop_at) {
    chop = chop_at;
    return *this;
  }

  /**
   * @brief implicit conversion operator
   *
   * @return an `IOManipRangeOpts` object initialized with current setting
   */
  operator IOManipRangeOpts() {
    IOManipRangeOpts range_opts;
    range_opts.sep = sep;
    range_opts.left = left;
    range_opts.right = right;
    range_opts.chop = chop;
    return range_opts;
  }
};

/**
 * @brief options class for customizing output of pointer-based container or array
 *
 * similar to `IOManipContainerOpts` byt intended for low-level data structure
 * like raw pointer or memory buffer that representation sequence of values
 */
struct IOManipPointerOpts {
  std::string sep = " ";
  std::string left = "[";
  std::string right = "]";
  realT chop = clara::chop;

  /**
   * @brief set the separator string used between element in formatted output
   *
   * @param separator new separator string to insert between element
   * @return reference to this object for method chaining
   */
  IOManipPointerOpts& set_sep(std::string separator) {
    sep = std::move(separator);
    return *this;
  }

  /**
   * @brief set left delimiter string used to enclose the container output
   *
   * @param left_delimiter new left delimiter string
   * @return reference to this object for method chaining
   */
  IOManipPointerOpts& set_left(std::string left_delimiter) {
    left = std::move(left_delimiter);
    return *this;
  }

  /**
   * @brief set the right delimiter string used to enclose the container output
   *
   * @param right_delimiter new right delimiter string
   * @return reference to this object for method chaining
   */
  IOManipPointerOpts& set_right(std::string right_delimiter) {
    right = std::move(right_delimiter);
    return *this;
  }

  /**
   * @brief set the chopping threshold for near-zero value
   *
   * value whose absolute magnitude is less thatn this threshold will be shown as
   * `0`, this helps reduce visual clutter caused by floating-point imprecision
   *
   * @param chop_at threshold value below which number are treated as zero
   * @return reference to this object for method chaining
   */
  IOManipPointerOpts& set_chop(realT chop_at) {
    chop = chop_at;
    return *this;
  }
};

/**
 * @brief options class for customizing output format of quantum state in dirac
 *        notation
 *
 * provide fine-grained control over how quantum state are displayed
 * particularly useful in quantum computing libraries where bra-ket notation is common
 */
struct IOManipDiracOpts {
  IOManipComplexOpts cplx_opts{};

  std::string plus_op = " +\n";
  std::string mul_op = " * ";
  bool amplitude_after = true;
  bool discard_zeros = true;

  /*
   * @brief set full set complex number formatting options
   *
   * allow replacing current complex formatting setting with a new `IOManipComplexOpts`
   * object, enabling fine-grained constrol over how amplitude are displayed
   *
   * @param complex_opts new complex formatting options to apply
   * @return reference to this object for method chaining
   */
  IOManipDiracOpts& set_complex_opts(IOManipComplexOpts complex_opts) {
    cplx_opts = std::move(complex_opts);
    return *this;
  }

  /**
   * @brief set the operator string used between different terms in dirac notation
   *
   * useful for changing how components are sparated visually
   * using commas, semicolon, or even no separator at all
   *
   * @param plus_operator new operator string to insert between terms
   * @return reference to this object for method chaining
   */
  IOManipDiracOpts& set_plus_op(std::string plus_operator) {
    plus_op = std::move(plus_operator);
    return *this;
  }

  /**
   * @brief set operator string used between amplitude and basis state
   *
   * this affect how the amplited is combined with the basis vector
   *
   * @param mul_operator new separator string to use between amplitude basis
   * @return refernece to this object for method chaining
   */
  IOManipDiracOpts& set_mul_op(std::string mul_operator) {
    mul_op = std::move(mul_operator);
    return *this;
  }

  /**
   * @brief control the order of amplitude and basis state in output
   *
   * @param show_amplitudes_after whether to print amplitude after basis state
   * @return reference to this object for method chaining
   */
  IOManipDiracOpts& set_amplitudes_after(bool show_amplitudes_after) {
    amplitude_after = show_amplitudes_after;
    return *this;
  }

  /**
   * @brief enable or disable printing of zero-amplitude component
   *
   * when enabled (default), components with amplitude below the chop threshold
   * are omitted entirely from output
   *
   * @param discarding_zeros whether to skip zero-amplitude term
   * @return reference to this object for method chaining
   */
  IOManipDiracOpts& set_discard_zeros(bool discarding_zeros) {
    discard_zeros = discarding_zeros;
    return *this;
  }
};

namespace internal {
/**
 * @brief constant defining the max number of terms allowed in dirac notation output
 *
 * using limit the output of basis states shown when printing quantum state
 * to avoid overwhelming console output, can be adjust
 */
constexpr idx maxn = 64;
}  // namespace internal

}  // namespace clara

#endif  // !OPTIONS_H_
