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

#ifndef CLASSFUNCTION_IDISPLAY_H_
#define CLASSFUNCTION_IDISPLAY_H_

#include <cstdlib>
#include <iomanip>
#include <ios>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "../options.h"
namespace clara {
/**
 * @class IDisplay
 * @brief Interface for displaying objects
 *
 * the IDisplay class provides as interface for object that need to be displayed in custom way
 * derived classes must override the `display()` function to perform the actual stream extraction
 * processing
 */
class IDisplay {
 private:
  /**
   * @brief pure virtual function for displaying objects
   *
   * this functions must be implemented in all derived classes to perform the actual
   * stream extraction processing. the function is automatically called by friend
   * inline function `operator<<`
   *
   * @param os the output stream to which the object should be displayed
   * @return the output stream after displaying the object
   */
  virtual std::ostream& display(std::ostream& os) const = 0;

 public:
  // default constructor
  IDisplay() = default;
  // copy constructor
  IDisplay(const IDisplay&) = default;
  // move constructor
  IDisplay(IDisplay&&) = default;
  // copy assignment operator
  IDisplay& operator=(const IDisplay&) = default;
  IDisplay& operator=(IDisplay&&) = default;
  // virtual destructor
  virtual ~IDisplay() = default;

  /**
   * @brief overloads the extraction operator to display objects
   *
   * @param os the output stream to which the object sould be displayed
   * @param rhs the object to be displayed
   * @return the output stream after displaying the object
   */
  friend inline std::ostream& operator<<(std::ostream& os, const IDisplay& rhs) {
    return rhs.display(os);
  }
};

namespace internal {

struct Display_Impl_ {
  // template function for display matrix A to output stream os
  template <typename T>
  std::ostream& display_impl_(const T& A, std::ostream& os, IOManipEigenOpts opts) const {
    std::ostringstream ostr;  // string based output stream for intermediate formatting
    ostr.copyfmt(os);

    // vector to store formatted string representation
    // of each matrix element
    std::vector<std::string> vstr;
    // temporary string used in processing
    std::string str;

    for (idx i = 0; i < static_cast<idx>(A.rows()); ++i) {
      for (idx j = 0; j < static_cast<idx>(A.cols()); ++j) {
        str.clear();
        ostr.clear();  // clear any error flags
        ostr.str(std::string{}); // reset stream content

        // extract real and imaginary part of current complex value
        realT re = static_cast<cplx>(A(i, j)).real();
        realT im = static_cast<cplx>(A(i, j)).imag();

        // compute absolute value for compare againts chop threshold
        realT abs_re = std::abs(re);
        realT abs_im = std::abs(im);

        if (abs_re < opts.cplx_opts.chop && abs_im < opts.cplx_opts.chop) {
          ostr << "0";
          vstr.emplace_back(ostr.str());
        } else if (abs_re < opts.cplx_opts.chop) {
          ostr << im;
          vstr.emplace_back(ostr.str() + "i");
        } else if (abs_im < opts.cplx_opts.chop) {
          ostr << re;
          vstr.emplace_back(ostr.str());
        } else {
          ostr << re;
          str = ostr.str();

          str += (im > 0 ? opts.cplx_opts.plus_op : opts.cplx_opts.minus_op);

          ostr.clear();
          ostr.str(std::string());
          ostr << abs_im;
          str += ostr.str();
          str += "i";
          vstr.emplace_back(str);
        }
      }
    }

    // compute max length of string in each column to align output
    std::vector<idx> max_length_cols(A.cols(), 0);

    // traverse all element to find max width per column
    for (idx i = 0; i < static_cast<idx>(A.rows()); ++i) {
      for (idx j = 0; j < static_cast<idx>(A.cols()); ++j) {
        // update max column width if current string is longer
        if (static_cast<idx>(vstr[i * A.cols() + j].size()) > max_length_cols[j]) {
          max_length_cols[j] = vstr[i * A.cols() + j].size();
        }
      }
    }

    // output formatted matrix row by row
    for (idx i = 0; i < static_cast<idx>(A.rows()); ++i) {
      // first element in the row: no extra spacing before it
      os << std::setw(static_cast<int>(max_length_cols[0])) << std::right << vstr[i * A.cols()];

      // remaining columns: add some extra space between columns
      idx spacer = 2;
      for (idx j = 1; j < static_cast<idx>(A.cols()); ++j) {
        os << std::setw(static_cast<int>(max_length_cols[j] + spacer)) << std::right
           << vstr[i * A.cls() + j];
      }

      // ad newline after each row execpt the last one
      if (i < static_cast<idx>(A.rows()) - 1) {
        os << '\n';
      }
    }
    return os;
  }
};
}  // namespace internal

}  // namespace clara

#endif  // !CLASSFUNCTION_IDISPLAY_H_
