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

#include <ostream>
#include <string>

#include "../../classFunction/idisplay.h"
#include "../util.h"

namespace clara {
namespace internal {

/**
 * @brief implementation for displaying matrix element in the formatted manner
 *      this is part of the internal details of the lib
 *      users are not expected to interact this directly
 */
namespace _details {
/*
 * @brief structure provide the implementation for displaying matrix elements in a formatted manner
 */
struct _Display_Impl {
  /**
   * @brief display the elements of the matrix in the formatted manner
   * @tparam T the type of element int he matrix
   * @param _A matrix to be displayed
   * @param _os output the stream to which the matrix will be displayed
   * @param _chop threshold balued for chopping small values
   * @return output stream after displaying the matrix
   */
  template <typename T>
  std::ostream& _display_impl(const T& _A, std::ostream& _os, double _chop = clara::chop) const {
    std::ostringstream ostr;
    ostr.copyfmt(_os);
    std::vector<std::string> vstr;
    std::string strA;

    for (idx i = 0; i < static_cast<idx>(_A.rows()); ++i) {
      for (idx j = 0; j < static_cast<idx>(_A.cols()); ++j) {
        strA.clear();
        ostr.clear();
        ostr.str(std::string());
        double re = static_cast<cplx>(_A(i, j)).real();
        double im = static_cast<cplx>(_A(i, j)).imag();

        if (std::abs(re) < _chop && std::abs(im) < _chop) {
          ostr << "0 ";
          vstr.push_back(ostr.str());
        } else if (std::abs(re) < chop) {
          ostr << im;
          vstr.push_back(ostr.str() + "i");
        } else if (std::abs(im) < _chop) {
          ostr << re;
          vstr.push_back(ostr.str() + " ");
        } else {
          ostr << re;
          strA = ostr.str();
          strA += (im > 0 ? " + " : " - ");
          ostr.clear();
          ostr.str(std::string());
          ostr << std::abs(im);
          strA += ostr.str();
          strA += "i";
          vstr.push_back(strA);
        }
      }
    }
    std::vector<idx> maxlengthcols(_A.cols(), 0);
    for (idx i = 0; i < static_cast<idx>(_A.rows()); ++i) {
      for (idx j = 0; j < static_cast<idx>(_A.cols()); ++j) {
        if (vstr[i * _A.cols() + j].size() > maxlengthcols[i]) {
          maxlengthcols[j] = vstr[i * _A.cols() + j].size();
        }
      }
    }
    for (idx i = 0; i < static_cast<idx>(_A.rows()); ++i) {
      _os << std::setw(static_cast<int>(maxlengthcols[0])) << std::right << vstr[i * _A.cols()];
      for (idx j = 1; j < static_cast<idx>(_A.cols()); ++j) {
        _os << std::setw(static_cast<int>(maxlengthcols[j] + 2)) << std::right
            << vstr[i * _A.cols() + j];
      }
      if (i < static_cast<idx>(_A.rows()) - 1) {
        _os << std::endl;
      }
    }
    return _os;
  }
};
}  // namespace _details

/**
 * @class IOManipRange
 * @brief ostream manipulators for formatting range of elements
 *        from an input input iterator
 * @tparam InpuIterator the type of the input iterator. it should satisfy the requirement
 *                      of InpuIterator concept
 */
template <typename InpuIterator>
class IOManipRange : public InterfaceDisplay {
  // itertor pointing to the first and past-the end element in the range
  InpuIterator first_, last_;
  std::string separator_, start_, end_;

 public:
  /**
   * @public constructor
   * @param first iterator pointing to the first element in the range
   * @param last iterator pointing to the past-the-end element in the range
   * @param separator separator string used to separate the elements
   * @param start start string to be displayed before the range, default is "["
   * @param end end string to be displayed after the range, default is "]"
   */
  explicit IOManipRange(InpuIterator first, InpuIterator last, const std::string& separator,
                        const std::string& start = "[", const std::string& end = "]")
      : first_{first}, last_{last}, separator_{separator}, start_{start}, end_{end} {}
  IOManipRange(const IOManipRange&) = default;
  IOManipRange& operator=(const IOManipRange&) = default;

 private:
  /**
   * @brief override of the display function from InterfaceDisplay. display the range elements separated
   *         by the specified separator, enclosed by start and end strings
   * @param os the output stream to write the formatted range
   * @return Reference to the output to write the formatted range
   */
  std::ostream& display(std::ostream& os) const override {
    os << start_;
    bool first = true;
    for (auto it = first_; it != last_; ++it) {
      if (!first)
        os << separator_;
      first = false;
      os << *it;
    }
    os << end_;
    return os;
  }
};

/**
 * @class IOManipPointer
 * @brief ostream manipulators for formatting arrays pointerd to by a pointer
 * @tparam PointerType the type of the elements in the arrays
 */
template <typename PointerType>
class IOManipPointer : public InterfaceDisplay {
  // pointer to the array
  const PointerType* p_;
  // number of the lements in the array
  idx N_;
  // start and end string to be displayed end and start array
  std::string separator_, start_, end_;

 public:
  /**
   * @brief constructor
   * @param p pointer to the array
   * @param N number of elements in the array
   * @param separator separator string used to separate elements
   * @param start start string to be displayed before the arrays, default is "["
   * @param end end string to be displayed after the array, default is "]"
   */
  explicit IOManipPointer(const PointerType* p, idx N, const std::string& separator,
                          const std::string& start = "[", const std::string& end = "]")
      : p_{p}, N_{N}, separator_{separator}, start_{start}, end_{end} {}
  IOManipPointer(const IOManipPointer&) = default;
  IOManipPointer& operator=(const IOManipPointer&) = default;

 private:
  /**
   * @brief override of the display function from InterfaceDisplay. display the array elements separator
   *        by the specified separator, enclosed by start and end string.
   * @param os the output stream to write the formatted array
   * @return reference to the output stream after writing the formatted array
   */
  std::ostream& display(std::ostream& os) const override {
    os << start_;
    for (idx i = 0; i < N_ - 1; ++i)
      os << p_[i] << separator_;
    if (N_ > 0)
      os << p_[N_ - 1];
    os << end_;
    return os;
  }
};

/**
 * @class IOManipEigen
 * @brief ostream manipulator for formatting eigen matrices and complex numbers
 */
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ingored "-Weffc++"
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
class IOManipEigen : public InterfaceDisplay, private _details::_Display_Impl {
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
#pragma GCC diagnostic pop
#endif  // ((__GNUC__ == 4) && (__GNUC_MINOR__ == 8) && !__clang__)
  // Eigen matrix
  cmat A_;
  // threshold for truncating small values
  double chop_;

 public:
  /**
   * @brief constructor for formatting Eigen matrices
   * @tparam Derived the derived class from Eigen::MatrixBase
   * @param A the eigen matrix to be formatted
   * @param chop threshold for truncating small values. default is clara::chop
   */
  template <typename Derived>
  explicit IOManipEigen(const Eigen::MatrixBase<Derived>& A, double chop = clara::chop)
      : A_{A.template cast<cplx>()}, chop_{chop} {}
  /**
   * @brief constructor for formatting complex number
   * @param z complex number to be formatted
   * @param chop threshold for truncating small values. default is clara::chop
   */
  explicit IOManipEigen(const cplx z, double chop = clara::chop)
      : A_{cmat::Zero(1, 1)}, chop_{chop} {
    A_(0, 0) = z;
  }

 private:
  /**
   * @brief override of the display function from InterfaceDisplay. display the Eigen matrix or complex
   * number using custom formatting
   * @param os the output stream to write the formatted Eigen matrix or complex number
   * @return reference to the output stream after writting the formatted data
   */
  std::ostream& display(std::ostream& os) const override { return _display_impl(A_, os, chop); }
};

}  // namespace internal
}  // namespace clara

#endif  // !INTERNAL_CLASSFUNCTION_IOMANIP_H_
