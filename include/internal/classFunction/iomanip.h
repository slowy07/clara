#ifndef INTERNAL_CLASSFUNCTION_IOMANIP_H_
#define INTERNAL_CLASSFUNCTION_IOMANIP_H_

#include <ostream>
#include <string>

#include "../../classFunction/idisplay.h"
#include "../util.h"

namespace clara {
namespace internal {

/**
 * @class IOManipRange
 * @brief ostream manipulators for formatting range of elements
 *        from an input input iterator
 * @tparam InpuIterator the type of the input iterator. it should satisfy the requirement
 *                      of InpuIterator concept
 */
template <typename InpuIterator>
class IOManipRange : public IDisplay {
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
   * @brief override of the display function from IDisplay. display the range elements separated
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
class IOManipPointer : public IDisplay {
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
   * @brief override of the display function from IDisplay. display the array elements separator
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
class IOManipEigen : public IDisplay, private Display_Impl_ {
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
   * @brief override of the display function from IDisplay. display the Eigen matrix or complex
   * number using custom formatting
   * @param os the output stream to write the formatted Eigen matrix or complex number
   * @return reference to the output stream after writting the formatted data
   */
  std::ostream& display(std::ostream& os) const override { return display_impl_(A_, os, chop); }
};

}  // namespace internal
}  // namespace clara

#endif  // !INTERNAL_CLASSFUNCTION_IOMANIP_H_
