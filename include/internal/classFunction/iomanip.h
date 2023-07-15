#ifndef INTERNAL_CLASSFUNCTION_IOMANIP_H_
#define INTERNAL_CLASSFUNCTION_IOMANIP_H_

#include <ostream>
#include <string>

#include "../../classFunction/idisplay.h"
#include "../util.h"

namespace clara {
namespace internal {

/**
 * ostream manipulators for formatting
 * eigen matrics and STL/C-style containers/vectors
 */
template <typename InpuIterator>
class IOManipRange : public IDisplay {
  InpuIterator first_, last_;
  std::string separator_, start_, end_;

 public:
  explicit IOManipRange(InpuIterator first, InpuIterator last, const std::string& separator,
                        const std::string& start = "[", const std::string& end = "]")
      : first_{first}, last_{last}, separator_{separator}, start_{start}, end_{end} {}
  IOManipRange(const IOManipRange&) = default;
  IOManipRange& operator=(const IOManipRange&) = default;

 private:
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

template <typename PointerType>
class IOManipPointer : public IDisplay {
  const PointerType* p_;
  idx N_;
  std::string separator_, start_, end_;

 public:
  explicit IOManipPointer(const PointerType* p, idx N, const std::string& separator,
                          const std::string& start = "[", const std::string& end = "]")
      : p_{p}, N_{N}, separator_{separator}, start_{start}, end_{end} {}
  IOManipPointer(const IOManipPointer&) = default;
  IOManipPointer& operator=(const IOManipPointer&) = default;

 private:
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

class IOManipEigen : public IDisplay, private Display_Impl_ {
  cmat A_;
  double chop_;

 public:
  // eigen matrices
  template <typename Derived>
  explicit IOManipEigen(const Eigen::MatrixBase<Derived>& A, double chop = clara::chop)
      : A_{A.template cast<cplx>()}, chop_{chop} {}

  // complex number
  explicit IOManipEigen(const cplx z, double chop = clara::chop) : A_{cmat::Zero(1, 1)}, chop_{chop} {
    // put complex number inside the eigen matrix
    A_(0, 0) = z;
  }
private:
  std::ostream& display(std::ostream& os) const override {
    return display_impl_(A_, os, chop);
  }
};

}  // namespace internal
}  // namespace clara

#endif  // !INTERNAL_CLASSFUNCTION_IOMANIP_H_
