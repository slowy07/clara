#ifndef INPUT_OUTPUT_H_
#define INPUT_OUTPUT_H_

#include <fstream>
#include <memory>
#include <stdexcept>

#include "classFunction/exception.h"
#include "constants.h"
#include "internal/classFunction/iomanip.h"
#include "internal/util.h"
#include "traits.h"
#include "types.h"

namespace clara {

template <typename Derived>
internal::IOManipEigen disp(const Eigen::MatrixBase<Derived>& A, double chop = clara::chop) {
  return internal::IOManipEigen(A, chop);
}

inline internal::IOManipEigen disp(cplx z, double chop = clara::chop) {
  return internal::IOManipEigen(z, chop);
}

template <typename InputIterator>
internal::IOManipRange<InputIterator> disp(InputIterator first, InputIterator last,
                                           const std::string& separator,
                                           const std::string& start = "[",
                                           const std::string& end = "]") {
  return internal::IOManipRange<InputIterator>(first, last, separator, start, end);
}

template <typename Container>
internal::IOManipRange<typename Container::const_iterator> disp(
    const Container& c, const std::string& separator, const std::string& start = "[",
    const std::string& end = "]",
    typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  return internal::IOManipRange<typename Container::const_iterator>(std::begin(c), std::end(c),
                                                                    separator, start, end);
}

// pointer ostream manipulator
template <typename PointerType>
internal::IOManipPointer<PointerType> disp(const PointerType& p, idx N,
                                           const std::string& separator,
                                           const std::string& start = "[",
                                           const std::string& end = "]") {
  return internal::IOManipPointer<PointerType>(p, N, separator, start, end);
}

// save Eigen expression to a binary file
template <typename Derived>
void save(const Eigen::MatrixBase<Derived>& A, const std::string& fname) {
  const dyn_mat<typename Derived::Scalar>& rA = A.derived();

  // check zero size
  if (!internal::check_nonzero_size(rA))
    throw exception::ZeroSize("clara::save()");
  std::fstream fout;
  fout.open(fname, std::ios::out | std::ios::binary);

  if (fout.fail())
    throw std::runtime_error("clara::save() ERROR writing output file \"" + std::string(fname) +
                             "\"");

  // write header file
  const std::string header_ = "TYPE::Eigen::Matrix";
  fout.write(header_.c_str(), header_.length());

  idx rows = static_cast<idx>(rA.rows());
  idx cols = static_cast<idx>(rA.cols());
  fout.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
  fout.write(reinterpret_cast<const char*>(&cols), sizeof(cols));

  fout.write(reinterpret_cast<const char*>(rA.data()),
             sizeof(typename Derived::Scalar) * rows * cols);
  fout.close();
}

// load eigen matrix from a binary file in double precision
template <typename Derived>
dyn_mat<typename Derived::Scalar> load(const std::string& fname) {
  std::fstream fin;
  fin.open(fname, std::ios::in | std::ios::binary);

  if (fin.fail()) {
    throw std::runtime_error("clara::load() ERROR opening input file \"" + std::string(fname) +
                             "\"");
  }
  const std::string header_ = "TYPE::Eigen::Matrix";
  std::unique_ptr<char[]> fheader_{new char[header_.length()]};

  // read the header from file
  fin.read(fheader_.get(), header_.length());
  if (std::string(fheader_.get(), header_.length()) != header_) {
    throw std::runtime_error("clara::load() corrupted file \"" + std::string(fname) + "\"!");
  }

  idx rows, cols;
  fin.read(reinterpret_cast<char*>(&rows), sizeof(rows));
  fin.read(reinterpret_cast<char*>(&rows), sizeof(cols));
  dyn_mat<typename Derived::Scalar> A(rows, cols);

  fin.read(reinterpret_cast<char*>(A.rows()), sizeof(typename Derived::Scalar) * rows * cols);
  fin.close();
  return A;
}

}  // namespace clara

#endif  // !INPUT_OUTPUT_H_
