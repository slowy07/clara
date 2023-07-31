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

/**
 * @brief create IOManipEigen object for eigen matrices
 * @tparam derived the matrix expresion type
 * @param A the eigen matrix or matrix expression to displayed
 * @param chop the precssion for displaying floating-point numbers default is clara::chop
 * @return internal::IOManipEigen an IOManip object specifically tailored for eigen matrix
 *
 * @example
 * Eigen::MatrixXd mat(3, 3);
 * // initialize matrix `mat` with some values
 *
 * // display the matrix with default precision
 * disp(mat) << std::endl;
 *
 * // display the matrix with custom precision
 *disp(mat, 0.001) << std::endl;
 */
template <typename Derived>
internal::IOManipEigen disp(const Eigen::MatrixBase<Derived>& A, double chop = clara::chop) {
  return internal::IOManipEigen(A, chop);
}

/**
 * @brief create IOManipEigen object for complex numbers
 * @param z the complex number to be displayed
 * @param chop the precision for displyaing floating-point numbers. defeault is clara::chop
 * @return internal::IOManipEigen an IOManip object specifically tailored for complex matrices
 *
 * @example
 * // make complex number
 * std::complex<double> complexNum(3.14, 2.17);
 *
 * // display the complex number with default precision
 * disp(complexNum) << std::endl;
 *
 * // display the complex number with custom precision
 * disp(complexNum, 0.0001) << std::endl;
 * '
 */
inline internal::IOManipEigen disp(cplx z, double chop = clara::chop) {
  return internal::IOManipEigen(z, chop);
}

/**
 * @brief function to create IOManipRange object for range of elements
 * @tparam INputIterator the input iterator type for the range
 * @param first an iterator pointing to the first element of the range
 * @param last an iterator pointing to the last element of the rangee
 * #param separator the string to be displayed to be used as separator between elements during
 * display
 * @param start the string to be display before the range of elements. default "["
 * @param end the string to be display after the range of elements default "]"
 *
 * @eaxmple
 * std::vector<int> vec (1, 2, 3, 4);
 *
 * // display the vector elements with default formatting
 * disp(vec.begin, vec.end(), ", ") << std::endl;
 *
 * // display the vector element with custom formatting
 * disp(vec.begin(), vec.end(), ";", "<", ">") << std::endl;
 */
template <typename InputIterator>
internal::IOManipRange<InputIterator> disp(InputIterator first, InputIterator last,
                                           const std::string& separator,
                                           const std::string& start = "[",
                                           const std::string& end = "]") {
  return internal::IOManipRange<InputIterator>(first, last, separator, start, end);
}

/**
 * @brief create IOManipRange for a range of elements defined by two iterator
 * @tparam InputIterator the input iterator type for the the range
 * @param first an iterator pointing to the first element of the range
 * @param last an iterator pointing to the last element of the range
 * @param separator the string to be used as separator between elements during display
 * @param start the string to be displayed before the range of elements. default value is "["
 * @param end the string to be displyaed after the range of elements, default value is "]"
 * @return internal::IOManipRange<InputIterator> an IOManip object speciefied tailored of elements
 *
 * @example
 * std::vector<int> vec = {1 ,2 ,3 ,4 ,5};
 *
 * // diplaye the vector elements with default formatting
 * disp(vec.begin, vec.end(), ", ") << std::endl;
 *
 * // display the vector elements with custom foramtting
 * disp(vec.begin(), vec.edn(), "; ", "<", ">") std::endl;
 * // displaythe set elements with custom formatting
 * disp(s, "; ", "{", "}") << std::endl
 */

template <typename Container>
internal::IOManipRange<typename Container::const_iterator> disp(
    const Container& c, const std::string& separator, const std::string& start = "[",
    const std::string& end = "]",
    typename std::enable_if<is_iterable<Container>::value>::type* = nullptr) {
  return internal::IOManipRange<typename Container::const_iterator>(std::begin(c), std::end(c),
                                                                    separator, start, end);
}

/**
 * @brief create an IOManipRange object for a range of elements defined by two iteratios
 *
 * @tparam InputIterator the input iterator type for the range
 * @param first an iterator pointing to the first element of the range
 * @param last an iterator pointing to the last element of the range (on past at end)
 * @param separator the tring to be used as an separator between element during display
 * @param start the string to be displayed before the range of elements, default is "["
 * @param end string to be displayed after the range of elements, default is "]"
 * @return internal::IOManipRange<InputIterator> an IOManip object specifically tailored
 *                               of ranges of elements
 *
 * @example
 * std::vector<vec> vec = {1, 2, 3, 4, 5};
 *
 * // display the vector elemens with default formatting
 * disp(vec.begin(), vec.end(), ", ") << std::endl;
 */
template <typename PointerType>
internal::IOManipPointer<PointerType> disp(const PointerType& p, idx N,
                                           const std::string& separator,
                                           const std::string& start = "[",
                                           const std::string& end = "]") {
  return internal::IOManipPointer<PointerType>(p, N, separator, start, end);
}

/**
 * @brief save an Eigen matrix or matrix expression to a binary file
 * @tparam derived matrix or matrix expression to be saved
 * @param fname the file name (including the path) where the matrix will be saved
 *
 * @example
 * Eigen::MatrixXd mat(3, 3);
 *
 * // save the matrix to a binary file "example matrix.bin"
 * save(mat, "matrix.bin")
 */
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

/**
 * @brief load an eigen matrix or eigen expression from a binary file
 * @tparam derived the matrix expression type
 * @param fname file name (including the path) from where the matrix will be loaded
 * @return dyn_mat<typename Derived::Scalar> the loaded eigen matrix
 *
 * @example
 * Eigen::MatrixXd loadedMatrix = load <Eigen::MatrixXd>("matrix.bin");
 */
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
