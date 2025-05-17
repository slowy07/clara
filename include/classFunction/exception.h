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

#ifndef CLASSFUNCTION_EXCEPTION_H_
#define CLASSFUNCTION_EXCEPTION_H_

#include <exception>
#include <optional>
#include <string>
namespace clara {

namespace exception {

/**
 * @brief base exception class that all other exceptions in the library
 *
 * this class provide common interface for throwing descriptive exceptions, including context
 * information if available
 */
class Exception : public std::exception {
 protected:
  std::string where_;                   // exception ocurred
  mutable std::string message_;         // cached error message return
  std::optional<std::string> context_;  // optional additional context for debug

 public:
  /**
   * @brief construct a new exception object
   *
   * @param where description of where the exception occured
   * @param context optional additional context or explanation about the error
   */
  explicit Exception(std::string where, std::optional<std::string> context = std::nullopt)
      : where_{std::move(where)}, message_{}, context_{std::move(context)} {}

  /**
   * @brief return the error message a C-sytle string
   *
   * this method override the standard what() method from std::exception
   * it constructs the full message on demand using the `where` location,
   * description(), and optional context
   */
  const char* what() const noexcept override {
    message_.clear();                    // rebuild message every time
    message_ += "[! : " + where_ + ']';  // prefix with location
    message_ += " " + description();     // append specific description

    if (context_.has_value()) {
      message_ += "[MESSAGE: " + context_.value() + ']';
    }
    return message_.c_str();
  }

  /**
   * @brief provide human readable description of the exception type
   *
   * must be overriden by derived classes to provided specific error message
   */
  virtual std::string description() const = 0;
};

// providing default implementation for the base class pure virtual function
inline std::string Exception::description() const { return "clara::exception::Exception"; }

/**
 * @brief generic unknown exception
 */
class Unknown : public Exception {
 public:
  std::string description() const override { return "Unknown exception"; }
  using Exception::Exception;
};

/**
 * @brief throw when an object has zero size, which may not be allowed
 */
class ZeroSize : public Exception {
 public:
  std::string description() const override { return "object has zero size"; }
  using Exception::Exception;
};

/**
 * @brief throw when a matrix is expected to be squared by isn't
 */
class MatrixNotSquare : public Exception {
 public:
  std::string description() const override { return "matrix is not square"; }
  using Exception::Exception;
};

/**
 * @brief throw when a matrix is expected to be a column vector but isn't
 */
class MatrixNotCvector : public Exception {
 public:
  std::string description() const override { return "matrix is not a column vector"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a matrix is expected to be a row vector but isn't
 */
class MatrixNotRvector : public Exception {
 public:
  std::string description() const override { return "matrix is not a row vector"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a matrix must be either square or a column vector
 */
class MatrixNotSquareNorCvector : public Exception {
 public:
  std::string description() const override { return "matrix is not square nor column vector"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a matrix must be eiter square or a row vector
 */
class MatrixNotSquareNorRvector : public Exception {
 public:
  std::string description() const override { return "matrix is not square nor row vector"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a matrix does not match a subsystem requirement
 */
class MatrixMismatchSubsys : public Exception {
 public:
  std::string description() const override { return "matrix mismatch subsystem"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a dimension value is invalid
 */
class DimsInvalid : public Exception {
 public:
  std::string description() const override { return "invalid dimension"; }
  using Exception::Exception;
};

class DimsNotEqual : public Exception {
 public:
  std::string description() const override { return "dimension not equal"; }
  using Exception::Exception;
};

/**
 * @brief throw when dimension do not match the size of a matrix
 */
class DimsMismatchMatrix : public Exception {
 public:
  std::string description() const override { return "dimension mismatch matrix size"; }
  using Exception::Exception;
};

/**
 * @brief thrown when dimension do not match the size of a column vector
 */
class DimsMismatchCvector : public Exception {
 public:
  std::string description() const override { return "dimension mismatch column vector size"; }
  using Exception::Exception;
};

/**
 * @brief throw when dimension do not match the size of a row vector
 */
class DimsMismatchRvector : public Exception {
 public:
  std::string description() const override { return "dimension mismatch row vector size"; }
  using Exception::Exception;
};

/**
 * @brief throw when dimension do not match the size of a general vector
 */
class DimsMismatchVector : public Exception {
 public:
  std::string description() const override { return "dimension mismatch vector size"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a subsystem doesn't match the expected dimension
 */
class SubsysMismatchDims : public Exception {
 public:
  std::string description() const override { return "subsystem mismatch dimension"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a permutation is invalid or malformed
 */
class PermInvalid : public Exception {
 public:
  std::string description() const override { return "invalid permutation"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a permutation doesn't match the expected dimension
 */
class PermMismatchDims : public Exception {
 public:
  std::string description() const override { return "permutation mismatch dimension"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a matrix is expected to be 2x2
 */
class NotQubitMatrix : public Exception {
 public:
  std::string description() const override { return "matrix is not 2 x 2"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a column vector is expected to be 2x1
 */
class NotQubitCvector : public Exception {
 public:
  std::string description() const override { return "column vector is not 2 x 1"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a row vector is expected to be 1x2
 */
class NotQubitRvector : public Exception {
 public:
  std::string description() const override { return "row vector is not 1 x 2"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a vector is expected to be either 2x1
 */
class NotQubitVector : public Exception {
 public:
  std::string description() const override { return "vector is not 2 x 1 nor 1 x 2"; }
  using Exception::Exception;
};

/**
 * @brief throw when a subsystem is expected to contain only qubits, but doesn't
 */
class NotQubitSubsys : public Exception {
 public:
  std::string description() const override { return "subystem are not qubits"; }
  using Exception::Exception;
};

/**
 * @brief throw when a system is expected to be bi-partite but isn't
 */
class NotBipartite : public Exception {
 public:
  std::string description() const override { return "not bi-partite"; }
  using Exception::Exception;
};

/**
 * @brief throw when a codewrd does not exist in a lookup or database
 */
class NoCodeword : public Exception {
 public:
  std::string description() const override { return "codeword does not exist"; }
  using Exception::Exception;
};

/**
 * @brief thrown when an argument is out of its valid range
 */
class OutOfRange : public Exception {
 public:
  std::string description() const override { return "argument out of range"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a requested element cannot be found
 */
class NotFound : public Exception {
 public:
  std::string description() const override { return "element not found"; }
  using Exception::Exception;
};

/**
 * @thrown when there is a mismatch in data types
 */
class TypeMismatch : public Exception {
 public:
  std::string description() const override { return "type mismatch"; }
  using Exception::Exception;
};

/**
 * @brief thrown when size of two objects are expected to match but don't
 */
class SizeMismatch : public Exception {
 public:
  std::string description() const override { return "size mismatch"; }
  using Exception::Exception;
};

/**
 * @thrown when an operation is not defined for a given type
 */
class UndefinedType : public Exception {
 public:
  std::string description() const override { return "not defined for this type"; }
  using Exception::Exception;
};

/**
 * @brief throw when attempting to re-measure a qubit that was already measured
 */
class QuditAlreadyMeasured : public Exception {
 public:
  std::string description() const override { return "qubit was already measured"; }
  using Exception::Exception;
};

/**
 * @brief throw when duplicate elements are detected in a system
 */
class Duplicates : public Exception {
 public:
  std::string description() const override { return "system contains duplicate"; }
  using Exception::Exception;
};

/**
 * @brief custom exception placeholder for user-defined errors
 */
class CustomException : public Exception {
 public:
  std::string description() const override { return "custom exception"; }
  using Exception::Exception;
};

/**
 * @brief thrown when a feature function is not yet implemented
 */
class NotImplemented : public Exception {
 public:
  std::string description() const override { return "not implemented"; }
  using Exception::Exception;
};

/**
 * @brief thrown when iterator is in an ivalid state or used incorrectly
 */
class InvalidIterator : public Exception {
 public:
  std::string description() const override { return "invalid iterator"; }
  using Exception::Exception;
};

/**
 * @brief throw when conditional expression or statement is invalid
 */
class InvalidConditional : public Exception {
 public:
  std::string description() const override { return "invalid conditional"; }
  using Exception::Exception;
};

/**
 * @brief throw when a circuit contains conditional but shouldn't
 */
class CircuitContainsConditionals : public Exception {
 public:
  std::string description() const override { return "circuit contains conditionals"; }
  using Exception::Exception;
};

}  // namespace exception
}  // namespace clara

#endif  // !CLASSFUNCTION_EXCEPTION_H_
