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
#include <string>
namespace clara {

namespace exception {

/**
 * @class Exception
 * @brief base class for custom exception handling in clara library
 *
 * this class derived from std::exception and servers as the base class for custom exxception
 * handling in the clara library. it provides functionality to construct informative error messages
 * that include the type of exception and the location where the exception occured
 */
class Exception : public std::exception {
 private:
  std::string where_;
  mutable std::string msg_;

 public:
  /**
   * @brief construct for the exception class
   * @param where the location where the exception osccured (usually the function name or class
   * name)
   */
  Exception(const std::string& where) : where_{where}, msg_{} {}

  /**
   * @brief function to get the exception messages
   * @return the exception message as an C-style string
   */
  virtual const char* what() const noexcept override {
    msg_.clear();
    msg_ += "IN ";
    msg_ += where_;
    msg_ += ": ";
    msg_ += this->type_description();
    msg_ += "!";
    return msg_.c_str();
  }
  /**
   * @brief virtual function to get a description of the specific type of exception
   * @return the type description of the exception as string
   */
  virtual std::string type_description() const = 0;
};

inline std::string Exception::type_description() const { return "clara::exception::Exception"; }

/**
 * @brief Unknown
 * @brief Exception throw when no other exception is suitable
 *
 * this class is derived from exception and represent an uknwon exception
 * it is thrown when no ther specific type is suitable for a particular error
 */
class Unknown : public Exception {
 public:
  std::string type_description() const override { return "UNKNOWN EXCEPTION"; }
  /**
   * @brief constructur for the unknown class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class ZeroSize
 * @brief Exception representing an object with zero size
 *
 * this calss is derived from Exception and represents an exception
 * where an object has zero size, which is invalid in the context of the operation
 */
class ZeroSize : public Exception {
 public:
  /**
   * @brief get the type description of the ZeroSize exception
   * @return type description of the ZeroSize exception as a string
   */
  std::string type_description() const override { return "Object has zero size"; }
  /**
   * @brief for the ZeroSize class
   * @param where the location where the exception occurred (usually the function name or class
   * name)
   */
  using Exception::Exception;
};

/**
 * @brief MatrixNotSquare
 * @brief Exception representing a non-square matrix
 *
 * this class is derived from Exception and represent an exception
 * where the matrix is not square, while the operation requires a square matrix
 */
class MatrixNotSquare : public Exception {
 public:
  /**
   * @brief get the type desciprtion of the MatrixNotSquare exception
   * @return type description of the MatrixNotSquare exception as string
   */
  std::string type_description() const override { return "Matrix is not square"; }
  /**
   * @brief constructor for the MatrixNotSquare class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class MatrixNotCvector
 * @brief Exception representing a non-column vector matrix
 *
 * this class is derived from exception and represent an exception
 * where the matrix is not a column vector, while the operation requires a column vector
 */
class MatrixNotCvector : public Exception {
 public:
  /**
   * @brief get type description of the MatrixNotCvector
   * @return the type description of the MatrixNotCvector exxception as a string
   */
  std::string type_description() const override { return "Matrix is not column vector"; }
  /**
   * @brief constructor for the MatrixNotCvector clss
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class MatrixNotRvector
 * @brief Exception representing a non-row vector matrix
 *
 * this class is derived from Exception and represent an exception
 * where the matrix is not a row vector, while the operation requires a row vector
 */
class MatrixNotRvector : public Exception {
 public:
  /**
   * @brief get the type description of the MatrixNotRvector exception
   * @return type description of the MatrixNotRvector exception as string
   */
  std::string type_description() const override { return "Matrix is not row vector"; }
  /**
   * @brief constructor for the MatrixNotRvector class
   * @param where the location where exception occured
   */
  using Exception::Exception;
};

/**
 * @class MatrixNotVector
 * @brief Exception representing a non-vector matrix
 *
 * this class is derived from Exception and represent an exception
 * wehere the matrix is not a vector, while the oepration requires a vector
 */
class MatrixNotVector : public Exception {
 public:
  /**
   * @brief get the description of the MatrixNotVector exception
   * @return the type description of MatrixNotVector exception as a string
   */
  std::string type_description() const override { return "Matrix is not vector"; }
  /**
   * @brief Construct for the MatrixNotVector class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class MatrixNotSquareNotCvector
 * @brief Exception representing a non-square matrix that is also not a column vector
 *
 * this class is derived from Exception and represent an Exception
 * where the matrix is not square and is also not a column vector
 * while the oepration requires a square matrix or column vectorx
 */
class MatrixNotSquareNotCvector : public Exception {
 public:
  /**
   * @brief get the type description of the MatrixNotSquareNotCvector exception
   * @return type description of the MatrixNotSquareNotCvector exception as a stringa
   */
  std::string type_description() const override { return "Matrix is not square not column vector"; }
  /**
   * @brief constructure for the MatrixNotCvector class
   * @param the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @brief clara::exception::MatrixNotSquareNotRvector
 * @brief matrix is not square not row vector exception
 */
class MatrixNotSquareNotRvector : public Exception {
 public:
  /**
   * @brief get the type description of the MatrixNotSquareNotRvector exception
   * @return type description of the MatrixNotSquareNotRvector exception as a string
   */
  std::string type_description() const override { return "Matrix is not square not row vector"; }
  /**
   * @brief constructor for the MatrixNotSquareNotRvector class
   * @param the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @brief  MatrixNotSquareNotVector
 * @brief Exception representing a non-square matrix that is also not a vector
 *
 * this class detived from Exception and represent an exception where the matrix is not
 * square and is also not a vector, while the operation requires a squre mtrix or a vector
 */
class MatrixNotSquareNotVector : public Exception {
 public:
  /**
   * @brief get the type description of the MatrixNotSquareNotVector exception
   * @return the type description of the MatrixNotSquareNotVector exception as string
   */
  std::string type_description() const override { return "Matrix is not square not vector"; }
  /**
   * @brief constructor fopr the MatrixNotSquareNotVector class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class MatrixMismatchSubsys
 * @brief Exception representing a matrix mismatch in a subsystem
 *
 * this class is derived from Exception and represents an exception
 * where there is a matrix mismatch in a subsystem of the application
 */
class MatrixMismatchSubsys : public Exception {
 public:
  /**
   * @brief get type description of the MatrixMismatchSubsys exception
   * @return type description of the MatrixMismatchSubsys exception as a string
   */
  std::string type_description() const override { return "Matrix mismatch subsystem"; }
  /**
   * @brief constructor for the MatrixMismatchSubsys class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class DimsInvalid
 * @brief Exception representing invalid dimensions
 *
 * this class is dervied from Exception an represents an exception
 * where the provided dimensions are invalid or out of range
 */
class DimsInvalid : public Exception {
 public:
  /**
   * @brief get the type descriptio of the DimsInvalid exception
   * @return tyhe type description the DimsInvalid exception as a string
   */
  std::string type_description() const override { return "Invalid dimension (s)"; }
  /**
   * @brief constructor for DimsInvalid class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class DimsNotEqual
 * @brief Exception representing dimensions that are not equal
 *
 * this class is derived from Exception and represents an exception
 * where the dimensions provides
 */
class DimsNotEqual : public Exception {
 public:
  /**
   * @brief get the type description of the DimsNotEqual exception
   * @return type description of the DimsNotEqual exception as string
   */
  std::string type_description() const override { return "Dimensional not equal"; }
  /**
   * @brief constructor for the DimsNotEqual class
   * @param where the location where exxception occurred
   */
  using Exception::Exception;
};

/**
 * @class DimsMismatchMatrix
 * @brief Exception representing a dimension mismatch in matrix size
 *
 * this class is derived from Exception and represent an exception
 * where the provided the dimension of matrices do not match, leading to an invalid
 * operation
 */
class DimsMismatchMatrix : public Exception {
 public:
  /**
   * @brief get the type description of the DimsMismatchMatrix exception
   * @return the type description of the DimsMismatchMatrix exxception as string
   */
  std::string type_description() const override { return "Dimension mismatch matrix size"; }
  /**
   * @brief constructor for the DimsMismatchMatrix class
   * @param where the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class DimsMismatchCvector
 * @brief Exception representing dimension mismatch in column vector size
 *
 * this class is detived from exception and represent an exception where
 * the dimension of the provided cloumn vector and the expected size
 * do not match, leading to an invalid operation
 */
class DimsMismatchCvector : public Exception {
 public:
  /**
   * @brief get the type description of the DimsMismatchCvector exception
   * @return the type description of the DimsMismatchCvector exception as string
   */
  std::string type_description() const override {
    return "Dimension(s) mismatch column vector size";
  }
  /**
   * @brief construct for the DimsMismatchCvector class
   * @param where the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class DimsMismatchRvector
 * @brief Exception representing dimension mismatch in reow vector size
 *
 * this class is derived from Exception and represent an exception
 * where the dimension of the prpvided row vector and the expected size do
 * not match, leading to an invalid operation
 */
class DimsMismatchRvector : public Exception {
 public:
  /**
   * @brief get the type description of the DimsMismatchRvector exception
   * @return the type description of the DimsMismatchRvector exception as string
   */
  std::string type_description() const override { return "Dimension(s) mismatch row vector size"; }
  /**
   * @brief constructor for DimsMismatchRvector class
   * @param where the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class DimsMismatchRvector
 * @brief Exception representing dimension in vector size
 *
 * this class dervied from Exception and represents an exception
 * where the dimension of the provided vector and the expected size
 * do not match, leading to an invalid operation
 */
class DimsMismatchVector : public Exception {
 public:
  /**
   * @brief get type description of the DimsMismatchVector exception
   * @return the type description of the DimsMismatchVector exception as string
   */
  std::string type_description() const override { return "Dimension(s) mismatch vector size"; }
  /**
   * @brief constructor for the DimsMismatchVector class
   * @param where the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class SubsysMismatchdims
 * @brief Exception representing subsystem mismatch in dimensions
 *
 * this class is derived from exception and represents an exception
 * where the dimension of the subsystem do not match as expired
 * leading to an invalid operations
 */
class SubsysMismatchdims : public Exception {
 public:
  /**
   * @brief get the type description of the SubsysMismatchdims exception
   * @return the type description of the SubsysMismatchdims exception as string
   */
  std::string type_description() const override { return "Subsytem mismatch dimensions"; }
  /**
   * @brief constructor for the SubsysMismatchdims class
   * @param the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class PermInvalid
 * @brief Exception representing an invalid permutation
 *
 * this class is derived from Exception and represent an exception
 * where an invalid permutation is provided, leading to an invalid operation
 */
class PermInvalid : public Exception {
 public:
  /**
   * @brief th type description of the PermInvalid exception
   * @return the type description of the PermInvalid exception as string
   */
  std::string type_description() const override { return "Invalid permutation"; }
  /**
   * @brief constructor for the PermInvalid class
   * @param where the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class clara::exception::PermMismatchDims
 * @brief permutateion mismatch dimension exception
 */
class PermMismatchDims : public Exception {
 public:
  std::string type_description() const override { return "Permutation mismatch dimensions"; }
  using Exception::Exception;
};

/**
 * @class PermMismatchDims
 * @brief Exception representing permutation mismatch in dimension
 *
 * this class is derived from Exception and represent an exception
 * where the dimension of the permutation do not match expected
 * leading to an invalid operation
 */
class NotQubitMatrix : public Exception {
 public:
  /**
   * @brief get the type description of the PermMismatchDims exception
   * @return the type description of the PermMismatchDims exception as string
   */
  std::string type_description() const override { return "Matrix is not 2 x 2"; }
  /**
   * @brief get the type description of the PermMismatchDims exception
   * @return the type description of the PermMismatchDims exceptin as string
   */
  using Exception::Exception;
};

/**
 * @class NotQubitCvector
 * @brief Exception representing  column vector that is not 2x1
 *
 * this class is derived from exception and represents as exception
 * where a column vector is provided, but its dimensions are not 2x1
 * which is required for a qubit representation
 */
class NotQubitCvector : public Exception {
 public:
  /**
   * @brief get the type description of NotQubitCvector exception
   * @return type description of the NotQubitCvector exception as string
   */
  std::string type_description() const override { return "Column vector is not 2 x 1"; }
  /**
   * @brief constructor for the not qubit NotQubitCvector class
   * @param where the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class NotQubitVector
 * @brief Exception representing a row vector that is not 1 x 2
 *
 * this class is derived from Exception and represent an exception
 * where a row vector is provided, but its dimensions are not 1x2.
 * which is required for a qubit representation
 */
class NotQubitRvector : public Exception {
 public:
  /**
   * @brief constructor for the NotQubitRvector class
   */
  std::string type_description() const override { return "Row vector is not 1 x 2"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::NotQubitVector
 * @brief vector is not 2 x 1 nor 1 x 2 exception
 *
 * this class is derived from Exception and represents an exception
 * where a vector is provided, but its dimensions are not 2x1 nor 1x2,
 * which are required for a qubit representation.
 */
class NotQubitVector : public Exception {
 public:
  std::string type_description() const override { return "Vector is not 2 x 1 nor 1 x 2"; }
  using Exception::Exception;
};

/**
 * @class NotQubitSubsys
 * @brief Exception representing subsystem that are not qubits
 *
 * this class is derived from exception and represent an exception
 * where the subsystem provided are not qubit as expected
 */
class NotQubitSubsys : public Exception {
 public:
  /**
   * @brief get the type description of NotQubitSubsys exception
   * @preturn the type description of the NotQubitSubsys exception as string
   */
  std::string type_description() const override { return "Subsytem are not qubits"; }
  /**
   * @brief constructor for the NotQubitSubsys class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class NotBipartite
 * @brief Exception representing a situation thtis not bipartite
 *
 * this class derived from exception and represents an exception
 * where the given situation is not bipartite as expected
 */
class NotBipartite : public Exception {
 public:
  /**
   * @brief get the type description of the NotBipartite exception
   * @return the type description of the NotQubitRvector exception as string
   */
  std::string type_description() const override { return "Not bi-bipartite"; }
  /**
   * @brief constructor for the NotBipartite class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class NoCodeword
 * @brief Exception representing the absence of a specific Codeword
 *
 * this class is derived from exception and represent an exception
 * where a specific codeword does not exist as expected
 */
class NoCodeword : public Exception {
 public:
  /**
   * @brief get the type description of the NoCodeword exception
   * @return the type description of the NoCodeword exception as a string
   */
  std::string type_description() const override { return "Codeword does not exist"; };
  /**
   * @brief constructor for the NoCodeword class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class OutOfRange
 * @brief Exception representing a parameter out of range
 *
 * this class is derived from Exception and represents an exception
 * where a parameter provided is out of the valid range
 */
class OutOfRange : public Exception {
 public:
  /**
   * @brief get the type description of the OutOfRange exception
   * @param the type description of the OutOfRange exception as string
   */
  std::string type_description() const override { return "Parameter out of range"; }
  /**
   * @brief constructor of the OutOfRange class
   * @param where the location where the exception occured
   */
  using Exception::Exception;
};

/**
 * @class TypeMismatch
 * @brief Exception representing a type mismatch
 * this class is derived from Exception and represent an exception
 * where a type mismatch is encountered
 */
class TypeMismatch : public Exception {
 public:
  /**
   * @brief get type description of the TypeMismatch exception
   * @return the type description of the TypeMismatch exxception as string
   */
  std::string type_description() const override { return "Type mismatch"; }
  /**
   * @brief constructor for the TypeMismatch class
   * @param where the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class SizeMismatch
 * @brief Exception representing a size mismatch
 *
 * this class is derived from Exception and represent an exception
 * where a size mismatch is encountered
 */
class SizeMismatch : public Exception {
 public:
  /**
   * @brief get the type description of the SizeMismatch exception
   * @return the type description of the SizeMismatch exception as a string
   */
  std::string type_description() const override { return "Size mismatch"; }
  /**
   * @brief construct for the SizeMismatch class
   * @param where the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class UndefinedType
 * @brief Exception representing an operation not defined for this type
 */
class UndefinedType : public Exception {
 public:
  /**
   * @brief get the type description of the UndefinedType exception
   * @return the type description of the UndefinedType exception as string
   */
  std::string type_description() const override { return "Not defined for this type"; }
  /**
   * @brief constructor for UndefinedType class
   * @param where the location where the exception occurred
   */
  using Exception::Exception;
};

/**
 * @class CustomException
 * @brief Custom exception with a user-defined message
 *
 * this class is derived from the Exception base class and represents a custom exception with a
 * user-defined message the user can provide custom messagetht will be included in the error
 * description
 */
class CustomException : public Exception {
  std::string what_{};
  /**
   * @brief get the type description of the CustomException
   * @return the type description of the CustomException as string
   */
  std::string type_description() const override { return "CUSTOM EXCEPTION " + what_; }

 public:
  /**
   * @brief constructor for the CustomExceptio class
   * @param where the location where the exception occurred
   * @param what the custom message describing the exception
   */
  CustomException(const std::string& where, const std::string& what)
      : Exception{where}, what_{what} {}
};

}  // namespace exception
}  // namespace clara

#endif  // !CLASSFUNCTION_EXCEPTION_H_
