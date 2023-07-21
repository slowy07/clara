#ifndef CLASSFUNCTION_EXCEPTION_H_
#define CLASSFUNCTION_EXCEPTION_H_

#include <exception>
#include <string>
namespace clara {

namespace exception {

class Exception : public std::exception {
 private:
  std::string where_;
  mutable std::string msg_;

 public:
  Exception(const std::string& where) : where_{where}, msg_{} {}

  virtual const char* what() const noexcept override {
    msg_.clear();
    msg_ += "IN ";
    msg_ += where_;
    msg_ += ": ";
    msg_ += this->type_description();
    msg_ += "!";
    return msg_.c_str();
  }
  virtual std::string type_description() const = 0;
};

inline std::string Exception::type_description() const { return "clara::exception::Exception"; }

/**
 * @class clara::exception::Unknwon
 * thrown when no other exception is suitable
 */
class Unknown : public Exception {
 public:
  std::string type_description() const override { return "UNKNOWN EXCEPTION"; }
  using Exception::Exception;
};

/**
 * @brief exception type description
 * @return object has zero size exception
 */

class ZeroSize : public Exception {
 public:
  std::string type_description() const override { return "Object has zero size"; }
  using Exception::Exception;
};

/**
 * @brief clara::exception::MatrixNotSquare
 * @brief matrix is not square exception
 */
class MatrixNotSquare : public Exception {
 public:
  std::string type_description() const override { return "Matrix is not square"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::MatrixNotCvector
 * @brief matrix is not column vector exception
 */

class MatrixNotCvector : public Exception {
 public:
  std::string type_description() const override { return "Matrix is not column vector"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::MatrixRowVector
 * @brief matrix is not a row vector exception
 */
class MatrixNotRvector : public Exception {
 public:
  std::string type_description() const override { return "Matrix is not row vector"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::MatrixxNotVector
 * @brief matrix is not vector Exception
 */
class MatrixNotVector : public Exception {
 public:
  std::string type_description() const override { return "Matrix is not vector"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::MatrixSquareNotCvector
 * @brief matrix is not column vector exception
 */
class MatrixNotSquareNotCvector : public Exception {
 public:
  std::string type_description() const override { return "Matrix is not square not column vector"; }
  using Exception::Exception;
};

/**
 * @brief clara::exception::MatrixNotSquareNotRvector
 * @biref matrix is not square not row vector exception
 */
class MatrixNotSquareNotRvector : public Exception {
 public:
  std::string type_description() const override { return "Matrix is not square not row vector"; }
  using Exception::Exception;
};

/**
 * @brief clara::Exception::MatrixNotSquareNotVector
 * @brief matrix is not square not vector exception
 */
class MatriNotSquareNotVector : public Exception {
 public:
  std::string type_description() const override { return "Matrix is not square not vector"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::MatrixMismatchSubsys
 * @brief matrix mismatch subsytem exception
 */
class MatrixMismatchSubsys : public Exception {
 public:
  std::string type_description() const override { return "Matrix mismatch subsystem"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::DimsInvalid
 * @brief invalid dimension exception
 */
class DimsInvalid : public Exception {
 public:
  std::string type_description() const override { return "Invalid dimension (s)"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::DimsNotEqual
 * @brief dimension not equal exception
 */
class DimsNotEqual : public Exception {
 public:
  std::string type_description() const override { return "Dimensional not equal"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::DimsMismatchMatrix
 * @brief dimension mismatch matrix size exception
 */
class DimsMismatchMatrix : public Exception {
 public:
  std::string type_description() const override { return "Dimension mismatch matrix size"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::DimsMismatchCvector
 * @brief Dimension(s) mismatch column vector size exception
 * product of the elements of std::vector<idx> of dimension is not equal to
 * the number of elements of the Eigen::Matrix
 */
class DimsMismatchCvector : public Exception {
 public:
  std::string type_description() const override {
    return "Dimension(s) mismatch column vector size";
  }
  using Exception::Exception;
};

/**
 * @class clara::exception::DimsMismatchRvector
 * @brief Dimension(s) mismatch row vector size exception
 * product of the elements of std::vector<idx> of dimension is not equal
 * to the number of the elements of the Eigen::Matrix
 */
class DimsMismatchRvector : public Exception {
 public:
  std::string type_description() const override { return "Dimension(s) mismatch row vector size"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::DimsMismatchVector
 * @brief Dimension mismatch vector size exception
 * prodcut of the element of std::vector<idx> of dimension is not equal
 */
class DimsMismatchVector : public Exception {
 public:
  std::string type_description() const override { return "Dimension(s) mismatch vector size"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::SubsysMismatchDims
 * @brief Subsystem mismatch dimension exception
 */
class SubsysMismatchdims : public Exception {
 public:
  std::string type_description() const override { return "Subsytem mismatch dimensions"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::PermInvalid
 * @brief invalid permutation exception
 */
class PermInvalid : public Exception {
 public:
  std::string type_description() const override { return "Invalid permutation"; }
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
 * @class clara::exception::NotQubitMatrix
 * @brief matrix is not 2 x 2 exception
 */
class NotQubitMatrix : public Exception {
 public:
  std::string type_description() const override { return "Matrix is not 2 x 2"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::NOtQubitCvector
 * @brief column vector is not 2 x 1 exception
 */
class NOtQubitCvector : public Exception {
 public:
  std::string type_description() const override { return "Column vector is not 2 x 1"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::NotQubitRvector
 * @brief row vector is not 1 x 2 exception
 */
class NotQubitRvector : public Exception {
 public:
  std::string type_description() const override { return "Row vector is not 1 x 2"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::NotQubitVector
 * @brief vector is not 2 x 1 nor 1 x 2 exception
 */
class NotQubitVector : public Exception {
 public:
  std::string type_description() const override { return "Vector is not 2 x 1 nor 1 x 2"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::NotQubitSubsys
 * @brief Subystem are not qubits exception
 */
class NotQubitSubsys : public Exception {
 public:
  std::string type_description() const override { return "Subsytem are not qubits"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::NotBipartite
 * @brief not bi-bipartite exception
 */
class NotBipartite : public Exception {
 public:
  std::string type_description() const override { return "Not bi-bipartite"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::NoCodeword
 * @brief codeword does not exist exception
 */
class NoCodeword : public Exception {
 public:
  std::string type_description() const override { return "Codeword does not exist"; };
  using Exception::Exception;
};

/**
 * @class clara::exception::OutOfRange
 * @brief paramtere out of range exception
 */
class OutOfRange : public Exception {
 public:
  std::string type_description() const override { return "Parameter out of range"; }
  using Exception::Exception;
};

class TypeMismatch : public Exception {
 public:
  std::string type_description() const override { return "Type mismatch"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::SizMismatch
 * @brief size mismatch exception
 */
class SizeMismatch : public Exception {
 public:
  std::string type_description() const override { return "Size mismatch"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::UndefinedType
 * @brief not defined for this type exception
 */
class UndefinedType : public Exception {
 public:
  std::string type_description() const override { return "Not defined for this type"; }
  using Exception::Exception;
};

/**
 * @class clara::exception::CustomException
 * @brief custom exception
 */
class CustomException : public Exception {
  std::string what_{};
  std::string type_description() const override { return "CUSTOM EXCEPTION " + what_; }

 public:
  CustomException(const std::string& where, const std::string& what)
      : Exception{where}, what_{what} {}
};

}  // namespace exception
}  // namespace clara

#endif  // !CLASSFUNCTION_EXCEPTION_H_
