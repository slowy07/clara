#ifndef CLASSFUNCTION_EXCEPTION_H_
#define CLASSFUNCTION_EXCEPTION_H_

#include <exception>
#include <string>
namespace clara {

class Exception : public std::exception {
 public:
  enum class Type {
    UNKNOWN_EXCEPTION = 1,
    ZERO_SIZE,
    MATRIX_NOT_SQUARE,
    MATRIX_NOT_CVECTOR,
    MATRIX_NOT_RVECTOR,
    MATRIX_NOT_VECTOR,
    MATRIX_NOT_SQUARE_OR_CVECTOR,
    MATRIX_NOT_SQUARE_OR_RVECTOR,
    MATRIX_NOT_SQUARE_OR_VECTOR,
    MATRIX_MISMATCH_SUBSYS,
    DIMS_INVALID,
    DIMS_NOT_EQUAL,
    DIMS_MISMATCH_MATRIX,
    DIMS_MISMATCH_CVECTOR,
    DIMS_MISMATCH_RVECTOR,
    DIMS_MISMATCH_VECTOR,
    SUBSYS_MISMATCH_DIMS,
    PERM_INVALID,
    PERM_MISMATCH_DIMS,
    NOT_QUBIT_MATRIX,
    NOT_QUBIT_CVECTOR,
    NOT_QUBIT_RVECTOR,
    NOT_QUBIT_VECTOR,
    NOT_QUBIT_SUBSYS,
    NOT_BIPARTITE,
    NO_CODEWORD,
    OUT_OF_RANGE,
    TYPE_MISMATCH,
    SIZE_MISMATCH,
    UNDEFINED_TYPE,
    CUSTOM_EXCEPTION
  };

  /**
   * construct exception
   */
  Exception(const std::string& where, [[maybe_unused]] const Type& type)
      : where_{where}, msg_{}, type_{}, custom_{} {
    construct_exception_msg_();
  }

  Exception(const std::string& where, const std::string& custom)
      : where_{where}, msg_{}, type_{Type::CUSTOM_EXCEPTION}, custom_{custom} {
    construct_exception_msg_();
  }

  virtual const char* what() const noexcept override { return msg_.c_str(); }

 private:
  std::string where_, msg_;
  Type type_;
  std::string custom_;

  void construct_exception_msg_() {
    msg_ += "IN ";
    msg_ += where_;
    msg_ += ": ";
    switch (type_) {
      case Type::UNKNOWN_EXCEPTION:
        msg_ += "UNKNOWN EXCEPTION!";
        break;
      case Type::ZERO_SIZE:
        msg_ += "Object has zero size!";
        break;
      case Type::MATRIX_NOT_SQUARE:
        msg_ += "Matrix is not square!";
        break;
      case Type::MATRIX_NOT_CVECTOR:
        msg_ += "Matrix is not column vector!";
        break;
      case Type::MATRIX_NOT_RVECTOR:
        msg_ += "Matrix is not row vector!";
        break;
      case Type::MATRIX_NOT_VECTOR:
        msg_ += "Matrix is not vector!";
        break;
      case Type::MATRIX_NOT_SQUARE_OR_CVECTOR:
        msg_ += "Matrix is not square nor column vector!";
        break;
      case Type::MATRIX_NOT_SQUARE_OR_RVECTOR:
        msg_ += "Matrix is not square nor row vector!";
        break;
      case Type::MATRIX_NOT_SQUARE_OR_VECTOR:
        msg_ += "Matrix is not square nor vector!";
        break;
      case Type::MATRIX_MISMATCH_SUBSYS:
        msg_ += "Matrix mismatch subsystems!";
        break;
      case Type::DIMS_INVALID:
        msg_ += "Invalid dimension(s)!";
        break;
      case Type::DIMS_NOT_EQUAL:
        msg_ += "Dimensions not equal!";
        break;
      case Type::DIMS_MISMATCH_MATRIX:
        msg_ += "Dimension(s) mismatch matrix size!";
        break;
      case Type::DIMS_MISMATCH_CVECTOR:
        msg_ += "Dimension(s) mismatch column vector!";
        break;
      case Type::DIMS_MISMATCH_RVECTOR:
        msg_ += "Dimension(s) mismatch row vector!";
        break;
      case Type::DIMS_MISMATCH_VECTOR:
        msg_ += "Dimension(s) mismatch vector!";
        break;
      case Type::SUBSYS_MISMATCH_DIMS:
        msg_ += "Subsystems mismatch dimensions!";
        break;
      case Type::PERM_INVALID:
        msg_ += "Invalid permutation!";
        break;
      case Type::PERM_MISMATCH_DIMS:
        msg_ += "Permutation mismatch dimensions!";
        break;
      case Type::NOT_QUBIT_MATRIX:
        msg_ += "Matrix is not 2 x 2!";
        break;
      case Type::NOT_QUBIT_CVECTOR:
        msg_ += "Column vector is not 2 x 1!";
        break;
      case Type::NOT_QUBIT_RVECTOR:
        msg_ += "Row vector is not 1 x 2!";
        break;
      case Type::NOT_QUBIT_VECTOR:
        msg_ += "Vector is not 2 x 1 nor 1 x 2!";
        break;
      case Type::NOT_QUBIT_SUBSYS:
        msg_ += "Subsystems are not qubits!";
        break;
      case Type::NOT_BIPARTITE:
        msg_ += "Not bi-partite!";
        break;
      case Type::NO_CODEWORD:
        msg_ += "Codeword does not exist!";
        break;
      case Type::OUT_OF_RANGE:
        msg_ += "Parameter out of range!";
        break;
      case Type::TYPE_MISMATCH:
        msg_ += "Type mismatch!";
        break;
      case Type::SIZE_MISMATCH:
        msg_ += "Size mismatch!";
        break;
      case Type::UNDEFINED_TYPE:
        msg_ += "Not defined for this type!";
        break;
      case Type::CUSTOM_EXCEPTION:
        msg_ += "CUSTOM EXCEPTION ";
        break;
    }
  }
};

}  // namespace clara

#endif  // !CLASSFUNCTION_EXCEPTION_H_
