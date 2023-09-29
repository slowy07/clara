#include <cmath>
#include <initializer_list>
#include <ios>
#include <iterator>
#include <list>
#include <vector>

#include "../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

class ComplexNumber {
 public:
  ComplexNumber(double real, double imag) : real_(real), imag_(imag) {}
  double real() const { return real_; }
  double image() const { return imag_; }

 private:
  double real_;
  double imag_;
};

TEST(clara_is_complex, AllTests) {
  // test if provided types are complex numbers
  EXPECT_TRUE(clara::is_complex<std::complex<double>>::value);
  EXPECT_TRUE(clara::is_complex<std::complex<int>>::value);

  // custom ComplexNumber class is not complex number
  EXPECT_FALSE(clara::is_complex<ComplexNumber>::value);
}

TEST(clara_is_iterable, AllTests) {
  // custom iterable class for testing
  class CustomIterable {
   public:
    using value_type = int;
  
    const int* begin() {return &data[0];}
    const int* end() {return &data[size]; }

    const int* begin() const {return &data[0];}
    const int* end() const {return &data[size];}
    
    CustomIterable(std::initializer_list<int> initList) : data(initList), size(initList.size()) {}

  private:
    std::vector<int> data;
    std::size_t size;
  };

  EXPECT_FALSE(clara::is_iterable<int>::value);
  // std::vector<int> is iterable
  EXPECT_TRUE(clara::is_iterable<std::vector<int>>::value);
  // std::list<double> is iterable
  EXPECT_TRUE(clara::is_iterable<std::list<double>>::value);
  // CustomIterable is iterable
  EXPECT_TRUE(clara::is_iterable<CustomIterable>::value);
}

TEST(clara_is_matrix_expression, AllTests) {
  dmat A, B;
  int x{}, y{}, z{};

  class NonMatrixCustomClass {
  public:
    int getValue() const {return value_;}
    void setValue(int value) {value_ = value;}
  private:
    int value_;
  };

  // scalar multiplication is a matrix expression
  EXPECT_TRUE(clara::is_matrix_expression<decltype(3 * A)>::value);
  // matrix addition is matrix expression
  EXPECT_TRUE(clara::is_matrix_expression<decltype(A + B)>::value);
  // x + y * z is not a matrix expression
  EXPECT_FALSE(clara::is_matrix_expression<decltype(x + y * z)>::value);
  // NonMatrixCustomClass is not a matrix expression
  EXPECT_FALSE(clara::is_matrix_expression<NonMatrixCustomClass>::value);
}
