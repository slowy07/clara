#include <cmath>
#include <initializer_list>
#include <iterator>
#include <list>
#include <vector>

#include "../../include/clara.h"
#include "gtest/gtest.h"

using namespace clara;

TEST(clara_is_complex, AllTests) {
  class ComplexNumber {
   public:
    ComplexNumber(double real, double imag) : real_(real), imag_(imag) {}

    double real() const { return real_; }
    double imag() const { return imag_; }

   private:
    double real_;
    double imag_;
  };
  EXPECT_TRUE(clara::is_complex<std::complex<double>>::value);
  EXPECT_TRUE(clara::is_complex<std::complex<int>>::value);
}

TEST(clara_is_iterable, AllTests) {
  class CustomIterable {
   public:
    using value_type = int;
    int* begin() { return &data[0]; }

    int* end() { return &data[size]; }

    CustomIterable(std::initializer_list<int> initList) : data(initList), size(initList.size()) {}

   private:
    std::vector<int> data;
    std::size_t size;
  };
  EXPECT_TRUE(clara::is_iterable<CustomIterable>::value);

  EXPECT_FALSE(clara::is_iterable<int>::value);
  EXPECT_TRUE(clara::is_iterable<std::vector<int>>::value);
  EXPECT_TRUE(clara::is_iterable<std::list<double>>::value);
  EXPECT_TRUE(clara::is_iterable<std::set<std::string>>::value);
}

TEST(clara_is_matrix_expression, AllTests) {
  dmat A, B;
  int x{}, y{}, z{};

  class NonMatrixCustomClass {
   public:
    int getValue() const { return value_; }
    void setValue(int value) { value_ = value; }

   private:
    int value_;
  };

  EXPECT_TRUE(clara::is_matrix_expression<decltype(3 * A)>::value);
  EXPECT_TRUE(clara::is_matrix_expression<decltype(A + B)>::value);
  EXPECT_FALSE(clara::is_matrix_expression<decltype(x + y * z)>::value);
  EXPECT_FALSE(clara::is_matrix_expression<NonMatrixCustomClass>::value);
}
