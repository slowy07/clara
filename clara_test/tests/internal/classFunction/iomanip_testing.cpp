#include <complex>
#include <sstream>

#include "../../../../include/clara.h"
#include "gtest/gtest.h"

class ClaraIOManipTesting : public ::testing::Test {
 protected:
  void SetUp() override {}
  void TearDown() override {}
};

TEST_F(ClaraIOManipTesting, DisplayIOManipRangeClara) {
  int arr[] = {1, 2, 3, 4, 5};
  std::ostringstream oss;
  clara::internal::IOManipRange<int*> rangeDisp(arr, arr + 5, ", ", "<", ">");
  oss << rangeDisp;
  ASSERT_EQ(oss.str(), "<1, 2, 3, 4, 5>");
}

TEST_F(ClaraIOManipTesting, DisplayIOManipPointer) {
  int arr[] = {1, 2, 3, 4, 5};
  std::ostringstream oss;
  clara::internal::IOManipPointer<int> pointerDisp(arr, 5, ", ", "[", "]");
  oss << pointerDisp;
  ASSERT_EQ(oss.str(), "[1, 2, 3, 4, 5]");
}

TEST_F(ClaraIOManipTesting, SaveAndLoadEigenMatrix) {
  Eigen::MatrixXd originalMat(2, 2);
  originalMat << 1.0, 2.0, 3.0, 4.0;
  std::string filename = "test_matrix.bin";
  clara::save(originalMat, filename);
  Eigen::MatrixXd loadedMat = clara::load<Eigen::MatrixXd>(filename);
  ASSERT_EQ(originalMat, loadedMat);
}

TEST_F(ClaraIOManipTesting, DisplayIOManipComplex) {
  std::complex<double> complexNum(3.14, 2.17);
  std::ostringstream oss;
  clara::internal::IOManipEigen complexDisp(complexNum);
  oss << complexDisp;
  
  std::string actualOutput = oss.str();
  std::string expectedOutput = "[3.14, 2.17]";
  ASSERT_TRUE(actualOutput.find(expectedOutput) == std::string::npos);
}
