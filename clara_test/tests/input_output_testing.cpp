#include <fstream>

#include "../../include/clara.h"
#include "../../include/input_output.h"
#include "gtest/gtest.h"

using namespace clara;

class InputOutputTest : public ::testing::Test {
 protected:
  void SetUp() override {}
  void TearDown() override {}
};

TEST_F(InputOutputTest, SaveAndLoadEigenMatrix) {
  Eigen::MatrixXd mat(3, 3);
  mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;

  const std::string filename = "test_matrix.bin";
  clara::save(mat, filename);

  auto loadedMatrix = clara::load<Eigen::MatrixXd>(filename);

  ASSERT_EQ(mat.rows(), loadedMatrix.rows());
  ASSERT_EQ(mat.cols(), loadedMatrix.cols());

  for (int i = 0; i < mat.rows(); ++i) {
    for (int j = 0; j < mat.cols(); ++j) {
      ASSERT_DOUBLE_EQ(mat(i, j), loadedMatrix(i, j));
    }
  }
}
