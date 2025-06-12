#include <gtest/gtest.h>
#include <sstream>
#include <beam/LinAlg/Matrix.hpp>

using beam::Matrix;
using beam::MatrixOrder;

TEST(MatrixTest, ConstructionAndSize) {
  Matrix<double> m(3, 4, MatrixOrder::RowMajorOrder);
  EXPECT_EQ(m.rows(), 3);
  EXPECT_EQ(m.cols(), 4);
  EXPECT_EQ(m.size(), 12);
}

TEST(MatrixTest, InitializerListConstructor) {
  Matrix<double> m{{1.0, 2.0}, {3.0, 4.0}};
  EXPECT_EQ(m.rows(), 2);
  EXPECT_EQ(m.cols(), 2);
  EXPECT_EQ(m.size(), 4);
  EXPECT_DOUBLE_EQ(m(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(m(1, 0), 3.0);
  EXPECT_DOUBLE_EQ(m(1, 1), 4.0);
}

TEST(MatrixTest, CopyConstructor) {
  Matrix<double> m1{{1.0, 2.0}, {3.0, 4.0}};
  Matrix<double> m2(m1);
  EXPECT_EQ(m2.rows(), m1.rows());
  EXPECT_EQ(m2.cols(), m1.cols());
  EXPECT_EQ(m2.size(), m1.size());
  EXPECT_DOUBLE_EQ(m2(1, 1), 4.0);
}

TEST(MatrixTest, CopyAssignment) {
  Matrix<double> m1{{1.0, 2.0}, {3.0, 4.0}};
  Matrix<double> m2(2, 2, MatrixOrder::RowMajorOrder);
  m2 = m1;
  EXPECT_EQ(m2.rows(), 2);
  EXPECT_EQ(m2.cols(), 2);
  EXPECT_EQ(m2.size(), 4);
  EXPECT_DOUBLE_EQ(m2(0, 1), 2.0);
}

TEST(MatrixTest, MoveConstructor) {
  Matrix<double> m1{{1.0, 2.0}, {3.0, 4.0}};
  Matrix<double> m2(std::move(m1));
  EXPECT_EQ(m2.rows(), 2);
  EXPECT_EQ(m2.cols(), 2);
  EXPECT_EQ(m2.size(), 4);
  EXPECT_DOUBLE_EQ(m2(0, 0), 1.0);
}

TEST(MatrixTest, MoveAssignment) {
  Matrix<double> m1{{1.0, 2.0}, {3.0, 4.0}};
  Matrix<double> m2(2, 2, MatrixOrder::RowMajorOrder);
  m2 = std::move(m1);
  EXPECT_EQ(m2.rows(), 2);
  EXPECT_EQ(m2.cols(), 2);
  EXPECT_EQ(m2.size(), 4);
  EXPECT_DOUBLE_EQ(m2(1, 0), 3.0);
}

TEST(MatrixTest, OutOfBoundsAccess) {
  Matrix<double> m(2, 2, MatrixOrder::RowMajorOrder);
  EXPECT_THROW(m(2, 0), std::out_of_range);
  EXPECT_THROW(m(0, 2), std::out_of_range);
}

TEST(MatrixTest, PrettyPrint) {
  Matrix<double> m{{1.0, 2.0}, {3.0, 4.0}};
  std::ostringstream oss;
  oss << m;
  std::string output = oss.str();
  // Check that expected values are in the output string.
  EXPECT_NE(output.find("1"), std::string::npos);
  EXPECT_NE(output.find("2"), std::string::npos);
  EXPECT_NE(output.find("3"), std::string::npos);
  EXPECT_NE(output.find("4"), std::string::npos);
}