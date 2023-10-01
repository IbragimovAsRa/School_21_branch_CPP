#include "test.h"

TEST(MatrixMethodTest, MethodSumMatrix1) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double mres[] = {9, 8, 8, 3};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix MRes(2, 2, mres);

  M1.SumMatrix(M2);
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodSumMatrix2) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double mres[] = {9, 8, 8, 3};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix MRes(2, 2, mres);

  M1 = M1 + M2;
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodSumMatrix3) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double mres[] = {9, 8, 8, 3};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix MRes(2, 2, mres);

  M1 += M2;
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodSubMatrix1) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double mres[] = {1, -4, -2, -1};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix MRes(2, 2, mres);

  M1.SubMatrix(M2);
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodSubMatrix2) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double mres[] = {1, -4, -2, -1};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix MRes(2, 2, mres);

  M1 = M1 - M2;
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodSubMatrix3) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double mres[] = {1, -4, -2, -1};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix MRes(2, 2, mres);

  M1 -= M2;
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodMulNumber1) {
  double m[] = {5, 2, 3, 1};
  double mres[] = {15, 6, 9, 3};

  S21Matrix M1(2, 2, m);
  S21Matrix MRes(2, 2, mres);

  M1.MulNumber(3);
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodMulNumber2) {
  double m[] = {5, 2, 3, 1};
  double mres[] = {15, 6, 9, 3};

  S21Matrix M1(2, 2, m);
  S21Matrix MRes(2, 2, mres);

  M1 = M1 * 3;
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodMulNumber3) {
  double m[] = {5, 2, 3, 1};
  double mres[] = {15, 6, 9, 3};

  S21Matrix M1(2, 2, m);
  S21Matrix MRes(2, 2, mres);

  M1 *= 3;
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodMulMatrix1) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double mres[] = {30, 34, 17, 20};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix MRes(2, 2, mres);

  M1.MulMatrix(M2);
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodMulMatrix2) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double mres[] = {30, 34, 17, 20};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix MRes(2, 2, mres);

  M1 = M1 * M2;
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodMulMatrix3) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double mres[] = {30, 34, 17, 20};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix MRes(2, 2, mres);

  M1 *= M2;
  EXPECT_EQ(M1.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodEqMatrix1) {
  double m1[] = {5, 2, 3, 1};
  double m2[] = {4, 6, 5, 2};
  double m3[] = {5, 2, 3, 1};

  S21Matrix M1(2, 2, m1);
  S21Matrix M2(2, 2, m2);
  S21Matrix M3(2, 2, m3);

  EXPECT_EQ(M1.EqMatrix(M3), 1);
  EXPECT_EQ(M1 == M3, 1);

  EXPECT_EQ(M1.EqMatrix(M2), 0);
  EXPECT_EQ(M1 == M2, 0);
}

TEST(MatrixMethodTest, MethodTranspose1) {
  double m[] = {5, 2, 3, 1};
  double mres[] = {5, 3, 2, 1};

  S21Matrix M1(2, 2, m);
  S21Matrix MRes(2, 2, mres);

  S21Matrix M2 = M1.Transpose();
  EXPECT_EQ(M2.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodCalcComplements1) {
  double m[] = {5, 2, 3, 1};
  double mres[] = {1, -3, -2, 5};

  S21Matrix M1(2, 2, m);
  S21Matrix MRes(2, 2, mres);

  S21Matrix M2 = M1.CalcComplements();
  EXPECT_EQ(M2.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodInverseMatrix1) {
  double m[] = {5, 2, 3, 1};
  double mres[] = {-1, 2, 3, -5};

  S21Matrix M1(2, 2, m);
  S21Matrix MRes(2, 2, mres);

  S21Matrix M2 = M1.InverseMatrix();
  EXPECT_EQ(M2.EqMatrix(MRes), 1);
}

TEST(MatrixMethodTest, MethodIndex1) {
  double m1[] = {5, 2, 3, 1};
  double m2[2][2] = {{5, 2}, {3, 1}};

  S21Matrix M1(2, 2, m1);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      EXPECT_EQ(abs((M1(i, j) - m2[i][j])) < EPS, 1);
    }
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
