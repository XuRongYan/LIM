//
// Created by 徐溶延 on 2020/6/14.
//
#include <vector>

#include <gtest/gtest.h>

#include "../utils/EigenUtils.h"

class EigenUtilsTest : public ::testing::Test {
public:
protected:
    void SetUp() override {
        Test::SetUp();
        Vt_.resize(5, 3);
        Vt_ << 1, 2, 3,
               4, 5, 6,
               7, 8, 9,
               10, 11, 12,
               13, 14, 15;
        V_.resize(15);
        V_ << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
        for (size_t i = 0; i < 4; i++) {
            triplets.emplace_back(i, i, 1);
        }
    }

public:
    Eigen::MatrixX3f Vt_;
    Eigen::VectorXf V_;
    std::vector<Eigen::Triplet<float>> triplets;
};

TEST_F(EigenUtilsTest, vt2vTest) {
    Eigen::VectorXf v = xry_mesh::vt2v<float>(Vt_);
    //判断大小
    ASSERT_EQ(v.rows(), V_.rows());
    ASSERT_EQ(v.cols(), V_.cols());
    //判断数值
    for (size_t i = 0; i < V_.rows(); i++) {
        ASSERT_FLOAT_EQ(v[i], V_[i]);
    }
}

TEST_F(EigenUtilsTest, v2vtTest) {
    Eigen::MatrixX3f vt = xry_mesh::v2vt<float>(V_, 3);
    //判断大小
    ASSERT_EQ(vt.rows(), Vt_.rows());
    ASSERT_EQ(vt.cols(), Vt_.cols());
    //判断数值
    for (size_t r = 0; r < Vt_.rows(); r++) {
        for (size_t c = 0; c < Vt_.cols(); c++) {
            ASSERT_FLOAT_EQ(vt(r, c), Vt_(r, c));
        }
    }
}

TEST_F(EigenUtilsTest, sparse2tripletTest) {
    Eigen::SparseMatrix<float> L_regular(4, 4), L_nonregular(4, 4);
    L_regular.setFromTriplets(triplets.begin(), triplets.end());
    L_nonregular.setFromTriplets(triplets.begin(), triplets.end() - 1);
    ASSERT_TRUE(xry_mesh::isSparseMatrixInvertible(L_regular));
    ASSERT_FALSE(xry_mesh::isSparseMatrixInvertible(L_nonregular));
}

TEST_F(EigenUtilsTest, v2dTo3dTest) {
    Eigen::Matrix2Xf A(2, 2);
    A << 1, 2,
         3, 4;
    Eigen::Matrix3Xf A3d = xry_mesh::v2dTo3d<float>(A);
    for (size_t r = 0; r < A.rows(); r++) {
        for (size_t c = 0; c < A.cols(); c++) {
            ASSERT_FLOAT_EQ(A(r, c), A3d(r, c));
        }
    }
    for (size_t i = 0; i < A.cols(); i++) {
        ASSERT_FLOAT_EQ(0, A3d(2, i));
    }
}
