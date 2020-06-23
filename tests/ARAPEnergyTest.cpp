//
// Created by 徐溶延 on 2020/6/15.
//
#include <SurfaceMesh/SurfaceMesh.h>
#include <SurfaceMesh/IO.h>
#include <gtest/gtest.h>

#include "../energy/ARAPEnergy.h"
#include "../utils/SurfaceMeshUtils.h"
#include "../utils/EigenUtils.h"
#include "../utils/FileIOUtils.h"

class ARAPEnergyTest : public ::testing::Test {
protected:
    void SetUp() override {
        Test::SetUp();
        Surface_Mesh::read_obj(mesh_, "plane_2d_2_2.obj");
        Eigen::Matrix3Xf V = xry_mesh::getGeometryMatrix<float>(mesh_);
        Eigen::MatrixX2f v2d = V.transpose().block(0, 0, V.cols(), 2);
        Eigen::Matrix3Xi F = xry_mesh::getTopoMatrix(mesh_);
        Eigen::VectorXf x = xry_mesh::vt2v<float>(v2d);
        xry_mesh::computeAreasDet<float>(mesh_);
        std::vector<float> areas = xry_mesh::faceProperty2StdVector<float>(mesh_, "f:areas");
        arapEnergy = xry_mesh::ARAPEnergy(x, v2d.transpose(), F, areas);
        arapEnergy.init();
    }

protected:
    Surface_Mesh::SurfaceMesh mesh_;
    xry_mesh::ARAPEnergy arapEnergy;
};

TEST_F(ARAPEnergyTest, rotateMatrixTest) {
    std::vector<Eigen::Matrix2f> vec_R = arapEnergy.getVecR();
    for (const auto &R : vec_R) {
        EXPECT_FLOAT_EQ(1.0, R.determinant());
    }
    Eigen::VectorXf x(2 * mesh_.n_vertices());
    x << 0, 0, 0, 1, 0, 2, 1, 0, 1, 1, 1, 2, 2, 0, 2, 1, 4, 4;
    arapEnergy.update(x);
    ASSERT_EQ(x, arapEnergy.getX());
    for (const auto &R : vec_R) {
        EXPECT_FLOAT_EQ(1.0, R.determinant());
    }
}

TEST_F(ARAPEnergyTest, valueTest) {
    float val = arapEnergy.value();
    ASSERT_FLOAT_EQ(val, 0);
}

TEST_F(ARAPEnergyTest, numericalJacobianTest) {
	Eigen::VectorXf J = arapEnergy.jacobian();
	Eigen::VectorXf J_num = arapEnergy.numericalJacobian(1e-6);
	const float err = (J - J_num).squaredNorm();
	dbg(err);
	ASSERT_TRUE(err < 1e-5);
}