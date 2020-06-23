//
// Created by 徐溶延 on 2020/6/22.
//

#include <gtest/gtest.h>
#include <SurfaceMesh/SurfaceMesh.h>
#include <SurfaceMesh/IO.h>
#include "../energy/BarrierEnergy.h"

class BarrierEnergyTest : public ::testing::Test {
protected:
    void SetUp() override {
        Test::SetUp();
        Surface_Mesh::read_obj(mesh_, "plane_2d_1_1.obj");
        Eigen::Matrix3Xf V = xry_mesh::getGeometryMatrix<float>(mesh_);
        Eigen::MatrixX2f v2d = V.transpose().block(0, 0, V.cols(), 2);
        Eigen::Matrix3Xi F = xry_mesh::getTopoMatrix(mesh_);
        Eigen::VectorXf x = xry_mesh::vt2v<float>(v2d);
        xry_mesh::computeAreasDet<float>(mesh_);
        std::vector<float> areas = xry_mesh::faceProperty2StdVector<float>(mesh_, "f:areas");
        barrierEnergy = xry_mesh::BarrierEnergy(x, v2d.transpose(), F, areas);
        barrierEnergy.init();
    }

protected:
    Surface_Mesh::SurfaceMesh mesh_;
    xry_mesh::BarrierEnergy barrierEnergy;
};

TEST_F(BarrierEnergyTest, phiTest) {
    ASSERT_FLOAT_EQ(barrierEnergy.phi(0), std::numeric_limits<float>::infinity());
    ASSERT_FLOAT_EQ(barrierEnergy.phi(barrierEnergy.getSJ()), 0);
    ASSERT_FLOAT_EQ(barrierEnergy.phi(50), 0);
}

TEST_F(BarrierEnergyTest, gTest) {
    ASSERT_FLOAT_EQ(barrierEnergy.g(0), 0);
    ASSERT_FLOAT_EQ(barrierEnergy.g(barrierEnergy.getSJ()), 1);
    ASSERT_FLOAT_EQ(barrierEnergy.gPrime(barrierEnergy.getSJ()), 0);
    ASSERT_FLOAT_EQ(barrierEnergy.gPrime2(barrierEnergy.getSJ()), 0);
}

TEST_F(BarrierEnergyTest, weight1Test) {
    auto w1 = barrierEnergy.getWeight1();
    for (auto w : w1) {
        ASSERT_FLOAT_EQ(w, -0.9796198);
    }
}

TEST_F(BarrierEnergyTest, weight2Test) {
    auto w2 = barrierEnergy.getWeight2();
    for (auto w : w2) {
        ASSERT_FLOAT_EQ(w, 5.59782917);
    }
}

TEST_F(BarrierEnergyTest, vertIdTest) {
    Eigen::Matrix<int, 6, 2> ids;
    Eigen::Matrix<float, 6, 2> v_ids, v_test_ids;
    Eigen::VectorXf x_test = barrierEnergy.getX();
    ids << 0, 2,
            1, 3,
            4, 4,
            5, 5,
            2, 6,
            3, 7;

    v_test_ids << x_test[ids(0, 0)], x_test[ids(0, 1)],
                  x_test[ids(1, 0)], x_test[ids(1, 1)],
                  x_test[ids(2, 0)], x_test[ids(2, 1)],
                  x_test[ids(3, 0)], x_test[ids(3, 1)],
                  x_test[ids(4, 0)], x_test[ids(4, 1)],
                  x_test[ids(5, 0)], x_test[ids(5, 1)];

    v_ids << 0, 0,
             0, 1,
             1, 1,
             0, 0,
             0, 1,
             1, 1;
    ASSERT_EQ(ids, barrierEnergy.getTriangleVertexIdx());
    ASSERT_EQ(v_test_ids, v_ids);
}

TEST_F(BarrierEnergyTest, valueTest) {
    Eigen::VectorXf my_x(8);
    my_x << 0, 0, 0, 1, 1, 0, 1, 1;
    float val1 = barrierEnergy.value();
    float val2 = barrierEnergy.value(barrierEnergy.getX());
    ASSERT_FLOAT_EQ(val1, 0.28572408);
    ASSERT_FLOAT_EQ(val1, val2);
}

TEST_F(BarrierEnergyTest, jacobianTest) {
	Eigen::VectorXf J_test(8);
	J_test.setZero();
	J_test[0] += -0.9796198 * -0.5;
	J_test[1] += -0.9796198 * -0.5;
	J_test[4] += -0.9796198 * 0.5;
	J_test[5] += -0.9796198 * 0;
	J_test[2] += -0.9796198 * 0;
	J_test[3] += -0.9796198 * 0.5;

	J_test[2] += -0.9796198 * -0.5;
	J_test[3] += -0.9796198 * 0;
	J_test[4] += -0.9796198 * 0;
	J_test[5] += -0.9796198 * -0.5;
	J_test[6] += -0.9796198 * 0.5;
	J_test[7] += -0.9796198 * 0.5;

	Eigen::VectorXf J = barrierEnergy.jacobian();
	for (size_t i = 0; i < J.size(); i++) {
		ASSERT_FLOAT_EQ(J[i], J_test[i]);
	}
}

TEST_F(BarrierEnergyTest, lambdaGradTest) {
	Eigen::VectorXf J_test1(8), J_test2(8);
	J_test1.setZero();
	J_test1[0] = -0.5;
	J_test1[1] = -0.5;
	J_test1[4] = 0.5;
	J_test1[5] = 0;
	J_test1[2] = 0;
	J_test1[3] = 0.5;
	J_test2.setZero();
	J_test2[2] = -0.5;
	J_test2[3] = 0;
	J_test2[4] = 0;
	J_test2[5] = -0.5;
	J_test2[6] = 0.5;
	J_test2[7] = 0.5;

	std::vector<Eigen::VectorXf> vec_J = barrierEnergy.lambdaJacobian();
	for (size_t i = 0; i < vec_J[0].size(); i++) {
		ASSERT_FLOAT_EQ(vec_J[0][i], J_test1[i]);
	}

	for (size_t i = 0; i < vec_J[1].size(); i++) {
		ASSERT_FLOAT_EQ(vec_J[1][i], J_test2[i]);
	}
}

TEST_F(BarrierEnergyTest, hessianTest) {
	Eigen::SparseMatrix<float> H_test(8, 8);
	H_test.setZero();
	std::vector<Eigen::VectorXf> vec_J = barrierEnergy.lambdaJacobian();
	auto w1 = barrierEnergy.getWeight1();
	auto w2 = barrierEnergy.getWeight2();
	Eigen::SparseMatrix<float> H1(8, 8), H2(8, 8);
	H1 = w2[0] * (vec_J[0] * vec_J[0].transpose()).sparseView();
	H_test += H1;
	H2 = w2[1] * (vec_J[1] * vec_J[1].transpose()).sparseView();
	H_test += H2;
	H_test.coeffRef(0, 5) += w1[0];
	H_test.coeffRef(5, 0) += w1[0];
	H_test.coeffRef(0, 3) += -w1[0];
	H_test.coeffRef(3, 0) += -w1[0];
	H_test.coeffRef(1, 4) += -w1[0];
	H_test.coeffRef(4, 1) += -w1[0];
	H_test.coeffRef(1, 2) += w1[0];
	H_test.coeffRef(2, 1) += w1[0];
	H_test.coeffRef(4, 3) += w1[0];
	H_test.coeffRef(3, 4) += w1[0];
	H_test.coeffRef(5, 2) += -w1[0];
	H_test.coeffRef(2, 5) += -w1[0];

	H_test.coeffRef(2, 5) += w1[1];
	H_test.coeffRef(5, 2) += w1[1];
	H_test.coeffRef(2, 7) += -w1[1];
	H_test.coeffRef(7, 2) += -w1[1];
	H_test.coeffRef(3, 4) += -w1[1];
	H_test.coeffRef(4, 3) += -w1[1];
	H_test.coeffRef(3, 6) += w1[1];
	H_test.coeffRef(6, 3) += w1[1];
	H_test.coeffRef(4, 7) += w1[1];
	H_test.coeffRef(7, 4) += w1[1];
	H_test.coeffRef(5, 6) += -w1[1];
	H_test.coeffRef(6, 5) += -w1[1];

	auto H = barrierEnergy.hessian();
	ASSERT_FLOAT_EQ((H - H_test).norm(), 0);
}

TEST_F(BarrierEnergyTest, numericalJacobianTest) {
	Eigen::VectorXf J = barrierEnergy.jacobian();
	Eigen::VectorXf J_num = barrierEnergy.numericalJacobian(1e-6);
	const float err = (J - J_num).squaredNorm();
	dbg(err);
	ASSERT_TRUE(err < 1e-2);
}

