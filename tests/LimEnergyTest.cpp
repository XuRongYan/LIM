//
// Created by 徐溶延 on 2020/6/22.
//
#include <gtest/gtest.h>
#include <SurfaceMesh/SurfaceMesh.h>
#include <SurfaceMesh/IO.h>

#include "../energy/LimEnergy.h"
#include "../energy/BarrierEnergy.h"
#include "../energy/ARAPEnergy.h"
#include "../utils/SurfaceMeshUtils.h"
#include "../utils/EigenUtils.h"
#include "../utils/FileIOUtils.h"

class LimEnergyTest : public ::testing::Test {
protected:
	void SetUp() override {
		Test::SetUp();
		Surface_Mesh::read_obj(mesh_, "plane_2d_4_4.obj");
		std::vector<std::pair<int, Eigen::VectorXf>> pos_constrain1 = xry_mesh::readPosFile<float>("pos_4_4.pts", 2);
		limEnergy = xry_mesh::LimEnergy(mesh_, pos_constrain1);
		limEnergy.init();
	}

protected:
	Surface_Mesh::SurfaceMesh mesh_;
	xry_mesh::LimEnergy limEnergy;
};

TEST_F(LimEnergyTest, valueTest) {
	float val = limEnergy.value();
	float val_x = limEnergy.value(limEnergy.getX());
	ASSERT_FLOAT_EQ(val, val_x);
}

TEST_F(LimEnergyTest, numericalJacobianTest) {
	Eigen::VectorXf J = limEnergy.jacobian();
	Eigen::VectorXf J_num = limEnergy.numericalJacobian(1e-4);
	//dbg(J);
	//dbg(J_num);
	const float err = (J - J_num).squaredNorm();
	dbg(err);
	EXPECT_TRUE(err < 1e-2);
}

TEST_F(LimEnergyTest, arapNumericalJacobianTest) {
	limEnergy.setEnableBarrierFunc(false);
	Eigen::VectorXf J = limEnergy.jacobian();
	Eigen::VectorXf J_num = limEnergy.numericalJacobian(1e-4);
	const float err = (J - J_num).squaredNorm();
	dbg(err);
	EXPECT_TRUE(err < 1e-2);
}

TEST_F(LimEnergyTest, paraNumericalJacobianTest) {
	limEnergy.setAlpha(1);
	limEnergy.setBeta(1);
	limEnergy.setEnableBarrierFunc(false);
	Eigen::VectorXf J = limEnergy.jacobian();
	Eigen::VectorXf J_num = limEnergy.numericalJacobian(1e-6);
	const float err = (J - J_num).squaredNorm();
	dbg(err);
	EXPECT_TRUE(err < 1e-2);
}