//
// Created by 徐溶延 on 2020/6/14.
//
#include <SurfaceMesh/SurfaceMesh.h>
#include <SurfaceMesh/IO.h>
#include <gtest/gtest.h>

#include "../energy/PosEnergy.h"
#include "../utils/SurfaceMeshUtils.h"
#include "../utils/EigenUtils.h"
#include "../utils/FileIOUtils.h"

class PosEnergyTest : public ::testing::Test {
protected:
    void SetUp() override {
        Test::SetUp();
        Surface_Mesh::read_obj(mesh_, "plane_2d_2_2.obj");
        Eigen::Matrix3Xf V = xry_mesh::getGeometryMatrix<float>(mesh_);
        Eigen::MatrixX2f v2d = V.transpose().block(0, 0, V.cols(), 2);
        Eigen::VectorXf x = xry_mesh::vt2v<float>(v2d);
        std::vector<std::pair<int, Eigen::VectorXf>> pos_constrains
                = xry_mesh::readPosFile<float>("pos_2_2.pts", 2);
        posEnergy = xry_mesh::PosEnergy(x, 2, pos_constrains, 0);
        posEnergy.init();
    }

protected:
    Surface_Mesh::SurfaceMesh mesh_;
    xry_mesh::PosEnergy posEnergy;
};

TEST_F(PosEnergyTest, PosValueTest) {
    Eigen::VectorXf v(2 * mesh_.n_vertices());
    v << 0, 0, 0, 1, 0, 2, 1, 0, 1, 1, 1, 2, 2, 0, 2, 1, 2, 2;
    float value = posEnergy.value();
    ASSERT_EQ(v, posEnergy.getX());
    ASSERT_FLOAT_EQ(0, value);
}