//
// Created by 徐溶延 on 2020/6/14.
//
#include <SurfaceMesh/IO.h>
#include <gtest/gtest.h>

#include "../utils/SurfaceMeshUtils.h"

class SurfaceMeshUtilsTest : public ::testing::Test {
public:
protected:
    void SetUp() override {
        Test::SetUp();
        Surface_Mesh::read_obj(mesh_, "plane_2d_2_2.obj");
    }

public:
    Surface_Mesh::SurfaceMesh mesh_;
};

TEST_F(SurfaceMeshUtilsTest, computeAreaTest) {
    xry_mesh::computeAreas<Surface_Mesh::Scalar>(mesh_);
    //判断是否成功添加面积条目
    bool exist = false;
    for (const auto &key : mesh_.face_properties()) {
        if (key == "f:areas") exist = true;
    }
    ASSERT_TRUE(exist);

    //判断面积大小是否计算正确
    auto areas = mesh_.get_face_property<Surface_Mesh::Scalar>("f:areas");
    for (const auto &f : mesh_.faces()) {
        ASSERT_FLOAT_EQ(areas[f], 0.5);
    }

    Eigen::Vector2f p1, p2, p3;
    p1 << 0, 0;
    p2 << 1, 0;
    p3 << 0, 1;
    ASSERT_EQ(xry_mesh::computeArea2D(p1, p2, p3), 0.5);
}

TEST_F(SurfaceMeshUtilsTest, getGeometryMatrixTest) {
    Eigen::Matrix3Xf V = xry_mesh::getGeometryMatrix<float>(mesh_);
    //判断size
    ASSERT_EQ(V.rows(), 3);
    ASSERT_EQ(V.cols(), mesh_.n_vertices());
    //判断数值
    Eigen::Matrix3Xf VTest(3, mesh_.n_vertices());
    VTest << 0, 0, 0, 1, 1, 1, 2, 2, 2,
            0, 1, 2, 0, 1, 2, 0, 1, 2,
            0, 0, 0, 0, 0, 0, 0, 0, 0;
    for (size_t r = 0; r < 3; r++) {
        for (size_t c = 0; c < mesh_.n_vertices(); c++) {
            ASSERT_FLOAT_EQ(V(r, c), VTest(r, c));
        }
    }
}

TEST_F(SurfaceMeshUtilsTest, getTopoMatrixTest) {
    Eigen::Matrix3Xi F = xry_mesh::getTopoMatrix(mesh_);
    //判断size
    ASSERT_EQ(F.rows(), 3);
    ASSERT_EQ(F.cols(), mesh_.n_faces());
    //判断数值
    Eigen::Matrix3Xf FTest(3, mesh_.n_faces());
    FTest << 0, 1, 1, 2, 3, 4, 4, 5,
            3, 3, 4, 4, 6, 6, 7, 7,
            1, 4, 2, 5, 4, 7, 5, 8;
    for (size_t r = 0; r < 3; r++) {
        for (size_t c = 0; c < mesh_.n_faces(); c++) {
            ASSERT_FLOAT_EQ(F(r, c), FTest(r, c));
        }
    }
}

TEST_F(SurfaceMeshUtilsTest, faceProperty2StdVectorTest) {
    xry_mesh::computeAreas<Surface_Mesh::Scalar>(mesh_);
    std::vector<Surface_Mesh::Scalar> areas = xry_mesh::faceProperty2StdVector<Surface_Mesh::Scalar>(mesh_, "f:areas");
    //判断面积大小是否计算正确
    for (const auto &area : areas) {
        ASSERT_FLOAT_EQ(area, 0.5);
    }
}

TEST_F(SurfaceMeshUtilsTest, ensurePositiveAreasTest) {
    xry_mesh::ensurePositiveAreas(mesh_);
    Surface_Mesh::write_obj(mesh_, "ensure_positive_test.obj");
    for (const auto &f : mesh_.faces()) {
        std::vector<Surface_Mesh::Point> points;
        for (auto v : mesh_.vertices(f)) {
            points.emplace_back(mesh_.position(v));
        }
        ASSERT_TRUE(xry_mesh::computeAreaDet(points[0], points[1], points[2]) > 0);
    }
}
