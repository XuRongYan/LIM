//
// Created by 徐溶延 on 2020/5/30.
//

#include "SurfaceMeshUtils.h"
#include <Eigen/Dense>

namespace xry_mesh {
    Eigen::Matrix3Xi getTopoMatrix(const Surface_Mesh::SurfaceMesh &mesh) {
        std::vector<int> face_id;
        for (const auto &f : mesh.faces()) {
            for (const auto &v : mesh.vertices(f)) {
                face_id.push_back(v.idx());
            }
        }
        return Eigen::Map<Eigen::Matrix3Xi>(face_id.data(), 3, face_id.size() / 3);
    }

    void ensurePositiveAreas(Surface_Mesh::SurfaceMesh &mesh) {
        for (const auto &f : mesh.faces()) {
            std::vector<Surface_Mesh::Point> points;
            std::vector<Surface_Mesh::SurfaceMesh::Vertex> vertices;
            for (auto v : mesh.vertices(f)) {
                vertices.emplace_back(v);
                points.emplace_back(mesh.position(v));
            }
            if (computeAreaDet(points[0], points[1], points[2]) < 0) {
                std::swap(vertices[1], vertices[2]);
            }
        }
    }

} // namespace xry_mesh