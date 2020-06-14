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

} // namespace xry_mesh