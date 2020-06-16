//
// Created by 徐溶延 on 2020/6/15.
//

#include "BarrierEnergy.h"

namespace xry_mesh {
    BarrierEnergy::BarrierEnergy() {}

    BarrierEnergy::BarrierEnergy(const Eigen::VectorXf &x,
                                 const Eigen::Matrix2Xf &V,
                                 const Eigen::Matrix3Xi &F) : BaseEnergy(x),
                                                              V_(V),
                                                              F_(F) {}

    float BarrierEnergy::value() {
        return 0;
    }

    float BarrierEnergy::value(const Eigen::VectorXf &x) {
        return 0;
    }

    void BarrierEnergy::init() {
        assert(!areas.empty());
        assert(F_.cols() == areas.size());
        computeW1();
        computeW2();
        computeEpsilon();
        computeVertexIdxMatrix();
    }

    Eigen::VectorXf BarrierEnergy::jacobian() const {
        return Eigen::VectorXf();
    }

    Eigen::SparseMatrix<float> BarrierEnergy::hessian() const {
        return Eigen::SparseMatrix<float>();
    }

    void BarrierEnergy::update(const Eigen::VectorXf &x) {
        // TODO 每次更新时更新三角形面积
        BaseEnergy::update(x);
    }

    void BarrierEnergy::computeW1() {
        weight1.clear();
        weight1.resize(F_.cols());
        for (size_t i = 0; i < F_.cols(); i++) {
            weight1[i] = -gPrime(areas[i]) / std::pow(g(areas[i]), 2);
        }
    }

    void BarrierEnergy::computeW2() {
        weight2.clear();
        weight2.resize(F_.cols());
        for (size_t i = 0; i < F_.cols(); i++) {
            weight2[i] = 2 * std::pow(gPrime(areas[i]), 2)
                         - gPrime2(areas[i]) * gPrime(areas[i]) / std::pow(g(areas[i]), 3);
        }
    }

    void BarrierEnergy::computeEpsilon() {
        float MIN = std::numeric_limits<float>::infinity();
        for (size_t i = 0; i < F_.cols(); i++) {
            auto area = areas[i];
            MIN = std::min(MIN, area);
        }
        epsilon_ = 1e-5 * MIN;
        if (enable_dbg) {
            dbg(epsilon_);
        }
    }

    void BarrierEnergy::computeVertexIdxMatrix() {
        triangle_vertex_idx_.resize(6, F_.cols());
        for (size_t i = 0; i < F_.cols(); i++) {
            for (size_t j = 0; j < F_.rows(); i++) {
                size_t row = 0;
                for (size_t k = 0; k < 2; k++, row++) {
                    triangle_vertex_idx_(row, i) = 2 * F_(j, i) + k;
                }
            }
        }
    }

} // namespace xry_mesh