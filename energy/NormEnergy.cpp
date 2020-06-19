//
// Created by 徐溶延 on 2020/6/14.
//

#include "NormEnergy.h"

namespace xry_mesh {
    NormEnergy::NormEnergy() {}

    NormEnergy::NormEnergy(const Eigen::VectorXf &x) : BaseEnergy(x) {}

    float NormEnergy::value() {
        return (A_ * x_ - b_).squaredNorm();
    }

    float NormEnergy::value(const Eigen::VectorXf &x) {
        assert(x.rows() == A_.cols());
        return (A_ * x - b_).squaredNorm();
    }

    void NormEnergy::init() {
        computeA();
        computeB();
        assert(A_.cols() == b_.rows());
    }

    Eigen::VectorXf NormEnergy::jacobian() const {
        Eigen::SparseMatrix<float> At = A_.transpose();
        return 2 * At * A_ * x_ - 2 * At * b_;
    }

    Eigen::SparseMatrix<float> NormEnergy::hessian() const {
        Eigen::SparseMatrix<float> At = A_.transpose();
        return 2 * At * A_;
    }

    const Eigen::SparseMatrix<float> &NormEnergy::getA() const {
        return A_;
    }

    const Eigen::VectorXf &NormEnergy::getB() const {
        return b_;
    }

    void NormEnergy::update(const Eigen::VectorXf &x) {
        BaseEnergy::update(x);
        computeA();
        computeB();
        assert(A_.cols() == b_.rows());
    }
}
