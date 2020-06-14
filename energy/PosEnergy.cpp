//
// Created by 徐溶延 on 2020/6/14.
//

#include "PosEnergy.h"

namespace xry_mesh {

    PosEnergy::PosEnergy() {

    }

    void PosEnergy::init() {
        NormEnergy::init();
        if (enable_dbg) {
            dbg(x_);
        }
    }

    void PosEnergy::computeA() {
        A_.resize(x_.size(), x_.size());
        std::vector<Eigen::Triplet<float>> triplets;
        for (const auto &pair : pos_constrains_) {
            for (size_t i = 0; i < dim; i++) {
                triplets.emplace_back(dim * pair.first + i, dim * pair.first + i, 1);
            }
        }
        A_.setFromTriplets(triplets.begin(), triplets.end());
        if (enable_dbg) {
            dbg(A_);
        }
    }

    void PosEnergy::computeB() {
        b_ = computeDi();
        if (enable_ti) {
            b_ = computeTi(b_);
        }
        if (enable_dbg) {

            dbg(b_);
        }
    }

    void PosEnergy::update(const Eigen::VectorXf &x) {
        NormEnergy::update(x);
    }

    PosEnergy::PosEnergy(const Eigen::VectorXf &x,
                         size_t dim,
                         const std::vector<std::pair<int, Eigen::VectorXf>> &posConstrains,
                         float mu) : NormEnergy(x),
                                     dim(dim),
                                     pos_constrains_(posConstrains),
                                     mu_(mu) {
        assert(!posConstrains.empty());
        assert(posConstrains[0].second.size() == dim);
    }

    Eigen::VectorXf PosEnergy::computeDi() {
        Eigen::VectorXf d_i(x_.size());
        d_i.setZero();
        for (const auto &ele : pos_constrains_) {
            for (size_t i = 0; i < dim; i++) {
                d_i[dim * ele.first + i] = ele.second[i];
            }
        }
        return d_i;
    }

    Eigen::VectorXf PosEnergy::computeTi(const Eigen::VectorXf &d_i) {
        assert(A_.outerSize() == dim * pos_constrains_.size());
        return A_ * x_ + 1.0 / (1 + std::pow(mu_, 2)) * (d_i - A_ * x_);
    }
} // namespace xry_mesh