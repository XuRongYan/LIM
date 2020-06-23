//
// Created by 徐溶延 on 2020/6/14.
//

#include "ARAPEnergy.h"
#include <fstream>
#include "../utils/EigenUtils.h"

namespace xry_mesh {
    ARAPEnergy::ARAPEnergy() {}

    ARAPEnergy::ARAPEnergy(const Eigen::VectorXf &x,
                           const Eigen::Matrix2Xf &V,
                           const Eigen::Matrix3Xi &F,
                           const std::vector<float> &areas) : NormEnergy(x),
                                                              V_(V),
                                                              F_(F),
                                                              ideal_areas_(areas) {}

    float ARAPEnergy::value() {
        return value(x_);
    }

    float ARAPEnergy::value(const Eigen::VectorXf &x) {
        float error = 0;
        Eigen::VectorXf tmp_x = x;
        Eigen::Matrix2Xf V = Eigen::Map<Eigen::Matrix2Xf>(tmp_x.data(), 2, tmp_x.size() / 2);
        for (size_t i = 0; i < vec_R_.size(); i++) {
            Eigen::Matrix2f m;
            for (size_t j = 0; j < 2; j++) {
                m.col(j) = (V.col(F_(j + 1, i)) - V.col(F_(0, i)));
            }
            error += ideal_areas_[i] * (m * inv_deltas_[i] - vec_R_[i]).squaredNorm();
        }

        return error;
    }

    void ARAPEnergy::init() {
        assert(ideal_areas_.size() == F_.cols());
        computeIdeals();
        computeInvDetails();
        computeA();
        localPhase();
        computeB();
    }

    void ARAPEnergy::update(const Eigen::VectorXf &x) {
        x_ = x;
        V_ = Eigen::Map<Eigen::Matrix2Xf>(x_.data(), 2, x.rows() / 2);
        localPhase();
        computeB();
    }

    void ARAPEnergy::computeA() {
        A_ = computeGradientOperator();
    }

    void ARAPEnergy::computeB() {
        b_ = recomposeR();
    }

    void ARAPEnergy::localPhase() {
        vec_R_.clear();
        vec_R_.resize(F_.cols());
        Eigen::Matrix2f m;
        for (size_t i = 0; i < F_.cols(); i++) {
            for (size_t j = 0; j < 2; j++) {
                m.col(j) = (V_.col(F_(j + 1, i)) - V_.col(F_(0, i)));
            }
            optimizeRotation(m * inv_deltas_[i], vec_R_[i]);
        }
    }

    void ARAPEnergy::optimizeRotation(const Eigen::Matrix2f &J, Eigen::Matrix2f &R) {
        const Eigen::Vector2d data(J(0, 0) + J(1, 1), J(0, 1) - J(1, 0));
        R << data[0], data[1], -data[1], data[0];
        const double r = data.norm();
        if (r < 1e-6) {
            R.setIdentity();
        } else {
            R /= r;
        }
    }

    Eigen::Matrix<float, 2, 3>
    ARAPEnergy::flattenTriangle(const Eigen::Vector2f &p1, const Eigen::Vector2f &p2, const Eigen::Vector2f &p3) {
        const double aa = (p2 - p1).squaredNorm();
        const double bb = (p3 - p2).squaredNorm();
        const double cc = (p1 - p3).squaredNorm();

        const double a = std::sqrt(aa);
        const double c = std::sqrt(cc);

        const double angle = std::acos((aa + cc - bb) / (2 * a * c));

        Eigen::Matrix<float, 2, 3> p2d;
        p2d << 0, a, c * std::cos(angle),
                0, 0, c * std::sin(angle);
        return p2d;
    }

    Eigen::Matrix<float, 2, 3>
    ARAPEnergy::useOriginalTriangles(const Eigen::Vector2f &p1, const Eigen::Vector2f &p2, const Eigen::Vector2f &p3) {
        Eigen::Matrix<float, 2, 3> p2d;
        p2d << p1.x(), p2.x(), p3.x(),
                p1.y(), p2.y(), p3.y();
        return p2d;
    }

    void ARAPEnergy::computeIdeals() {
        ideal_elems_.clear();
        for (size_t i = 0; i < F_.cols(); i++) {
            Eigen::Matrix<float, 2, 3> ideal = useOriginalTriangles(V_.col(F_(0, i)),
                                                                    V_.col(F_(1, i)),
                                                                    V_.col(F_(2, i)));
            ideal_elems_.emplace_back(ideal);
        }
    }

    void ARAPEnergy::computeInvDetails() {
        inv_deltas_.clear();
        assert(ideal_elems_.size() == F_.cols());
        for (size_t i = 0; i < F_.cols(); i++) {
            Eigen::Matrix2f m;
            Eigen::Matrix<float, 2, 3> ideal = ideal_elems_[i];
            for (size_t j = 0; j < 2; j++) {
                m.col(j) = ideal.col(j + 1) - ideal.col(0);
            }
            inv_deltas_.emplace_back(m.inverse());
        }
    }

    Eigen::VectorXf ARAPEnergy::recomposeR() {
        Eigen::VectorXf R(4 * F_.cols());
        R.setZero();
        for (size_t i = 0; i < F_.cols(); i++) {
            recomposeRi(R, i);
        }
        return R;
    }

    void ARAPEnergy::recomposeRi(Eigen::VectorXf &R, size_t i) {
        for (size_t r = 0; r < 2; r++) {
            for (size_t c = 0; c < 2; c++) {
                R[4 * i + 2 * r + c] = sqrt(ideal_areas_[i]) * vec_R_[i](c, r);
            }
        }
    }

    Eigen::SparseMatrix<float> ARAPEnergy::computeGradientOperator() {
        std::vector<Eigen::Triplet<float>> triplets;
        for (size_t i = 0, col = 0; i < F_.cols(); i++) {
            for (size_t j = 0; j < 2; j++, col++) {
                float sum = 0;
                for (size_t k = 0; k < 2; k++) {
                    const float val = sqrt(ideal_areas_[i]) * inv_deltas_[i](k, j);
                    sum -= val;
                    for (size_t t = 0; t < 2; t++) {
                        triplets.emplace_back(2 * F_(k + 1, i) + t, 2 * col + t, val);
                    }
                }
                for (size_t t = 0; t < 2; t++) {
                    triplets.emplace_back(2 * F_(0, i) + t, 2 * col + t, sum);
                }
            }
        }
        Eigen::SparseMatrix<float> Gt(2 * V_.cols(), 4 * F_.cols());
        Gt.setFromTriplets(triplets.begin(), triplets.end());
        return Gt.transpose();
    }

    const Eigen::Matrix2Xf &ARAPEnergy::getV() const {
        return V_;
    }

    const Eigen::Matrix3Xi &ARAPEnergy::getF() const {
        return F_;
    }

    const std::vector<float> &ARAPEnergy::getIdealAreas() const {
        return ideal_areas_;
    }

    const std::vector<Eigen::Matrix2f> &ARAPEnergy::getVecR() const {
        return vec_R_;
    }

    const std::vector<Eigen::Matrix<float, 2, 3>> &ARAPEnergy::getIdealElems() const {
        return ideal_elems_;
    }

    const std::vector<Eigen::Matrix2f> &ARAPEnergy::getInvDeltas() const {
        return inv_deltas_;
    }

} // namespace xry_mesh