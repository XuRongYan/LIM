#include "arap.h"
#include <cmath>
#include <iostream>
#include <Eigen/Sparse>

namespace jy_mesh {
    int optimal_rotation(const Eigen::Matrix2d &J, Eigen::Matrix2d &R) {
        const Eigen::Vector2d data(J(0, 0) + J(1, 1), J(0, 1) - J(1, 0));
        R << data[0], data[1], -data[1], data[0];
        const double r = data.norm();
        if (r < 1e-6) {
            R.setIdentity();
        } else {
            R /= r;
        }
        return 1;
    }

    double flatten_triangle(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
                            const Eigen::Vector2d &p2, Eigen::Matrix<double, 2, 3> &p2d) {

        const double aa = (p1 - p0).squaredNorm();
        const double bb = (p2 - p1).squaredNorm();
        const double cc = (p0 - p2).squaredNorm();

        const double a = std::sqrt(aa);
        const double c = std::sqrt(cc);

        const double angle = std::acos((aa + cc - bb) / (2 * a * c));

        p2d << 0, a, c * std::cos(angle),
                0, 0, c * std::sin(angle);

        return a * p2d(1, 2);
    }

    void init_para_info(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
                        Eigen::VectorXd &sqrt_areas, std::vector<Eigen::Matrix2d> &inv,
                        Eigen::SparseMatrix<double> &Gt) {

        Eigen::Matrix2d m;
        Eigen::Matrix<double, 2, 3> t;
        sqrt_areas.resize(F.cols());
        inv.reserve(F.cols());
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(F.cols() * 6);
        for (size_t i = 0, col = 0; i < F.cols(); ++i) {
            const double A = flatten_triangle(V.col(F(0, i)), V.col(F(1, i)), V.col(F(2, i)), t);
            for (size_t j = 0; j < 2; ++j) {
                m.col(j) = t.col(j + 1) - t.col(0);
            }
            sqrt_areas[i] = sqrt(A);
            inv.push_back(m.inverse());
            for (size_t j = 0; j < 2; ++j, ++col) {
                double sum = 0;
                for (size_t k = 0; k < 2; ++k) {
                    const double val = inv[i](k, j) * sqrt_areas[i];
                    sum -= val;
                    triplets.push_back(Eigen::Triplet<double>(F(k + 1, i), col, val));
                }
                triplets.push_back(Eigen::Triplet<double>(F(0, i), col, sum));
            }

        }
        //build the matrix
        Gt.resize(V.cols(), 2 * F.cols());
        Gt.setFromTriplets(triplets.begin(), triplets.end());
    }


    int arap_deformation(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F,
                         const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
                         Eigen::Matrix2Xd &Vout, size_t max_iter) {
        Eigen::MatrixX2d x = V.transpose();

        //initilize
        Eigen::VectorXd sqrt_areas, areas(F.cols());
        std::vector<Eigen::Matrix2d> R(F.cols()), inv;
        Eigen::SparseMatrix<double> Gt;
        init_para_info(V, F, sqrt_areas, inv, Gt);
        for (size_t i = 0; i < sqrt_areas.size(); ++i) {
            areas[i] = sqrt_areas[i] * sqrt_areas[i];
        }
        Eigen::SparseMatrix<double> L = Gt * Gt.transpose();
        const double w = 1e8;
        for (size_t i = 0; i < bc.size(); ++i) {
            L.coeffRef(bc[i].first, bc[i].first) += w;
        }
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(L);

        size_t iter = 0;
        double err0 = -1, err1 = 0;
        Eigen::Matrix2d J(2, 2);
        Eigen::Matrix2d m;
        Eigen::MatrixX2d b(2 * F.cols(), 2);
        while (fabs(err1 - err0) > 1e-6 && iter < max_iter) {
            err0 = err1;
            //local_step(nodes);
            for (size_t i = 0; i < F.cols(); ++i) {
                for (size_t j = 0; j < 2; ++j) {
                    m.col(j) = x.row(F(j + 1, i)) - x.row(F(0, i));
                }
                optimal_rotation(m * inv[i], J);
                R[i] = J;
            }
            //global_step(nodes);
            for (size_t i = 0; i < F.cols(); ++i) {
                for (size_t j = 0; j < 2; ++j) {
                    b.row(2 * i + j) = sqrt_areas[i] * R[i].col(j);
                }
            }
            Eigen::MatrixX2d Gtb = Gt * b;
            for (size_t i = 0; i < bc.size(); ++i) {
                Gtb.row(bc[i].first) += w * bc[i].second;
            }
            x = llt.solve(Gtb); //GtGx = Gtb
            //err1 = energy_value(nodes);
            err1 = 0;
            for (size_t i = 0; i < F.cols(); ++i) {
                for (size_t j = 0; j < 2; ++j) {
                    m.col(j) = x.row(F(j + 1, i)) - x.row(F(0, i));
                }
                err1 += areas[i] * (m * inv[i] - R[i]).squaredNorm();
            }
            std::cout << ++iter << ":" << err1 << std::endl;
        }
        Vout = x.transpose();
        return 0;
    }
}