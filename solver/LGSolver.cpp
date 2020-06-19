//
// Created by 徐溶延 on 2020/6/15.
//

#include "LGSolver.h"
#include "../utils/SurfaceMeshUtils.h"
#include "../utils/EigenUtils.h"

LGSolver::LGSolver() {}

LGSolver::LGSolver(const Surface_Mesh::SurfaceMesh &mesh,
                   const std::vector<std::pair<int, Eigen::VectorXf>> &posConstrains,
                   size_t maxIter) : max_iter(maxIter),
                                     mesh_(mesh),
                                     pos_constrains_(posConstrains) {

}

void LGSolver::init() {
    assert(!mesh_.empty());
    Eigen::Matrix3Xf V3d = xry_mesh::getGeometryMatrix<float>(mesh_);
    V_ = V3d.block(0, 0, 2, V3d.cols());
    F_ = xry_mesh::getTopoMatrix(mesh_);
    x_ = xry_mesh::vt2v<float>(V_.transpose());
    xry_mesh::computeAreas<float>(mesh_);
    areas = xry_mesh::faceProperty2StdVector<float>(mesh_, "f:areas");

    arap_energy_ = xry_mesh::ARAPEnergy(x_, V_, F_, areas);
    arap_energy_.init();

    pos_energy_ = xry_mesh::PosEnergy(x_, 2, pos_constrains_, 0);
    pos_energy_.init();
    arap_A = arap_energy_.getA();
    pos_A = pos_energy_.getA();
    L_ = arap_A.transpose() * arap_A + alpha * pos_A;
    llt_.compute(L_);
    b_.resize(x_.size());
    assert(llt_.info() == Eigen::Success);
}

Eigen::VectorXf LGSolver::solve() {
    float err0 = -1, err1 = 0;
    size_t iter = 0;
    dbg("start solve...");
    while (fabs(err1 - err0) > 1e-6 && iter < max_iter) {
        err0 = err1;
        subSolve();
        err1 = computeError();
        dbg("..............................................");
        dbg(iter);
        dbg(err1);
        iter++;
    }
    return x_;
}

void LGSolver::subSolve() {
    b_.setZero();
    b_ = arap_A.transpose() * arap_energy_.getB() + alpha * pos_energy_.getB();
    x_ = llt_.solve(b_);
    arap_energy_.update(x_);
    pos_energy_.update(x_);
}

float LGSolver::computeError() {
    const float arap_error = arap_energy_.value();
    dbg(arap_error);
    const float pos_error = pos_energy_.value();
    dbg(pos_error);
    return arap_error + pos_error;
}
