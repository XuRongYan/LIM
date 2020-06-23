//
// Created by 徐溶延 on 2020/6/16.
//

#include "LMSolver.h"


LMSolver::LMSolver() {}

LMSolver::LMSolver(const Surface_Mesh::SurfaceMesh &mesh,
                   const std::vector<std::pair<int, Eigen::VectorXf>> &posConstrains,
                   size_t maxIter) : max_iter_(maxIter),
                                     mesh_(mesh),
                                     pos_constrains_(posConstrains) {}

void LMSolver::init() {
    lim_energy_ = xry_mesh::LimEnergy(mesh_, pos_constrains_);
    lim_energy_.setAlpha(alpha_);
    lim_energy_.setBeta(beta_);
    lim_energy_.setSigma(sigma_);
    lim_energy_.setR(r_);
    lim_energy_.setT(t_);
    lim_energy_.setSJ(s_j_);
    lim_energy_.setEnableBarrierFunc(enable_barrier_func_);
    lim_energy_.setEnableUpdateAlpha(enable_update_alpha_);
    lim_energy_.init();
    lim_energy_.setMu(mu_);
    J_ = lim_energy_.jacobian();
    H_ = lim_energy_.hessian();
    x_ = lim_energy_.getX();

}

Eigen::VectorXf LMSolver::solve() {
    float err0 = -1, err1 = 0;
    size_t iter = 0;
    dbg("start solve...");
    while (fabs(err1 - err0) > 1e-6 && iter < max_iter_) {
        err0 = err1;
        x_ = subSolve();
        err1 = computeError();
        dbg("..............................................");
        dbg(iter);
        dbg(err1);
        iter++;
    }
    return x_;
}

Eigen::VectorXf LMSolver::subSolve() {

    if (enable_mu_search_) {
        mu_ = lineSearchMu(mu_);
        dbg(mu_);
        lim_energy_.setMu(mu_);
    }


    Eigen::VectorXf p_i;
    p_i = solvePi(p_i);
    if (enable_sigma_search_) {
        sigma_ = lineSearchSigma(p_i, sigma_);
        dbg(sigma_);
        lim_energy_.setSigma(sigma_);
    }

    auto x_i = x_ - sigma_ * p_i;
    lim_energy_.update(x_i);
    J_ = lim_energy_.jacobian();
    H_ = lim_energy_.hessian();
    //compareNumerical(x_, x_i);
    return x_i;
}

Eigen::VectorXf LMSolver::solvePi(Eigen::VectorXf &p_i) {
    Eigen::SparseMatrix<float> I(H_.rows(), H_.cols());
    I.setIdentity();
    Eigen::SparseMatrix<float> M = H_ + mu_ * I;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<float>> llt(M);
    p_i = llt.solve(J_);
    //dbg(p_i);
    return p_i;
}

float LMSolver::lineSearchSigma(const Eigen::VectorXf &p_i, float sigma_i) {
    dbg("line search sigma...");
    Eigen::VectorXf V_i = x_ - sigma_i * p_i;
    if (lim_energy_.value(V_i) >= lim_energy_.value(x_)) {
        dbg("energy increase");
        while (lim_energy_.value(V_i) > lim_energy_.value(x_) && sigma_i > 1e-10 && sigma_i <= sigma_max_) {
            sigma_i *= 0.5;
            V_i = x_ - sigma_i * p_i;
        }
    } else {
        dbg("energy decrease");
        while (lim_energy_.value(V_i) < lim_energy_.value(x_) && sigma_i <= 0.5 * sigma_max_ && sigma_i > 1e-10) {
            sigma_i *= 2;
            V_i = x_ - sigma_i * p_i;
        }
    }
    if (sigma_i < 1e-10) sigma_i = 0;
    return sigma_i;
}

float LMSolver::lineSearchMu(float mu_i) {
    dbg("line search mu...");
    if (xry_mesh::isSparseMatrixInvertible(H_)) return 0;
    Eigen::SparseMatrix<float> I(H_.rows(), H_.cols());
    I.setIdentity();
    Eigen::SparseMatrix<float> matrix = H_ + mu_i * I;
    if (xry_mesh::isSparseMatrixInvertible(matrix)) {
        dbg("matrix invertible");
        while (xry_mesh::isSparseMatrixInvertible(matrix) && mu_i > 1e-6 && mu_i <= mu_max_) {
            mu_i *= 0.5;
            matrix = H_ + mu_i * I;
        }
    } else {
        dbg("matrix not invertible");
        while (!xry_mesh::isSparseMatrixInvertible(matrix) && mu_i <= 0.5 * mu_max_) {
        	if (mu_i == 0) mu_i = 1e-5;
            mu_i *= 2;
            matrix = H_ + mu_i * I;
        }
    }
    return mu_i;
}

float LMSolver::computeError() {
    return lim_energy_.value();
}

void LMSolver::compareNumerical(const Eigen::VectorXf &x,
                                const Eigen::VectorXf &x_i) {
    const float val0 = lim_energy_.value(x);
    const float eps = 1e-3;
    Eigen::VectorXf f_i(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        Eigen::VectorXf tmp_x_i = x;
        tmp_x_i[i] += eps;
        const float val_i = lim_energy_.value(tmp_x_i);
        f_i[i] = val_i;
    }
}

void LMSolver::setParameters(float alpha, float beta, float sigma_max, float mu_max, float r, float t, float s_j) {
    setAlpha(alpha);
    setBeta(beta);
    setSigmaMax(sigma_max);
    setSigma(sigma_max);
    setMuMax(mu_max);
    setMu(mu_max);
    setR(r);
    setT(t);
    setSJ(s_j);
}

float LMSolver::getMuMax() const {
    return mu_max_;
}

void LMSolver::setMuMax(float muMax) {
    mu_max_ = muMax;
}

float LMSolver::getSigmaMax() const {
    return sigma_max_;
}

void LMSolver::setSigmaMax(float sigmaMax) {
    sigma_max_ = sigmaMax;
}

float LMSolver::getAlpha() const {
    return alpha_;
}

void LMSolver::setAlpha(float alpha) {
    alpha_ = alpha;
}

float LMSolver::getBeta() const {
    return beta_;
}

void LMSolver::setBeta(float beta) {
    beta_ = beta;
}

float LMSolver::getSJ() const {
    return s_j_;
}

void LMSolver::setSJ(float sJ) {
    s_j_ = sJ;
}

float LMSolver::getMu() const {
    return mu_;
}

void LMSolver::setMu(float mu) {
    mu_ = mu;
}

float LMSolver::getSigma() const {
    return sigma_;
}

void LMSolver::setSigma(float sigma) {
    sigma_ = sigma;
}

float LMSolver::getR() const {
    return r_;
}

void LMSolver::setR(float r) {
    r_ = r;
}

float LMSolver::getT() const {
    return t_;
}

void LMSolver::setT(float t) {
    t_ = t;
}

size_t LMSolver::getMaxIter() const {
    return max_iter_;
}

void LMSolver::setMaxIter(size_t maxIter) {
    max_iter_ = maxIter;
}

bool LMSolver::isEnableMuSearch() const {
    return enable_mu_search_;
}

void LMSolver::setEnableMuSearch(bool enableMuSearch) {
    enable_mu_search_ = enableMuSearch;
}

bool LMSolver::isEnableSigmaSearch() const {
    return enable_sigma_search_;
}

void LMSolver::setEnableSigmaSearch(bool enableSigmaSearch) {
    enable_sigma_search_ = enableSigmaSearch;
}

bool LMSolver::isEnableBarrierFunc() const {
    return enable_barrier_func_;
}

void LMSolver::setEnableBarrierFunc(bool enableBarrierFunc) {
    enable_barrier_func_ = enableBarrierFunc;
}

bool LMSolver::isEnableUpdateAlpha() const {
    return enable_update_alpha_;
}

void LMSolver::setEnableUpdateAlpha(bool enableUpdateAlpha) {
    enable_update_alpha_ = enableUpdateAlpha;
}
