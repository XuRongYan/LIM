//
// Created by 徐溶延 on 2020/6/16.
//

#ifndef LIM_LMSOLVER_H
#define LIM_LMSOLVER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "dbg.h"
#include "../utils/EigenUtils.h"
#include "../energy/LimEnergy.h"


class LMSolver {
public:
    LMSolver();

    LMSolver(const Surface_Mesh::SurfaceMesh &mesh,
             const std::vector<std::pair<int, Eigen::VectorXf>> &posConstrains,
             size_t maxIter);

    void setParameters(float alpha, float beta, float sigma_max, float mu_max, float r, float t, float s_j);

    void init();

    Eigen::VectorXf solve();

    float getMuMax() const;

    void setMuMax(float muMax);

    float getSigmaMax() const;

    void setSigmaMax(float sigmaMax);

    float getAlpha() const;

    void setAlpha(float alpha);

    float getBeta() const;

    void setBeta(float beta);

    float getSJ() const;

    void setSJ(float sJ);

    float getMu() const;

    void setMu(float mu);

    float getSigma() const;

    void setSigma(float sigma);

	float getMuMin() const;

	void setMuMin(float muMin);

	float getMuCurr() const;

	void setMuCurr(float muCurr);

	float getSigmaMin() const;

	void setSigmaMin(float sigmaMin);

	float getSigmaCurr() const;

	void setSigmaCurr(float sigmaCurr);

	float getR() const;

    void setR(float r);

    float getT() const;

    void setT(float t);

    size_t getMaxIter() const;

    void setMaxIter(size_t maxIter);

    bool isEnableMuSearch() const;

    void setEnableMuSearch(bool enableMuSearch);

    bool isEnableSigmaSearch() const;

    void setEnableSigmaSearch(bool enableSigmaSearch);

    bool isEnableBarrierFunc() const;

    void setEnableBarrierFunc(bool enableBarrierFunc);

    bool isEnableUpdateAlpha() const;

    void setEnableUpdateAlpha(bool enableUpdateAlpha);

private:
    float mu_max_{1.0};
    float mu_min_{1e-8};
    float mu_curr_{1.0};
    float sigma_max_{1.0};
    float sigma_min_{1e-8};
	float sigma_curr_{1.0};
    float alpha_{1e16};
    float beta_{0.1};
    float s_j_{1.0};
    float mu_{1.0};
    float sigma_{1.0};
    float r_{1000};
    float t_{1e16};
    size_t max_iter_{100};
    bool enable_mu_search_{true};
    bool enable_sigma_search_{true};
    bool enable_barrier_func_{false};
    bool enable_update_alpha_{true};

    Surface_Mesh::SurfaceMesh mesh_;
    std::vector<std::pair<int, Eigen::VectorXf>> pos_constrains_;
    Eigen::VectorXf J_;
    Eigen::SparseMatrix<float> H_;
    xry_mesh::LimEnergy lim_energy_;
    Eigen::VectorXf x_;


    Eigen::VectorXf subSolve();

    Eigen::VectorXf solvePi(Eigen::VectorXf &p_i);

    float lineSearchSigma(const Eigen::VectorXf &p_i, float sigma_i);

    float lineSearchMu(float mu_i);

    float computeError();

    void compareNumerical(const Eigen::VectorXf &x,
                          const Eigen::VectorXf &x_i);


};


#endif //LIM_LMSOLVER_H
