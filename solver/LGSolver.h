//
// Created by 徐溶延 on 2020/6/15.
//

#ifndef LIM_LGSOLVER_H
#define LIM_LGSOLVER_H

#include <SurfaceMesh/SurfaceMesh.h>

#include "../energy/ARAPEnergy.h"
#include "../energy/PosEnergy.h"

class LGSolver {
public:
    LGSolver();

    LGSolver(const Surface_Mesh::SurfaceMesh &mesh,
             const std::vector<std::pair<int, Eigen::VectorXf>> &posConstrains,
             size_t maxIter);

    void init();

    Eigen::VectorXf solve();

    void subSolve();

private:
    float alpha = 1e8;
    size_t max_iter = 100;
    Surface_Mesh::SurfaceMesh mesh_;
    Eigen::Matrix2Xf V_;
    Eigen::Matrix3Xi F_;
    std::vector<std::pair<int, Eigen::VectorXf>> pos_constrains_;
    Eigen::SparseMatrix<float> L_;
    Eigen::VectorXf b_;
    Eigen::VectorXf x_;
    xry_mesh::ARAPEnergy arap_energy_;
    xry_mesh::PosEnergy pos_energy_;
    std::vector<float> areas;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<float>> llt_;
    Eigen::SparseMatrix<float> arap_A, pos_A;

    float computeError();

};


#endif //LIM_LGSOLVER_H
