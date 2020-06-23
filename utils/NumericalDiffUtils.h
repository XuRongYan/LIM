//
// Created by 徐溶延 on 2020/6/21.
//

#ifndef LIM_NUMERICALDIFFUTILS_H
#define LIM_NUMERICALDIFFUTILS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace xry_mesh {

    /**
     * 前向差分公式求解数值微分
     * @param fx_1 f(x + h)
     * @param fx_0 f(x)
     * @param x1 x + h
     * @param x0 x
     * @return
     */
    Eigen::VectorXf numericalGrad(const Eigen::VectorXf &fx_1,
                                  float fx_0,
                                  float eps);

} // namespace xry_mesh


#endif //LIM_NUMERICALDIFFUTILS_H
