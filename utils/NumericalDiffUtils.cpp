//
// Created by 徐溶延 on 2020/6/21.
//

#include <dbg.h>

#include "NumericalDiffUtils.h"


namespace xry_mesh {
    Eigen::VectorXf numericalGrad(const Eigen::VectorXf &fx_1,
                                  float fx_0,
                                  float eps) {
        assert(eps != 0);
        Eigen::VectorXf grad(fx_1.size());
        for (size_t i = 0; i < fx_1.size(); i++) {
            grad[i] = (fx_1[i] - fx_0) / eps;
        }
        //dbg(grad);
        return grad;
    }
} // namespace xry_mesh