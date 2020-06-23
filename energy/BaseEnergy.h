//
// Created by 徐溶延 on 2020/6/14.
//

#ifndef LIM_BASEENERGY_H
#define LIM_BASEENERGY_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "dbg.h"

namespace xry_mesh {
    /**
     * 能量函数基类
     */
    class BaseEnergy {
    public:
        BaseEnergy();

        explicit BaseEnergy(const Eigen::VectorXf &x);

        virtual float value() = 0;

        virtual float value(const Eigen::VectorXf &x) = 0;

        virtual void init() = 0;

        virtual Eigen::VectorXf jacobian() const = 0;

        virtual Eigen::SparseMatrix<float> hessian() const = 0;

		virtual /**
         * 数值方法计算梯度矩阵
         * @return
         */
        Eigen::VectorXf numericalJacobian(float esp = 1e-4);

        virtual void update(const Eigen::VectorXf &x);

        const Eigen::VectorXf &getX() const;

        void setX(const Eigen::VectorXf &x);

    protected:
        const static bool enable_dbg = false;
        Eigen::VectorXf x_;                 //求解变量（固定为列向量）
    };
} // namespace xry_mesh


#endif //LIM_BASEENERGY_H
