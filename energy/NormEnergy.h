//
// Created by 徐溶延 on 2020/6/14.
//

#ifndef LIM_NORMENERGY_H
#define LIM_NORMENERGY_H

#include "BaseEnergy.h"

namespace xry_mesh {
    /**
     * 二范数能量函数
     */
    class NormEnergy : public BaseEnergy {
    public:
        NormEnergy();

        NormEnergy(const Eigen::VectorXf &x);

        float value() override;

        virtual float value(const Eigen::VectorXf &x) override;

        void init() override;

        Eigen::VectorXf jacobian() const override;

        Eigen::SparseMatrix<float> hessian() const override;

        virtual void update(const Eigen::VectorXf &x) override;

        const Eigen::SparseMatrix<float> &getA() const;

        const Eigen::VectorXf &getB() const;

    protected:
        virtual void computeA() = 0;

        virtual void computeB() = 0;

    protected:
        Eigen::SparseMatrix<float> A_;
        Eigen::VectorXf b_;
    };
}


#endif //LIM_NORMENERGY_H
