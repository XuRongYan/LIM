//
// Created by 徐溶延 on 2020/6/14.
//

#ifndef LIM_POSENERGY_H
#define LIM_POSENERGY_H

#include "NormEnergy.h"

namespace xry_mesh {
    class PosEnergy : public NormEnergy{
    public:
        PosEnergy();

        PosEnergy(const Eigen::VectorXf &x,
                  size_t dim,
                  const std::vector<std::pair<int, Eigen::VectorXf>> &posConstrains,
                  float mu);

        void init() override;

        void update(const Eigen::VectorXf &x) override;

        const std::vector<std::pair<int, Eigen::VectorXf>> &getPosConstrains() const;

        void setPosConstrains(const std::vector<std::pair<int, Eigen::VectorXf>> &posConstrains);

        size_t getDim() const;

        void setDim(size_t dim);

        bool isEnableTi() const;

        void setEnableTi(bool enableTi);

        float getMu() const;

        void setMu(float mu);



    protected:
        void computeA() override;

        void computeB() override;

    public:
        float value(const Eigen::VectorXf &x) override;

        float value() override;


    private:
        std::vector<std::pair<int, Eigen::VectorXf>> pos_constrains_;
        size_t dim = 2;
        bool enable_ti = false;
        float mu_ = 0;

        /**
         * 计算常规位置约束b
         * @return
         */
        Eigen::VectorXf computeDi();

        /**
         * 计算改进的位置约束
         * @return
         */
        Eigen::VectorXf computeTi(const Eigen::VectorXf &d_i);

        Eigen::VectorXf computeTi(const Eigen::VectorXf &d_i, const Eigen::VectorXf &x);



    };
} // namespace xry_mesh


#endif //LIM_POSENERGY_H
