//
// Created by 徐溶延 on 2020/6/15.
//

#ifndef LIM_BARRIERENERGY_H
#define LIM_BARRIERENERGY_H

#include <SurfaceMesh/SurfaceMesh.h>

#include "../energy/BaseEnergy.h"
#include "../utils/SurfaceMeshUtils.h"
#include "../utils/EigenUtils.h"

namespace xry_mesh {
    class BarrierEnergy : public BaseEnergy {
    public:
        BarrierEnergy();

        BarrierEnergy(const Eigen::VectorXf &x, const Eigen::Matrix2Xf &v, const Eigen::Matrix3Xi &f);

        float value() override;

        float value(const Eigen::VectorXf &x) override;

        void init() override;

        Eigen::VectorXf jacobian() const override;

        Eigen::SparseMatrix<float> hessian() const override;

        void update(const Eigen::VectorXf &x) override;

    private:
        float s_j_ = 1.0;
        float epsilon_;
        Eigen::Matrix2Xf V_;
        Eigen::Matrix3Xi F_;
        Eigen::Matrix<int, 6, Eigen::Dynamic> triangle_vertex_idx_;
        std::vector<float> weight1, weight2;
        std::vector<float> areas;

        /**
         * 计算每个三角形三个顶点对应求解矩阵的位置
         */
        void computeVertexIdxMatrix();

        /**
         * 计算epsilon，具体见论文
         */
        void computeEpsilon();

        /**
         * 计算其中一个权值，见论文
         */
        void computeW1();

        /**
         * 计算其中一个权值，见论文
         */
        void computeW2();

        /**
        * 分段函数phi，具体见论文
        * @param x
        * @return
        */
        inline float phi(float x) {
            if (x <= 0) return std::numeric_limits<float>::infinity();
            else if (x > 0 && x < s_j_) return 1.0 / g(x) - 1;
            else return 0;
        }

        /**
         * 多项式函数g，具体见论文
         * @return
         */
        inline float g(float x) const {
            return 1.0 / std::pow(s_j_, 3) * std::pow(x, 3)
                   + 3.0 / std::pow(s_j_, 2) * std::pow(x, 2)
                   + 3.0 / s_j_ * x;
        }

        /**
        * g(x)一阶导
        * @return
        */
        inline float gPrime(float x) const {
            return 3.0 / std::pow(s_j_, 3) * std::pow(x, 2)
                   + 6.0 / std::pow(s_j_, 2) * x
                   + 3.0 / s_j_;
        }

        /**
        * g(x)二阶导
        * @param x
        * @return
        */
        inline float gPrime2(float x) const {
            return 6.0 / std::pow(s_j_, 3) * x
                   + 6.0 / std::pow(s_j_, 2);
        }

        /**
        * 面积约束函数c
        * @param p1
        * @param p2
        * @param p3
        * @return
        */
        inline float c(const Surface_Mesh::Point &p1,
                       const Surface_Mesh::Point &p2,
                       const Surface_Mesh::Point &p3) const {
            return xry_mesh::computeArea<float>(p1, p2, p3) - epsilon_;
        };

    };
} // namespace xry_mesh


#endif //LIM_BARRIERENERGY_H
