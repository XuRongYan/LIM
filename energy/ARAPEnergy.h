//
// Created by 徐溶延 on 2020/6/14.
//

#ifndef LIM_ARAPENERGY_H
#define LIM_ARAPENERGY_H

#include "NormEnergy.h"

namespace xry_mesh {
    class ARAPEnergy : public NormEnergy {
    public:
        ARAPEnergy();

        ARAPEnergy(const Eigen::VectorXf &x,
                   const Eigen::Matrix2Xf &V,
                   const Eigen::Matrix3Xi &F,
                   const std::vector<float> &areas);

        float value() override;

        float value(const Eigen::VectorXf &x) override;

        void init() override;

        void update(const Eigen::VectorXf &x) override;

        const Eigen::Matrix2Xf &getV() const;

        const Eigen::Matrix3Xi &getF() const;

        const std::vector<float> &getIdealAreas() const;

        const std::vector<Eigen::Matrix2f> &getVecR() const;

        const std::vector<Eigen::Matrix<float, 2, 3>> &getIdealElems() const;

        const std::vector<Eigen::Matrix2f> &getInvDeltas() const;

    protected:
        void computeA() override;

        void computeB() override;

    private:
        Eigen::Matrix2Xf V_;
        Eigen::Matrix3Xi F_;
        std::vector<float> ideal_areas_;
        std::vector<Eigen::Matrix2f> vec_R_;
        std::vector<Eigen::Matrix<float, 2, 3>> ideal_elems_;
        std::vector<Eigen::Matrix2f> inv_deltas_;

        /**
         * 更新local旋转矩阵
         */
        void localPhase();

        /**
         * 计算ideal elements
         */
        void computeIdeals();

        void computeInvDetails();

        /**
         * locally compute rotation matrix of a triangle
         * @param J PQ^{-1}
         * @param R rotation matrix
         * @return 1 if success, 0 if failed
         */
        void optimizeRotation(const Eigen::Matrix2f &J,
                              Eigen::Matrix2f &R);

        /**
         * flatten triangle to plane to get ideal elements
         * @param p1
         * @param p2
         * @param p3
         */
        Eigen::Matrix<float, 2, 3> flattenTriangle(const Eigen::Vector2f &p1,
                                                    const Eigen::Vector2f &p2,
                                                    const Eigen::Vector2f &p3);

        /**
         * 原网格的三角形直接作为ideal
         * @return
         */
        Eigen::Matrix<float, 2, 3> useOriginalTriangles(const Eigen::Vector2f &p1, const Eigen::Vector2f &p2, const Eigen::Vector2f &p3);

        /**
         * 计算梯度算子矩阵
         * @return
         */
        Eigen::SparseMatrix<float> computeGradientOperator();

        /**
         * 将vector形式旋转矩阵转换为求解向量
         * @return
         */
        Eigen::VectorXf recomposeR();

        /**
         * 重组单个旋转矩阵
         * @param R
         * @param i
         */
        void recomposeRi(Eigen::VectorXf &R, size_t i);
    };


} // namespace xry_mesh


#endif //LIM_ARAPENERGY_H
