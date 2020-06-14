//
// Created by 徐溶延 on 2020/6/1.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_EIGENUTILS_H
#define LOCALLYINJECTIVEMAPPINGS_EIGENUTILS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <dbg.h>

namespace xry_mesh {

    /**
     * 将vt矩阵转换为v矩阵
     * @tparam Scalar
     * @param vt = {{x1, y1, z1};{x2, y2, z2}; ...  ;{xn, yn, zn}}
     * @return v = {x1, y1, z1, x2, y2, z2, ..., xn, yn, zn}
     */
    template<class Scalar>
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    vt2v(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &vt) {
        const int dim = vt.cols();
        //dbg(vt.size());
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> v(dim * vt.rows());
        for (size_t i = 0; i < vt.rows(); i++) {
            for (size_t j = 0; j < dim; j++) {
                v(dim * i + j, 0) = vt(i, j);
            }
        }
        return v;
    }

    /**
     * 将v矩阵转换为vt矩阵
     * @tparam Scalar
     * @param v = {x1, y1, z1, x2, y2, z2, ..., xn, yn, zn}
     * @return vt = {{x1, y1, z1};{x2, y2, z2}; ...  ;{xn, yn, zn}}
     */
    template<class Scalar>
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    v2vt(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &v, int dim) {
        assert(v.rows() % dim == 0);
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> vt(v.rows() / dim, dim);
        for (size_t i = 0; i < v.rows(); i += dim) {
            for (size_t j = 0; j < dim; j++) {
                vt(i / dim, j) = v(i + j, 0);
            }
        }

        return vt;
    }

    template<class Scalar>
    std::vector<Eigen::Triplet<Scalar>> sparse2triplet(const Eigen::SparseMatrix<Scalar> &A) {
        std::vector<Eigen::Triplet<Scalar>> triplets;
        for (size_t i = 0; i < A.outerSize(); i++) {
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(A, i); it; it++) {
                triplets.emplace_back(it.row(), it.col(), it.value());
            }
        }
        return triplets;
    }



    template<typename Scalar>
    bool isSparseMatrixInvertible(Eigen::SparseMatrix<Scalar> &matrix) {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<Scalar>> llt(matrix);
        return llt.info() == Eigen::Success;
    }

    template<typename Scalar>
    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> v2dTo3d(const Eigen::Matrix<Scalar, 2, Eigen::Dynamic> &V2d) {
        Eigen::Matrix<Scalar, 3, Eigen::Dynamic> V(3, V2d.cols());
        for (size_t i = 0; i < V.cols(); i++) {
            V.col(i) << V2d(0, i), V2d(1, i), 0;
        }
        return V;
    }

} // namespace xry_mesh

#endif //LOCALLYINJECTIVEMAPPINGS_EIGENUTILS_H
