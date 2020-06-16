//
// Created by 徐溶延 on 2020/6/7.
//

#ifndef LOCALLYINJECTIVEMAPPINGS_FILEIOUTILS_H
#define LOCALLYINJECTIVEMAPPINGS_FILEIOUTILS_H
#include <vector>
#include <fstream>

#include <Eigen/Dense>
#include <dbg.h>

namespace xry_mesh {
    /**
     * 读取位置约束文件
     * @tparam Scalar
     * @param filename
     * @return
     */
    template<typename Scalar>
    std::vector<std::pair<int, Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>> readPosFile(const std::string &filename, int dim) {
        std::vector<std::pair<int, Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>> res;
        std::ifstream ifs(filename);
        if (ifs.fail()) {
            dbg("can not open", filename);
            return res;
        }
        int id;
        Scalar x, y, z;
        while (ifs >> id >> x >> y >> z) {
            Eigen::Matrix<Scalar, Eigen::Dynamic, 1> point(dim, 1);
            point(0, 0) = x;
            point(1, 0) = y;
            if (dim == 3) point(2, 0) = z;
            res.emplace_back(id, point);
        }
        ifs.close();
        return res;
    }

    template<typename Scalar>
    std::vector<std::pair<size_t , Eigen::Matrix<Scalar, 2, 1>>> readPosFile2D(const std::string &filename, int dim = 2) {
        assert(dim == 2);
        std::vector<std::pair<size_t , Eigen::Matrix<Scalar, 2, 1>>> res;
        std::ifstream ifs(filename);
        if (ifs.fail()) {
            dbg("can not open", filename);
            return res;
        }
        int id;
        Scalar x, y, z;
        while (ifs >> id >> x >> y >> z) {
            Eigen::Matrix<Scalar, Eigen::Dynamic, 1> point(dim, 1);
            point(0, 0) = x;
            point(1, 0) = y;
            //if (dim == 3) point(2, 0) = z;
            res.emplace_back(id, point);
        }
        ifs.close();
        return res;
    }
} // namespace xry_mesh


#endif //LOCALLYINJECTIVEMAPPINGS_FILEIOUTILS_H
