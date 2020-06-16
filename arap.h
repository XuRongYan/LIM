#ifndef SINGLE_PATCH_PARAMETERIZATION_H
#define SINGLE_PATCH_PARAMETERIZATION_H

#include <vector>
#include <Eigen/Dense>

namespace jy_mesh {
	int arap_deformation(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F,
		const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
		Eigen::Matrix2Xd &Vout, size_t max_iter = 1000);
}

#endif
