//
// Created by 徐溶延 on 2020/6/15.
//

#include "BarrierEnergy.h"

namespace xry_mesh {
	BarrierEnergy::BarrierEnergy() {}

	BarrierEnergy::BarrierEnergy(const Eigen::VectorXf &x,
								 const Eigen::Matrix2Xf &V,
								 const Eigen::Matrix3Xi &F,
								 const std::vector<float> &areas) : BaseEnergy(x),
																	V_(V),
																	F_(F),
																	areas_(areas) {}

	float BarrierEnergy::value() {
		float res = 0;
		dbg(s_j_);
		float MIN_AREA = std::numeric_limits<float>::infinity();
		for (size_t i = 0; i < F_.cols(); i++) {
			float c_val = c(i);
			MIN_AREA = std::min(c_val, MIN_AREA);
			res += phi(c_val);
		}
		dbg(MIN_AREA);
		return res;
	}

	float BarrierEnergy::value(const Eigen::VectorXf &x) {
		auto tmp_x = x;
		Eigen::Matrix2Xf V2d = Eigen::Map<Eigen::Matrix2Xf>(tmp_x.data(), 2, x.size() / 2);
		float res = 0;
		//根据给出的点计算三角形面积
		for (size_t i = 0; i < F_.cols(); i++) {
			std::vector<Eigen::Vector2f> points;
			for (size_t j = 0; j < 3; j++) {
				Eigen::Vector2f vec = V2d.col(F_(j, i));
				points.emplace_back(vec);
			}
			res += phi(c(points[0], points[1], points[2]));
		}
		return res;
	}

	void BarrierEnergy::init() {
		assert(!areas_.empty());
		assert(F_.cols() == areas_.size());
		computeVertexIdxMatrix();
		computeEpsilon();
		computeW1();
		computeW2();
	}

	Eigen::VectorXf BarrierEnergy::jacobian() const {
		Eigen::VectorXf J(x_.size());
		J.setZero();
		for (size_t f_idx = 0; f_idx < F_.cols(); f_idx++) {
			Eigen::Matrix<int, 6, 1> indices = triangle_vertex_idx_.col(f_idx);
			Eigen::Vector2f A(x_[indices[0]], x_[indices[1]]);
			Eigen::Vector2f B(x_[indices[2]], x_[indices[3]]);
			Eigen::Vector2f C(x_[indices[4]], x_[indices[5]]);
			float term = weight1[f_idx];
			J[indices[0]] += 0.5 * term * (B[1] - C[1]);
			J[indices[1]] += 0.5 * term * (C[0] - B[0]);
			J[indices[2]] += 0.5 * term * (C[1] - A[1]);
			J[indices[3]] += 0.5 * term * (A[0] - C[0]);
			J[indices[4]] += 0.5 * term * (A[1] - B[1]);
			J[indices[5]] += 0.5 * term * (B[0] - A[0]);
		}
		return J;
	}

	std::vector<Eigen::VectorXf> BarrierEnergy::lambdaJacobian() const {
		std::vector<Eigen::VectorXf> vec_J;

		for (size_t f_idx = 0; f_idx < F_.cols(); f_idx++) {
			Eigen::VectorXf J(x_.size());
			J.setZero();
			Eigen::Matrix<int, 6, 1> indices = triangle_vertex_idx_.col(f_idx);
			Eigen::Vector2f A(x_[indices[0]], x_[indices[1]]);
			Eigen::Vector2f B(x_[indices[2]], x_[indices[3]]);
			Eigen::Vector2f C(x_[indices[4]], x_[indices[5]]);
			J[indices[0]] += 0.5 * (B[1] - C[1]);
			J[indices[1]] += 0.5 * (C[0] - B[0]);
			J[indices[2]] += 0.5 * (C[1] - A[1]);
			J[indices[3]] += 0.5 * (A[0] - C[0]);
			J[indices[4]] += 0.5 * (A[1] - B[1]);
			J[indices[5]] += 0.5 * (B[0] - A[0]);
			vec_J.emplace_back(J);
		}
		return vec_J;
	}

	Eigen::SparseMatrix<float> BarrierEnergy::hessian() const {
		Eigen::SparseMatrix<float> H(x_.size(), x_.size());
		H.setZero();
		std::vector<Eigen::VectorXf> vec_J = lambdaJacobian();
		// part1
		for (size_t f_idx = 0; f_idx < F_.cols(); f_idx++) {
			auto J = vec_J[f_idx];
			Eigen::Matrix<int, 6, 1> indices = triangle_vertex_idx_.col(f_idx);
			float term = weight2[f_idx];
			H += term * (J * J.transpose()).sparseView();
		}
		const int non_flip_hessian_idx[6][2] = {{0, 3},
												{0, 5},
												{1, 2},
												{1, 4},
												{2, 5},
												{3, 4}};

		// part2
		for (size_t f_idx = 0; f_idx < F_.cols(); f_idx++) {
			Eigen::Matrix<int, 6, 1> indices = triangle_vertex_idx_.col(f_idx);
			float term = weight1[f_idx];
			H.coeffRef(indices[non_flip_hessian_idx[0][0]], indices[non_flip_hessian_idx[0][1]]) += term;
			H.coeffRef(indices[non_flip_hessian_idx[0][1]], indices[non_flip_hessian_idx[0][0]]) += term;
			H.coeffRef(indices[non_flip_hessian_idx[1][0]], indices[non_flip_hessian_idx[1][1]]) += -term;
			H.coeffRef(indices[non_flip_hessian_idx[1][1]], indices[non_flip_hessian_idx[1][0]]) += -term;
			H.coeffRef(indices[non_flip_hessian_idx[2][0]], indices[non_flip_hessian_idx[2][1]]) += -term;
			H.coeffRef(indices[non_flip_hessian_idx[2][1]], indices[non_flip_hessian_idx[2][0]]) += -term;
			H.coeffRef(indices[non_flip_hessian_idx[3][0]], indices[non_flip_hessian_idx[3][1]]) += term;
			H.coeffRef(indices[non_flip_hessian_idx[3][1]], indices[non_flip_hessian_idx[3][0]]) += term;
			H.coeffRef(indices[non_flip_hessian_idx[4][0]], indices[non_flip_hessian_idx[4][1]]) += term;
			H.coeffRef(indices[non_flip_hessian_idx[4][1]], indices[non_flip_hessian_idx[4][0]]) += term;
			H.coeffRef(indices[non_flip_hessian_idx[5][0]], indices[non_flip_hessian_idx[5][1]]) += -term;
			H.coeffRef(indices[non_flip_hessian_idx[5][1]], indices[non_flip_hessian_idx[5][0]]) += -term;
		}
		return H;
	}

	void BarrierEnergy::update(const Eigen::VectorXf &x) {
		BaseEnergy::update(x);
		V_ = Eigen::Map<Eigen::Matrix2Xf>(x_.data(), 2, x_.size() / 2);
		updateAreas();
		computeW1();
		computeW2();
	}

	void BarrierEnergy::computeW1() {
		weight1.clear();
		weight1.resize(F_.cols());
		for (size_t i = 0; i < F_.cols(); i++) {
			weight1[i] = -gPrime(c(i)) / std::pow(g(c(i)), 2);
		}
	}

	void BarrierEnergy::computeW2() {
		weight2.clear();
		weight2.resize(F_.cols());
		for (size_t i = 0; i < F_.cols(); i++) {
			weight2[i] = (2 * std::pow(gPrime(c(i)), 2) - gPrime2(c(i)) * g(c(i)))
						 / std::pow(g(c(i)), 3);
		}
	}

	void BarrierEnergy::computeEpsilon() {
		float MIN = std::numeric_limits<float>::infinity();
		for (size_t i = 0; i < F_.cols(); i++) {
			auto area = areas_[i];
			MIN = std::min(MIN, area);
		}
		epsilon_ = 1e-5 * MIN;
		if (enable_dbg) {
			dbg(epsilon_);
		}
	}

	void BarrierEnergy::computeVertexIdxMatrix() {
		triangle_vertex_idx_.resize(6, F_.cols());
		for (size_t i = 0; i < F_.cols(); i++) {
			size_t row = 0;
			for (size_t j = 0; j < F_.rows(); j++) {
				for (size_t k = 0; k < 2; k++, row++) {
					triangle_vertex_idx_(row, i) = 2 * F_(j, i) + k;
				}
			}
		}
	}

	void BarrierEnergy::updateAreas() {
		for (size_t i = 0; i < F_.cols(); i++) {
			std::vector<Eigen::Vector2f> points;
			for (size_t j = 0; j < 3; j++) {
				Eigen::Vector2f vec = V_.col(F_(j, i));
				points.emplace_back(vec);
			}
			areas_[i] = xry_mesh::computeArea2D<float>(points[0], points[1], points[2]);
		}
	}

	float BarrierEnergy::getSJ() const {
		return s_j_;
	}

	void BarrierEnergy::setSJ(float sJ) {
		s_j_ = sJ;
	}

	const std::vector<float> &BarrierEnergy::getWeight1() const {
		return weight1;
	}

	const std::vector<float> &BarrierEnergy::getWeight2() const {
		return weight2;
	}

	const std::vector<float> &BarrierEnergy::getAreas() const {
		return areas_;
	}

	const Eigen::Matrix<int, 6, Eigen::Dynamic> &BarrierEnergy::getTriangleVertexIdx() const {
		return triangle_vertex_idx_;
	}


} // namespace xry_mesh