//
// Created by 徐溶延 on 2020/6/17.
//

#include "LimEnergy.h"

namespace xry_mesh {

	LimEnergy::LimEnergy() {}

	LimEnergy::LimEnergy(const Surface_Mesh::SurfaceMesh &mesh,
						 const std::vector<std::pair<int, Eigen::VectorXf>> &posConstrains) : mesh_(mesh),
																							  pos_constrains_(
																									  posConstrains) {}

	float LimEnergy::value() {
		float val = 0;
		float arap_err = arap_energy_.value();
		dbg(arap_err);
		val += arap_err;
		float pos_err = alpha_ * pos_energy_.value();
		dbg(pos_err);
		val += pos_err;
		if (enable_barrier_func_) {
			float barrier_err = beta_ * barrier_energy_.value();
			dbg(barrier_err);
			val += barrier_err;
		}
		return val;
	}

	float LimEnergy::value(const Eigen::VectorXf &x) {
		float val = 0;
		float arap_err = arap_energy_.value(x);
		//dbg(arap_err);
		val += arap_err;
		float pos_err = alpha_ * pos_energy_.value(x);
		//dbg(pos_err);
		val += pos_err;
		if (enable_barrier_func_) {
			float barrier_err = beta_ * barrier_energy_.value(x);
			//dbg(barrier_err);
			val += barrier_err;
		}
		return val;
	}

	void LimEnergy::init() {
		assert(!mesh_.empty());
		Eigen::Matrix3Xf V3d = xry_mesh::getGeometryMatrix<float>(mesh_);
		V_ = V3d.block(0, 0, 2, V3d.cols());
		F_ = xry_mesh::getTopoMatrix(mesh_);
		x_ = xry_mesh::vt2v<float>(V_.transpose());
		xry_mesh::computeAreasDet<float>(mesh_);
		areas_ = xry_mesh::faceProperty2StdVector<float>(mesh_, "f:areas");

		arap_energy_ = xry_mesh::ARAPEnergy(x_, V_, F_, areas_);
		arap_energy_.init();

		pos_energy_ = xry_mesh::PosEnergy(x_, 2, pos_constrains_, mu_);
		pos_energy_.init();

		if (enable_barrier_func_) {
			barrier_energy_ = xry_mesh::BarrierEnergy(x_, V_, F_, areas_);
			barrier_energy_.setSJ(s_j_);
			barrier_energy_.init();
		}
		if (enable_update_alpha_) {
			computeAlpha();
		}

	}

	Eigen::VectorXf LimEnergy::jacobian() const {
		Eigen::VectorXf res(x_.size());
		res.setZero();
		res += arap_energy_.jacobian();
		res += alpha_ * pos_energy_.jacobian();
		if (enable_barrier_func_) {
			res += beta_ * barrier_energy_.jacobian();
		}
		return res;
	}

	Eigen::VectorXf LimEnergy::numericalJacobian(float esp) {
		return BaseEnergy::numericalJacobian(esp);
//		Eigen::VectorXf res(x_.size());
//		res.setZero();
//		res += arap_energy_.numericalJacobian(esp);
//		res += alpha_ * pos_energy_.numericalJacobian(esp);
//		if (enable_barrier_func_) {
//			res += beta_ * barrier_energy_.numericalJacobian(esp);
//		}
//		return res;
	}

	Eigen::SparseMatrix<float> LimEnergy::hessian() const {
		Eigen::SparseMatrix<float> res(x_.size(), x_.size());
		res.setZero();
		res += arap_energy_.hessian();
		res += alpha_ * pos_energy_.hessian();
		if (enable_barrier_func_) {
			res += beta_ * barrier_energy_.hessian();
		}
		return res;
	}

	void LimEnergy::update(const Eigen::VectorXf &x) {
		BaseEnergy::update(x);
		V_ = Eigen::Map<Eigen::Matrix2Xf>(x_.data(), 2, x_.size() / 2);
		arap_energy_.update(x);
		pos_energy_.setMu(mu_);
		pos_energy_.update(x);
		if (enable_barrier_func_) {
			barrier_energy_.update(x);
		}
		if (enable_update_alpha_) {
			computeAlpha();
		}
	}

	bool LimEnergy::isEnableBarrierFunc() const {
		return enable_barrier_func_;
	}

	void LimEnergy::setEnableBarrierFunc(bool enableBarrierFunc) {
		enable_barrier_func_ = enableBarrierFunc;
	}

	bool LimEnergy::isEnableUpdateAlpha() const {
		return enable_update_alpha_;
	}

	void LimEnergy::setEnableUpdateAlpha(bool enableUpdateAlpha) {
		enable_update_alpha_ = enableUpdateAlpha;
	}

	float LimEnergy::getAlpha() const {
		return alpha_;
	}

	void LimEnergy::setAlpha(float alpha) {
		alpha_ = alpha;
	}

	float LimEnergy::getBeta() const {
		return beta_;
	}

	void LimEnergy::setBeta(float beta) {
		beta_ = beta;
	}

	float LimEnergy::getSigmaMax() const {
		return sigma_max_;
	}

	void LimEnergy::setSigmaMax(float sigmaMax) {
		sigma_max_ = sigmaMax;
	}

	float LimEnergy::getSigma() const {
		return sigma_;
	}

	void LimEnergy::setSigma(float sigma) {
		sigma_ = sigma;
	}

	float LimEnergy::getMuMax() const {
		return mu_max_;
	}

	void LimEnergy::setMuMax(float muMax) {
		mu_max_ = muMax;
	}

	float LimEnergy::getMu() const {
		return mu_;
	}

	void LimEnergy::setMu(float mu) {
		mu_ = mu;
		pos_energy_.setMu(mu);
	}

	float LimEnergy::getGamma() const {
		return gamma_;
	}

	void LimEnergy::setGamma(float gamma) {
		gamma_ = gamma;
	}

	float LimEnergy::getR() const {
		return r_;
	}

	void LimEnergy::setR(float r) {
		r_ = r;

	}

	float LimEnergy::getT() const {
		return t_;
	}

	void LimEnergy::setT(float t) {
		t_ = t;
	}

	const std::vector<float> &LimEnergy::getAreas() const {
		return areas_;
	}

	void LimEnergy::setAreas(const std::vector<float> &areas) {
		LimEnergy::areas_ = areas;
	}

	float LimEnergy::getSJ() const {
		return s_j_;
	}

	void LimEnergy::setSJ(float sJ) {
		s_j_ = sJ;
	}
} // namespace xry_mesh