//
// Created by 徐溶延 on 2020/6/17.
//

#ifndef LIM_LIMENERGY_H
#define LIM_LIMENERGY_H

#include "BaseEnergy.h"
#include "ARAPEnergy.h"
#include "PosEnergy.h"
#include "BarrierEnergy.h"


namespace xry_mesh {
	class LimEnergy : public BaseEnergy {
	public:
		LimEnergy();

		LimEnergy(const Surface_Mesh::SurfaceMesh &mesh,
				  const std::vector<std::pair<int, Eigen::VectorXf>> &posConstrains);

		float value() override;

		float value(const Eigen::VectorXf &x) override;

		void init() override;

		Eigen::VectorXf jacobian() const override;

		Eigen::VectorXf numericalJacobian(float esp) override;

		Eigen::SparseMatrix<float> hessian() const override;

		void update(const Eigen::VectorXf &x) override;

		bool isEnableBarrierFunc() const;

		void setEnableBarrierFunc(bool enableBarrierFunc);

		bool isEnableUpdateAlpha() const;

		void setEnableUpdateAlpha(bool enableUpdateAlpha);

		float getAlpha() const;

		void setAlpha(float alpha);

		float getBeta() const;

		void setBeta(float beta);

		float getSigmaMax() const;

		void setSigmaMax(float sigmaMax);

		float getSigma() const;

		void setSigma(float sigma);

		float getMuMax() const;

		void setMuMax(float muMax);

		float getMu() const;

		void setMu(float mu);

		float getGamma() const;

		void setGamma(float gamma);

		float getR() const;

		void setR(float r);

		float getT() const;

		void setT(float t);

		const std::vector<float> &getAreas() const;

		void setAreas(const std::vector<float> &areas);

		float getSJ() const;

		void setSJ(float sJ);

	private:
		bool enable_barrier_func_ = true;
		bool enable_update_alpha_ = true;
		float alpha_{1e8};
		float beta_{0.1};
		float sigma_max_{1.0};
		float sigma_{1.0};
		float mu_max_{1.0};
		float mu_{1.0};
		float gamma_{1.0};
		float r_{1000};
		float t_{1e16};
		float s_j_{1.0};
		ARAPEnergy arap_energy_;
		PosEnergy pos_energy_;
		BarrierEnergy barrier_energy_;
		Surface_Mesh::SurfaceMesh mesh_;
		Eigen::Matrix2Xf V_;
		Eigen::Matrix3Xi F_;
		std::vector<float> areas_;
		std::vector<std::pair<int, Eigen::VectorXf>> pos_constrains_;

		inline void computeAlpha() {
			computeGamma();
			alpha_ = std::min(std::max(alpha_, gamma_), t_);
			dbg(alpha_);
		}

		inline void computeGamma() {
			gamma_ = r_ * (arap_energy_.value() + beta_ * barrier_energy_.value()) / pos_energy_.value();
		}
	};
} // namespace xry_mesh


#endif //LIM_LIMENERGY_H
