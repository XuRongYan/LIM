//
// Created by 徐溶延 on 2020/6/14.
//

#include "BaseEnergy.h"


namespace xry_mesh {
    BaseEnergy::BaseEnergy() {}

    BaseEnergy::BaseEnergy(const Eigen::VectorXf &x) : x_(x) {}

    void BaseEnergy::update(const Eigen::VectorXf &x) {
        x_ = x;
    }

    const Eigen::VectorXf &BaseEnergy::getX() const {
        return x_;
    }

    void BaseEnergy::setX(const Eigen::VectorXf &x) {
        update(x);
    }

	Eigen::VectorXf BaseEnergy::numericalJacobian(float esp) {
    	Eigen::VectorXf J(x_.size());
    	J.setZero();
		const float val = value();
		for (size_t f_idx = 0; f_idx < x_.size(); f_idx++) {
			Eigen::VectorXf x_i = x_;
			x_i[f_idx] += esp;
			const float val_i = value(x_i);
			J[f_idx] = (val_i - val) / esp;
		}
		return J;
	}
}


