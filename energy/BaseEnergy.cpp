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
}


