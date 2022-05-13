//
// Created by lpaulin on 20/05/19.
//

#include <functional>
#include "multivariateGaussian.h"


#include "../Math/myMath.h"


Gaussian::Gaussian(const VecXDynamic& m, const MatXDynamic sigma, bool normalize):mu(m), dim(mu.dim()) {

    invSigma = sigma;


    if (normalize) {
        VecXDynamic eigenV(dim);
        double lambda = invSigma.highestEigenPair(eigenV);
//std::cout << std::endl << lambda << std::endl;
        double alpha = std::sqrt(1. / lambda);

        normFactor = 1. / gradient(alpha * eigenV + mu).norm();
    } else {
        normFactor = 1.;// 1. / std::sqrt(std::pow(2 * M_PI, m.dim()) * sigma.det());
    }

}

VecXDynamic Gaussian::gradient(const VecXDynamic& x) const{

    return - ((invSigma * (x - mu)) * eval(x));

}
