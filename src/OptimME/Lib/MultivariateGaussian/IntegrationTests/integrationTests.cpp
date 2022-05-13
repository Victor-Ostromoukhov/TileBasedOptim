//
// Created by lpaulin on 25/03/19.
//
#include <cmath>
#include <iostream>
#include "integrationTests.h"
#include <algorithm>


double gaussian(double x2, double v){
    return std::exp(-0.5 * x2 / (v*v)) / (v * std::sqrt(2 * M_PI)) ;
}

double gaussian(const VecXDynamic& p, const VecXDynamic& m, double v){
    VecXDynamic diff;
    if (m.dim() == 0){
        diff = p;
    } else{
        diff = p - m;
    }
    return gaussian(diff.norm2(), v);
}

double gaussianIntegration(const std::vector<VecXDynamic>& points, const VecXDynamic& m, double v){
    std::function<double(const VecXDynamic&)> f = [&m, &v](const VecXDynamic& p){ return gaussian(p, m, v); };
    return integration(points, f);
}
