//
// Created by lpaulin on 25/03/19.
//

#ifndef STATTOOLSSAMPLING_INTEGRATIONTESTS_H
#define STATTOOLSSAMPLING_INTEGRATIONTESTS_H

#include <vector>
#include <functional>
#include "../../Math/VecXDynamic.h"
#include "../../Math/MatXDynamic.h"
#include <algorithm>

template <typename VECTYPE>
double integration(const std::vector<VECTYPE>& points, std::function<double(const VECTYPE&)> f){

    double val =0.;
    for (int i = 0; i < (int)points.size(); i++) {
      val += f(points[i]);
    }
    val /= (double)points.size();

    return val;
}

template <typename VECTYPE>
double integrationPointVariation(const VECTYPE& oldPoint,const VECTYPE& newPoint, std::function<double(const VECTYPE&)> f,double oldValue,int nbpts){
  return oldValue - f(oldPoint)/nbpts + f(newPoint)/nbpts;
}

#endif //STATTOOLSSAMPLING_INTEGRATIONTESTS_H
