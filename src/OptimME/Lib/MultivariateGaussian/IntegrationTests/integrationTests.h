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
// int bnOp = 0;
//     static const size_t nbThreads = 1;
//
//     std::vector<double> res(nbThreads, 0.);
//
//     std::function<double(const VECTYPE&)> fcopy[nbThreads];
//     for (size_t i = 0; i < nbThreads; ++i) {
//         fcopy[i] = f;
//     }
//
// #pragma omp parallel for shared(res)
//     for (int i = 0; i < nbThreads; ++i){
//         for (size_t j = i * points.size() / nbThreads; j < std::min(points.size(), (i+1) * points.size() / nbThreads); ++j){
//             res[i] += fcopy[i](points[j]);
//             bnOp += 1;
//         }
//     }
//
//     double val = 0.;
//     for (size_t j = 0; j < nbThreads; ++j){
//         val += res[j];
//         bnOp += 1;
//     }
    double val =0.;
    for (int i = 0; i < points.size(); i++) {
      // std::cout << "Les mioens" <<f(points[i]) << '\n';
      // std::cout << "Ceci est le pt passé à la fonction qui pas fctone :)"  << points[i][0] << "   " << points[i][1] << '\n';

      val += f(points[i]);
    }
    // std::cout << "Evolution de accumulator "  << val << '\n';

    // std::cout << val << " / " << points.size() << '\n';
    val /= (double)points.size();

    return val;
}

template <typename VECTYPE>
double integrationPointVariation(const VECTYPE& oldPoint,const VECTYPE& newPoint, std::function<double(const VECTYPE&)> f,double oldValue,int nbpts){
	std::cout << newPoint[0] << " " << newPoint[1] <<  " -> " << f(newPoint) << std::endl;
  return oldValue - f(oldPoint)/nbpts + f(newPoint)/nbpts;
}

#endif //STATTOOLSSAMPLING_INTEGRATIONTESTS_H
