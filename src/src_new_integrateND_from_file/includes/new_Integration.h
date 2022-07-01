//
// Created by lois on 18/09/2020.
//

#ifndef SOBOL_GAUSSIANS_H
#define SOBOL_GAUSSIANS_H

#include <cmath>
#include <vector>
#include <iostream>
#include "VecX.h"
#include "Random.h"
#include "drand48.h"

struct t_Heaviside2D {double integral; double muDiscotinuity[2]; double normal[2]; };
struct t_Heaviside3D {double integral; double muDiscotinuity[3]; double normal[3]; };
struct t_Heaviside4D {double integral; double muDiscotinuity[4]; double normal[4]; };
struct t_Heaviside5D {double integral; double muDiscotinuity[5]; double normal[5]; };
struct t_Heaviside6D {double integral; double muDiscotinuity[6]; double normal[6]; };


struct t_GaussianStruct1D {double integral; double mu[1] ; double mxCInv[1 * 1] ; };
struct t_GaussianStruct2D {double integral; double mu[2] ; double mxCInv[2 * 2] ; };
struct t_GaussianStruct3D {double integral; double mu[3] ; double mxCInv[3 * 3] ; };
struct t_GaussianStruct4D {double integral; double mu[4] ; double mxCInv[4 * 4] ; };
struct t_GaussianStruct5D {double integral; double mu[5] ; double mxCInv[5 * 5] ; };
struct t_GaussianStruct6D {double integral; double mu[6] ; double mxCInv[6 * 6] ; };
struct t_GaussianStruct8D {double integral; double mu[8] ; double mxCInv[8 * 8] ; };
struct t_GaussianStruct10D {double integral; double mu[10] ; double mxCInv[10 * 10] ; };
struct t_GaussianStruct12D {double integral; double mu[12] ; double mxCInv[12 * 12] ; };

#include "../Integrands/Heaviside2D_nIntegrands524288_testSet.hpp"
#include "../Integrands/SoftEllipses2D_testSet_4x4x128x128_testSet.hpp"
//#include "../Integrands/SoftEllipses2D_2x256x256_optimSet.hpp"

//extern t_Heaviside2D tab_Heaviside2D[262144] ;
extern t_Heaviside3D tab_Heaviside3D[16384] ;
extern t_Heaviside4D tab_Heaviside4D[16384] ;
extern t_Heaviside5D tab_Heaviside5D[16384] ;
extern t_Heaviside6D tab_Heaviside6D[16384] ;

inline double getMultivariateGaussian(const int nDims, const double *pt_ND, const double *mu, const double *mxCInv) {
    double accumulator = 0.;
    for (int row = 0; row < nDims; row++) {
        for (int col = 0; col < nDims; col++)
        	accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
    }
    return exp(-.5 * accumulator);
}	// getMultivariateGaussian

inline double getHardMultivariateGaussian(const int nDims, const double *pt_ND, const double *mu, const double *mxCInv) {
    double accumulator = 0.;
    for (int row = 0; row < nDims; row++) {
        for (int col = 0; col < nDims; col++)
        	accumulator += (pt_ND[row] - mu[row]) * (pt_ND[col] - mu[col]) * mxCInv[col*nDims+row];
    }
    return (exp(-0.5 * accumulator) <= 0.6065306597126334236038 ? 0 : 1);
}	// getMultivariateGaussian

//getGaussian1D[pt1D_,mu1D_,k_] := Quiet[Exp[-.5  ((pt1D-mu1D)/k)^2]]
//getSoftRectanglesND[ptND_,muND_,ktab_] := Product[getGaussian1D[ptND[[iDim]],muND[[iDim]],ktab[[iDim]] ],{iDim,Length[ptND]}]
//
//getRectangles1D[pt1D_,mu1D_,k_] := If[getGaussian1D[pt1D,mu1D,k] > 0.6065306597126334236038, 1, 0]  (* E^(-1/2) // N *)
//getRectanglesND[ptND_,muND_,ktab_] := Product[getRectangles1D[ptND[[iDim]],muND[[iDim]],ktab[[iDim]] ],{iDim,Length[ptND]}]

inline double getGaussian1D(const double pt1D, const double mu1D, const double sigma1D) {
	double tmp = (pt1D-mu1D)/sigma1D;
    return exp(-.5 * tmp*tmp);
}	// getMultivariateGaussian

inline double getRectanle1D(const double pt1D, const double mu1D, const double sigma1D) {
	double tmp = (pt1D-mu1D)/sigma1D;
    return (exp(-0.5 * tmp*tmp) <= 0.6065306597126334236038 ? 0 : 1);
}	// getMultivariateGaussian


inline double getSoftRectanglesND(const int nDims, const double *pt_ND, const double *mu, const double *sigmaND) {
    double accumulator = 1;
    for (int iDim = 0; iDim < nDims; iDim++) {
        	accumulator *= getGaussian1D(pt_ND[iDim],mu[iDim],sigmaND[iDim]);
    }
    return accumulator;
}	// getMultivariateGaussian

inline double getRectanglesND(const int nDims, const double *pt_ND, const double *mu, const double *sigmaND) {
    double accumulator = 1;
    for (int iDim = 0; iDim < nDims; iDim++) {
        	accumulator *= getRectanle1D(pt_ND[iDim],mu[iDim],sigmaND[iDim]);
    }
    return accumulator;
}	// getMultivariateGaussian

template <typename Vec>
inline double calculate_mse(const std::vector<Vec>& points,
                            const int integrandType = 1,
                            const int nintegrands = 1024,
                            const int seed=-1)
{
    RNG rng;
    if (seed != -1){
        rng.seed(seed);
    } else {
        rng.seed((1u<<31u) * drand48());
    }

    double mse_accumulator = 0.;
    int nDims = points[0].dim();
    // integration over 1K integrands over 16K
//#pragma omp parallel for schedule(dynamic) reduction(+:mse_accumulator)
#pragma omp parallel for reduction(+:mse_accumulator)
    for (int iintegrands = 0; iintegrands < nintegrands; iintegrands ++) {
        double integral, *mu, *mxCInv, *sigma,  *muDiscotinuity, *normal;
        int integrand_index;
		double *muleft,*muright,*mubottom,*mutop,*vleft,*vright,*vbottom,*vtop;

		integrand_index = iintegrands;
        double accumulator = 0.;
        switch (integrandType) {
		case 1 : // HeavisideND
            	// getHeaviside2D[{x_,y_},{muDiscotinuity_,normVector_},mulFactor_:1] := If[({x,y}-muDiscotinuity).normVector > 0, 1, 0]
        		switch(nDims) {
            	case 2 : // Heaviside2D
                  	integral = tab_Heaviside2D[integrand_index].integral;
                  	muDiscotinuity = tab_Heaviside2D[integrand_index].muDiscotinuity;
                  	normal = tab_Heaviside2D[integrand_index].normal;
                  	break;
//            	case 3 : // Heaviside3D
//                  	integral = tab_Heaviside3D[integrand_index].integral;
//                  	muDiscotinuity = tab_Heaviside3D[integrand_index].muDiscotinuity;
//                  	normal = tab_Heaviside3D[integrand_index].normal;
//            		break;
//            	case 4 : // Heaviside4D
//                  	integral = tab_Heaviside4D[integrand_index].integral;
//                  	muDiscotinuity = tab_Heaviside4D[integrand_index].muDiscotinuity;
//                  	normal = tab_Heaviside4D[integrand_index].normal;
//                  	break;
//            	case 5 : // Heaviside5D
//                  	integral = tab_Heaviside5D[integrand_index].integral;
//                  	muDiscotinuity = tab_Heaviside5D[integrand_index].muDiscotinuity;
//                  	normal = tab_Heaviside5D[integrand_index].normal;
//            		break;
//            	case 6 : // Heaviside6D
//                  	integral = tab_Heaviside6D[integrand_index].integral;
//                  	muDiscotinuity = tab_Heaviside6D[integrand_index].muDiscotinuity;
//                  	normal = tab_Heaviside6D[integrand_index].normal;
//            		break;
                default: std::cerr << "nDims=" << nDims << ": not implemented." << std::endl;
                      exit(1);
                      break;
        		}
                for (int ipt = 0; ipt < points.size() ; ipt++) {
                	double dotProduct = 0;
                	Vec pt = points[ipt];
                     for (int idim = 0; idim < nDims ; idim++) {
                       	double integrand_mu_componenet = muDiscotinuity[idim];
                       	double integrand_normal_componenet = normal[idim];
                       	dotProduct += (pt[idim] - integrand_mu_componenet) * integrand_normal_componenet;
                     }
                     if (dotProduct > 0 ) accumulator += 1.;
                }
                break;
    	case 2 :	// SoftEllipses
               	// getMultivariate2D[{x_,y_}, {mu_,mxCInv_},mulFactor_:1] := Quiet[ 1./mulFactor Exp[-.5 ({x,y}-mu).mxCInv.({x,y}-mu)] ]
        		switch(nDims) {
            	case 2 : // Gauss 2D
                	integral = tab_SoftEllipses2D[integrand_index].integral;
                    mu = tab_SoftEllipses2D[integrand_index].mu;
                    mxCInv = tab_SoftEllipses2D[integrand_index].mxCInv;
                  	break;
                default: std::cerr << "nDims=" << nDims << ": not implemented." << std::endl; exit(1); break;
        		}
                for (const Vec & point : points) {
                	accumulator += getMultivariateGaussian(nDims, point.data(), mu, mxCInv);
                }
                break;
        default: std::cout << "integrandType " << integrandType << " not implemented." << std::endl;
                exit(1);
                break;
        }
        accumulator /= (double)points.size();
        double mse = (integral-accumulator)*(integral-accumulator);
        mse_accumulator += mse;
    }
    return mse_accumulator / (double)(nintegrands);
} // calculate_mse


#endif //SOBOL_GAUSSIANS_H
