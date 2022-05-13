//
// Created by lpaulin on 20/05/19.
//

#ifndef STATTOOLSSAMPLING_MULTIVARIATEGAUSSIAN_H
#define STATTOOLSSAMPLING_MULTIVARIATEGAUSSIAN_H

#include <functional>
#include "IntegrationTests/integrationTests.h"

  // Classe de la gaussienne
class Gaussian {
    VecXDynamic mu;
    MatXDynamic invSigma;
    double normFactor;
    int dim;
public:

    Gaussian(const VecXDynamic& m, const MatXDynamic sigma, bool normalize = false);

    template <typename VECTYPE>
    inline double eval(const VECTYPE& p) const{
        VecXDynamic pShift(p);
        pShift -= mu;
        double result = double((transpose(pShift) * invSigma) * pShift);
        // std::cout << "Voilà pour vérifier" << '\n';
        // std::cout << "mu   " << mu[0]  << "  " << mu[1] << '\n';
        // std::cout << "point   "  << p[0] <<  "   "<< p[1] << '\n';
        // std::cout << "invSigma 1 > " << invSigma.access(0,0) << " 2 > " << invSigma.access(0,1) << " 3 > " << invSigma.access(1,0) << " 4 > " << invSigma.access(1,1) << '\n';
        // std::cout << "/* message */" <<  exp(-0.5 * result) << '\n';)
        std::cout << ((normFactor * exp(-0.5 * result)) <= 0.6065306597126334236038 ? 0 : 1) << '\n';
        return ((normFactor * exp(-0.5 * result)) <= 0.6065306597126334236038 ? 0 : 1);
    }

    template <typename VECTYPE>
    inline double eval(const VECTYPE& p,const int integrandType) const{
      double toReturn = 0.0;
      double result = 0.0;
      VecXDynamic pShift(p);
      pShift -= mu;
      switch (integrandType) {
        case 1:
          result = double((transpose(pShift) * invSigma) * pShift);
          // std::cout << "Voilà pour vérifier" << '\n';
          // std::cout << "mu   " << mu[0]  << "  " << mu[1] << '\n';
          // std::cout << "point   "  << p[0] <<  "   "<< p[1] << '\n';
          // std::cout << "invSigma 1 > " << invSigma.access(0,0) << " 2 > " << invSigma.access(0,1) << " 3 > " << invSigma.access(1,0) << " 4 > " << invSigma.access(1,1) << '\n';
          // std::cout << "/* message */" <<  exp(-0.5 * result) << '\n';)
          // std::cout << ((normFactor * exp(-0.5 * result)) <= 0.6065306597126334236038 ? 0 : 1) << '\n';
          toReturn = ((normFactor * exp(-0.5 * result)) <= 0.6065306597126334236038 ? 0 : 1);
        break;
        case 2:
            result = double((transpose(pShift) * invSigma) * pShift);
            toReturn = (normFactor * exp(-0.5 * result));
        break;
        case 3:
          result = 1.0;
          toReturn = 1.0;
          for (int i = 0; i < p.dim(); i++) {
            result = ( (pShift[i] / invSigma.access(0,i)) * (pShift[i] / invSigma.access(0,i)) );
            toReturn *= (exp(-.5 * result) <= 0.6065306597126334236038 ? 0 : 1);
          }
        break;
        case 4:
            result = 1.0;
            toReturn = 1.0;
            for (int i = 0; i < p.dim(); i++) {
              result = ( (pShift[i] / invSigma.access(0,i)) * (pShift[i] / invSigma.access(0,i)) );
              toReturn *= exp(-.5 * result);
            }
        break;
        default:
          std::cerr << "The chosen integrand type hasn't been implemented yet, so stay tune." << '\n';
        break;
      }
      return toReturn;
    }

    VecXDynamic gradient(const VecXDynamic& x) const;

};


// Fonction qui permet de calculer l'intégration d'une gaussienne dont on nous a passé les paramètres

template <typename VECTYPE>
inline double multivariateGaussianIntegration(const std::vector<VECTYPE>& points, const VecXDynamic& m, const MatXDynamic& sigma,int integrandType){

    Gaussian g(m, sigma, false);
    std::function<double(const VECTYPE&)> f = [&g,&integrandType](const VECTYPE& point){ return g.eval(point,integrandType); }; // Création de f,tel que f(xi) = yi

    return integration(points, f);

}

template <typename VECTYPE>
inline double multivariateGaussianPointValue(const VECTYPE& pointParam,const VecXDynamic& m, const MatXDynamic& sigma){

  Gaussian g(m, sigma, false);

  std::function<double(const VECTYPE&)> f = [&g](const VECTYPE& point){ return g.eval(point); };

  return f(pointParam);
}

template <typename VECTYPE>
inline double multivariateGaussianIntegrationPointModif(const VECTYPE& oldPoint,const VECTYPE& newPoint,const VecXDynamic& m, const MatXDynamic& sigma,double oldValue,int nbpts,int integrandType){

  Gaussian g(m, sigma, false);

  std::function<double(const VECTYPE&)> f = [&g,&integrandType](const VECTYPE& point){ return g.eval(point,integrandType); };
  return integrationPointVariation(oldPoint,newPoint,f,oldValue,nbpts);
}



#endif //STATTOOLSSAMPLING_MULTIVARIATEGAUSSIAN_H
