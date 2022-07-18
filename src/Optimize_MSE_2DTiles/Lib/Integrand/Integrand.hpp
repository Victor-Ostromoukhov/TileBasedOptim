#include <functional>
#include "../Math/VecXDynamic.h"
#include "../Math/MatXDynamic.h"



class Integrand {
  public:
    VecXDynamic shift;
    MatXDynamic sigma;
    int dim;

    Integrand(const VecXDynamic& m, const MatXDynamic& sig):shift(m), dim(shift.dim()) {
      sigma = sig;
    }

    template <typename VECTYPE>
    inline double eval(const VECTYPE& p,int integrandType) const {
      double toReturn = 0.0;
      double result = 0.0;

      VecXDynamic pShift(2);
      for (int i = 0; i < p.dim(); i++) {
          pShift[i] = p[i];
      }
      pShift -= shift;
      switch (integrandType) {
        case 1: // HeaviSide
        result = 0.0;
        toReturn = 0.0;
        for (int idim = 0; idim < p.dim() ; idim++) {

              result += pShift[idim] * sigma.access(0,idim);
           }
        if (result > 0 ) toReturn = 1.;
        break;
        case 2: // SoftEllipses
            result = double((transpose(pShift) * sigma) * pShift);
            toReturn = (exp(-0.5 * result));
        break;
        case 3: // HardEllipses
            result = double((transpose(pShift) * sigma) * pShift);
            toReturn = ((exp(-0.5 * result)) <= 0.6065306597126334236038 ? 0 : 1);
        break;
        case 4: // SoftRectangles
            result = 1.0;
            toReturn = 1.0;
            for (int i = 0; i < p.dim(); i++) {
              result = ( (pShift[i] / sigma.access(0,i)) * (pShift[i] / sigma.access(0,i)) );
              toReturn *= exp(-.5 * result);
            }
        break;
        case 5: // HardRectangles
          result = 1.0;
          toReturn = 1.0;
          for (int i = 0; i < p.dim(); i++) {
            result = ( (pShift[i] / sigma.access(0,i)) * (pShift[i] / sigma.access(0,i)) );
            toReturn *= (exp(-.5 * result) <= 0.6065306597126334236038 ? 0 : 1);
          }
        break;
        default:
          std::cerr << "The chosen integrand type hasn't been implemented yet, so stay tuned." << '\n';
        break;
      }
      return toReturn;
    }
};

// Fonction qui permet de calculer l'intégration d'une gaussienne dont on nous a passé les paramètres

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

template <typename VECTYPE>
inline double multivariateGaussianIntegration(const std::vector<VECTYPE>& points, const VecXDynamic& m, const MatXDynamic& sigma,int integrandType){

    Integrand integrand(m, sigma);
    std::function<double(const VECTYPE&)> f = [&integrand,&integrandType](const VECTYPE& point){ return integrand.eval(point,integrandType); }; // Création de f,tel que f(xi) = yi

    return integration(points, f);

}

template <typename VECTYPE>
inline double multivariateGaussianPointValue(const VECTYPE& pointParam,const VecXDynamic& m, const MatXDynamic& sigma){

  Integrand integrand(m, sigma);

  std::function<double(const VECTYPE&)> f = [&integrand](const VECTYPE& point){ return integrand.eval(point); };

  return f(pointParam);
}

template <typename VECTYPE>
inline double multivariateGaussianIntegrationPointModif(const VECTYPE& oldPoint,const VECTYPE& newPoint,const VecXDynamic& m, const MatXDynamic& sigma,double oldValue,int nbpts,int integrandType){

  Integrand integrand(m, sigma);

  std::function<double(const VECTYPE&)> f = [&integrand,&integrandType](const VECTYPE& point){ return integrand.eval(point,integrandType); };
  return integrationPointVariation(oldPoint,newPoint,f,oldValue,nbpts);
}
