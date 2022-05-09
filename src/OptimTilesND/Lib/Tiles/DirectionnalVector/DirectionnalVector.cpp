#include "DirectionnalVector.hpp"

DirectionnalVector::DirectionnalVector(){
  dimension = 0;
  importantDimension = 0;
  importanceValue = 0.0;
}

DirectionnalVector::DirectionnalVector(int dim,int dimImp,double val){
  dimension = dim;
  importantDimension = dimImp;
  importanceValue = val;
}

DirectionnalVector::~DirectionnalVector(){

}

void DirectionnalVector::times(double valToMultiplyBy){
  importanceValue *= valToMultiplyBy;
}

int DirectionnalVector::getDim(){
  return dimension;
};
int DirectionnalVector::getImpDim(){
  return importantDimension;
};
double DirectionnalVector::getImpVal(){
  return importanceValue;
};

std::string DirectionnalVector::printDebug(){
  std::string toRet = "";
  toRet += "DirectionnalVector " + std::to_string(dimension) + "-D  :";
      for (int i = 0; i < dimension; i++) {
        if (i == importantDimension -1) {
          toRet += " " + std::to_string(i+1) + "-D ==> " + std::to_string(importanceValue) + ";";
        }else{
          toRet += " " + std::to_string(i+1) +  "-D ==> 0;";
        }
      }
      toRet += "\n";
    return toRet;
}
