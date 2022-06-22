#include <string>
#include <fstream>
#include <cstdarg>
#include <vector>
#include "array"
#include "../Math/VecXDynamic.h"

template<int dimension>
class Tiles
{
  public:
    std::string beforeSamplingPoint;
    VecX<dimension> samplingPoint;
    VecX<dimension> previousRefPoint;
    std::array<VecX<dimension>,dimension> previousDirectionnalVectors;
    std::string afterPreviousDirectionnalVectors;

    inline Tiles();
    inline ~Tiles();

    inline  VecX<dimension>* getSamplingPoint();
    inline void setSamplingPoint(VecX<dimension>);

    inline VecX<dimension>* getPreviousRefPoint();
    inline int getDim();


    inline VecX<dimension> getPreviousDirectionnalVector(int) const;
    inline void setPreviousDirectionnalVector(int,VecX<dimension>);

    inline std::string getBeforeSP();
    inline void setBeforeSP(std::string);

    inline std::string getAfterDV();
    inline void setAfterDV(std::string);


    inline void printDebug();
    friend std::ostream& operator<<(std::ostream& os, Tiles t) {
      // before, all points, all directionnal vectores after
      os << t.getBeforeSP()  << *(t.getSamplingPoint()) << '\t' <<*(t.getPreviousRefPoint()) << '\t';
      for (int i = 0; i < t.getDim(); i++) {
        os << t.getPreviousDirectionnalVector(i+1) << " ";
      }
      os << t.getAfterDV() << "\n";
      return os;
    }
};

// template<int dimension>
// inline void Tiles<dimension>::vectorPointSumRand(VecX<dimension>* point,VecX<dimension> refPoint,double valToMultiplyBy,int indice,DirectionnalVector dV){
//   if (dimension == 2) {
//     int delta;
//     if (dV.getImpDim() == 1) {
//       delta = indice / 8;
//     }else{
//       delta = indice % 8;
//     }
//     (*point).coefs[dV.getImpDim()-1] = (refPoint[dV.getImpDim()-1] +dV.getImpVal()*(delta + valToMultiplyBy) / 8);
//   }else{
//     if (dimension == 3) {
//       int delta;
//       if (dV.getImpDim() == 1) {
//         delta = indice / 16;
//       }else{
//         if (dV.getImpDim() == 2) {
//           delta = ((indice/4)%4);
//         }else{
//           delta = indice % 4;
//         }
//       }
//       (*point)[dV.getImpDim()-1] = (refPoint[dV.getImpDim()-1] + dV.getImpVal()*(delta + valToMultiplyBy) / 4);
//     }else{
//       (*point)[dV.getImpDim()-1] = (refPoint[dV.getImpDim()-1]+dV.getImpVal()*valToMultiplyBy);
//     }
//   }
// }

template<int dimension>
Tiles<dimension> parseLine(int,std::string);

template<int dimension>
std::vector<Tiles<dimension>>* importTiles(std::string,int,std::vector<std::string>*);

template<int dimension>
inline Tiles<dimension>::Tiles(void){}

template<int dimension>
inline Tiles<dimension>::~Tiles(){

};

template<int dimension>
inline void Tiles<dimension>::setSamplingPoint(VecX<dimension> newSP){
  samplingPoint = newSP;
}

template<int dimension>
inline VecX<dimension>* Tiles<dimension>::getSamplingPoint(){
  return &samplingPoint;
}

template<int dimension>
inline int Tiles<dimension>::getDim(){
  return dimension;
}

template<int dimension>
inline VecX<dimension>* Tiles<dimension>::getPreviousRefPoint(){
  return &previousRefPoint;
}

template<int dimension>
inline VecX<dimension> Tiles<dimension>::getPreviousDirectionnalVector(int Dim) const{
  return previousDirectionnalVectors[Dim-1];
}

template<int dimension>
inline void Tiles<dimension>::setPreviousDirectionnalVector(int index,VecX<dimension> dV){
  previousDirectionnalVectors[index] = dV;
}

template<int dimension>
inline void Tiles<dimension>::printDebug(){
  std::cout << "Tuile " << dimension << "-D : \n Before Sampling Point : " << beforeSamplingPoint << "\t  Sampling Point : " << samplingPoint.printDebug() << "\t previousRefPoint : " << previousRefPoint.printDebug();
  for (int i = 0; i < dimension; i++) {
    std::cout << getPreviousDirectionnalVector(i+1).printDebug();
  }
   std::cout << " afterPreviousDirectionnalVectors : " << afterPreviousDirectionnalVectors;
  }

template<int dimension>
inline std::string Tiles<dimension>::getBeforeSP(){
    return beforeSamplingPoint;
}

template<int dimension>
inline void Tiles<dimension>::setBeforeSP(std::string bef){
  beforeSamplingPoint = bef;
}

template<int dimension>
inline std::string Tiles<dimension>::getAfterDV(){
    return afterPreviousDirectionnalVectors;
}

template<int dimension>
inline void Tiles<dimension>::setAfterDV(std::string aFDV){
  afterPreviousDirectionnalVectors = aFDV;
}

template<int dimension>
Tiles<dimension> parseLine(std::string  lineToParse){
  std::string afterPreviousDirectionnalVectors = "";
  double val;
  Tiles<dimension> tileToRet;
  size_t pos = 0;
  std::string token;
  // ======= beforeSamplingPoint ======= //
      // === TileType === //
      // pos = lineToParse.find("\t");
      // tileToRet.setBeforeSP(tileToRet.getBeforeSP()+lineToParse.substr(0,pos)+"\t");
      // lineToParse.erase(0, pos + 1);

      // === Structural ID === //
      pos = lineToParse.find("\t");
      tileToRet.setBeforeSP(tileToRet.getBeforeSP()+lineToParse.substr(0,pos)+"\t");
      lineToParse.erase(0, pos + 1);

      for (int sPComponent = 0; sPComponent < dimension; sPComponent++) {
        pos = lineToParse.find("\t");
        tileToRet.getSamplingPoint()->coefs[sPComponent] = std::stod(lineToParse.substr(0,pos));
        lineToParse.erase(0, pos + 1);
      }

      for (int sPComponent = 0; sPComponent < dimension; sPComponent++) {
        pos = lineToParse.find("\t");
        tileToRet.getPreviousRefPoint()->coefs[sPComponent] = std::stod(lineToParse.substr(0,pos));
        lineToParse.erase(0, pos + 1);
      }

      for (int dVNb = 0; dVNb < dimension; dVNb++) {
        for (int dVComponent = 0; dVComponent < dimension; dVComponent++) {
          pos = lineToParse.find("\t");
          val = std::stod(lineToParse.substr(0,pos));
          tileToRet.previousDirectionnalVectors[dVNb][dVComponent] = val;
          lineToParse.erase(0, pos + 1);
        }
      }


      tileToRet.setAfterDV(/*tileToRet.getAfterDV()+lineToParse.substr(0,std::string::npos)*/"");
      lineToParse.erase(0, std::string::npos);
      // std::cout << tileToRet << '\n';
      return tileToRet;
};

template<int dimension>
std::vector<Tiles<dimension>>* importTiles(std::string fileToRead,int nbpts,std::vector<std::string>* restOfTheDocument){
  std::ifstream tilesFile(fileToRead, std::ios::in);
  std::vector<Tiles<dimension>>* vectorOfTiles = new std::vector<Tiles<dimension>>;

  if (tilesFile.is_open()) {
    std::string line;
    for (int currentPoint = 0; currentPoint < nbpts; currentPoint++) {
        std::getline(tilesFile,line);
        vectorOfTiles->push_back(parseLine<dimension>(line));
    }
    while (std::getline(tilesFile,line)) {
      restOfTheDocument->push_back(line);
    }
    tilesFile.close();
  }else{
    std::cerr << "Unable to open file\n";
  }
  return vectorOfTiles;
}
