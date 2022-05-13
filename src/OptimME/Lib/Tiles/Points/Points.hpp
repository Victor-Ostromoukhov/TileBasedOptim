#include <iostream>
#include <vector>
#include <array>
template<int dimension>
class Points
{
  public:
    std::array<double,dimension> dimensionnalArray;
    inline Points();
    inline int getDim();
    inline std::string printDebug();
    friend std::ostream& operator<<(std::ostream& os, Points p) {
        for (int i = 0; i < p.getDim(); i++) {
          os << p.dimensionnalArray[i] << "\t";
        }
        return os;
  }
};
template<int dimension>
inline Points<dimension>::Points(void){
  for (int i = 0; i < dimension; i++) {
    dimensionnalArray[i] = 0.0;
  }
}

template<int dimension>
inline int Points<dimension>::getDim(){
  return dimension;
};

template<int dimension>
inline std::string Points<dimension>::printDebug(){
  std::string toRet = "";
  toRet +=  "Point " + std::to_string(dimension) + "-D :";
  for (int i = 0; i < dimension; i++) {
    toRet+= " " + std::to_string(i+1) + "-D ==> " + std::to_string(dimensionnalArray[i]) + ";";
  }
  toRet +="\n";
  return toRet;
}
