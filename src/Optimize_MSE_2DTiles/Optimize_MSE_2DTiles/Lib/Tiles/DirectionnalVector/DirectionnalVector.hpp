#include <iostream>
class DirectionnalVector
{
  private:
    int dimension;
    int importantDimension;
    double importanceValue;
  public:
    DirectionnalVector();
    DirectionnalVector(int,int,double);
    ~DirectionnalVector();
    void times(double);
    int getDim();
    int getImpDim();
    double getImpVal();
    std::string printDebug();
    friend std::ostream& operator<<(std::ostream& os, DirectionnalVector dV) {
      for (int i = 0; i < dV.getDim(); i++) {
        if (i == dV.getImpDim()-1) {
          os << dV.getImpVal() << "\t";
        }else{
          os << "0 \t";
        }
      };
      return os;
  }
};
