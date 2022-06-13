#include <iostream>
#include <vector>

class PointTernaire{
  public:
    int index;
    std::string xco;
    std::string yco;
    PointTernaire();
    PointTernaire(int,std::string,std::string);
    PointTernaire(int,int,int,int);
    ~PointTernaire();
    double getXcoAsDouble();
    double getYcoAsDouble();
    std::string getXcoAsConsommable();
    std::string getYcoAsConsommable();
    int ternaryStringToDecimalInt(std::string);
    std::string printAsDouble(bool);
    friend std::ostream& operator<<(std::ostream& os, PointTernaire& p) {
            os << std::to_string(p.index) << "\t"<< std::to_string(p.ternaryStringToDecimalInt((p.xco))) << "\t"<< std::to_string(p.ternaryStringToDecimalInt((p.yco)));
            return os;
        }
};

void exportPoints(std::vector<PointTernaire>* ,std::string);
std::string paddingNumber(std::string ,int ,bool);
std::string decimalIntToTernaryString(int);
