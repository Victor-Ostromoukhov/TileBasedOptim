#include <iostream>
#include <vector>
#include "PointTernaire/PointTernaire.hpp"

class Tiles{
  public:
    PointTernaire samplingPoint;
    std::string afterSamplingPoint;
    Tiles();
    Tiles(PointTernaire,std::string);
    ~Tiles();
    friend std::ostream& operator<<(std::ostream& os, Tiles& t) {
            os << t.samplingPoint << "\t"<< t.afterSamplingPoint;
            return os;
        }
};

std::vector<Tiles>* importTiles(std::string,int);
void exportTiles(std::vector<Tiles>* ,std::string);

std::vector<PointTernaire>* extractSamplingPoints(std::vector<Tiles>*);
