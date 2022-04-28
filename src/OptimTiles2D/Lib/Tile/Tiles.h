#include <string>
#include "Points_2D/Points_2D.h"
class Tiles{
  private:
    std::string beforeSamplingPoint;
    Points_2D samplingPoint;
    Points_2D previousRefPoint;
    Points_2D previousv1;
    Points_2D previousv2;
    std::string afterPreviousv2;
  public:
    Tiles();
    Tiles(std::string,Points_2D,Points_2D,Points_2D,Points_2D,std::string);
    ~Tiles();
    Points_2D getSamplingPoint();
    void setSamplingPoint(Points_2D);
    Points_2D getPreviousRefPoint();
    Points_2D getPreviousv1();
    Points_2D getPreviousv2();
    friend std::ostream& operator<<(std::ostream& os, Tiles& t) {

            os << t.beforeSamplingPoint << "\t" << t.samplingPoint<< "\t" << t.previousRefPoint << "\t"<< t.previousv1<< "\t" << t.previousv2 << "\t"<< t.afterPreviousv2;
            return os;
        }
};
