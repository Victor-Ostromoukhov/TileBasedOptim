#include <iostream>
#include <vector>

class ConfigurationMB
{
public:
  enum InteType {SoftEllipse = 2,Heaviside = 1};
  struct OneConfig{
    InteType inteType;
    bool previous;
    int batch;
    std::vector<int> dimToOpt;
    int tileDim;
  };
  int base;
  int nDims;
  int m;
  std::vector<OneConfig> ConfigVector;

  ConfigurationMB();
  ~ConfigurationMB();
  void initByFile(std::string);
  void printParameters();
};
