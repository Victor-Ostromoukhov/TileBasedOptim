#include <iostream>
#include <string>
#include <vector>

class ConfigurationConfigFile
{
public:
  int noctave;
  bool owen_permut_flag; 
  int depth;
  int seed;
  bool dbg_flag;
  std::string output_fname;
  std::string input_matrices;
  int nbPtsToOptimize;
  int thread;
  int nIter;
  int nbReal;
  std::string output_points;
  int writingInterval;


  ConfigurationConfigFile();
  ~ConfigurationConfigFile();

  void initByFile(std::string);
  void printParameters();
};
