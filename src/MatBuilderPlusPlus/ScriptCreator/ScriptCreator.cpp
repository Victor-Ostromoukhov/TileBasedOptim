#include <string>
#include "include/Configurations/ConfigurationConfigFile/ConfigurationConfigFile.hpp"
#include "include/Configurations/ConfigurationMB/ConfigurationMB.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstring>
#include <random>
//#include <omp.h>
#include <chrono>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <ios>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <vector>
#include <cstring>
#include <fstream>
#include "../include/CLI11.hpp"
#include "MatrixTools.h"
#include "MatrixSamplerClass.h"
#include "Scrambling.h"
#include <filesystem>
namespace fs = std::filesystem;

void initializeConfigFiles(ConfigurationConfigFile* configFile, std::string pathToConfigFile, ConfigurationMB* configMB, std::string pathToMBProfile){
  configFile->initByFile(pathToConfigFile);
  configMB->initByFile(pathToMBProfile);
}

void writeScript(ConfigurationConfigFile* configFile, ConfigurationMB* configMB, std::string scriptLocation){
  std::ofstream os;
  os.open(scriptLocation,std::ios_base::trunc);

  std::string list ="lst=(3 ";

  int val = 3;
  while (val*3 <= configFile->nbPtsToOptimize) {
    list += std::to_string(val * 2)+" ";
    val *= 3;
    list += std::to_string(val)+" ";
  }
  list+=")\n";

  os << "#!/bin/bash \n\n";

  os << "nbthreads="+std::to_string(configFile->thread)+"\n";// -t configurationConfigFile->nbthreads

  os << "nIterations="+std::to_string(configFile->nIter)+"\n";// -n configurationConfigFile->nIter

  os << "ifname=\""+configFile->output_points+"\"\n";// -i configurationConfigFile->output_points

  os << "ofname=\""+configFile->output_points+"\"\n";// -o configurationConfigFile->output_points

  os << "innoctave="+std::to_string(configFile->noctave)+"\n";//  --innoctave configurationConfigFile->noctave

  os << "writingInterval="+std::to_string(configFile->writingInterval)+"\n";  // --writingInterval configurationConfigFile->writingInterval

  os << "base="+std::to_string(configMB->base)+"\n";

  os << "nbTot="+std::to_string(configFile->nbPtsToOptimize)+"\n";

  os << list;

  os << "lst_length=${#lst[@]}\n\n\n\n";

  for (int nbOfOptim = 0; nbOfOptim < (int) configMB->ConfigVector.size(); nbOfOptim++) {
    os << "integrandType="+std::to_string(configMB->ConfigVector.at(nbOfOptim).inteType)+"\n"; // --integrandType configMB->ConfigVector->at(i).inteType

    os << "nItegrandsPerIteration="+std::to_string(configMB->ConfigVector.at(nbOfOptim).batch)+"\n";// -g configMB->ConfigVector->at(i).batch

    if (configMB->ConfigVector.at(nbOfOptim).previous) {
      os << "previous=true\n";
    }else{
      os << "previous=false\n";
    }

    os << "for level in $(eval echo {1..$((${lst_length} - 1))})\ndo\n";

    os << " to=${lst[$level]}\n limit=$((${lst[$((${level} - 1 ))]} + 1))\n";

    os << "   ~/bin/OptimMSE2DPipeLine -t ${nbthreads} -n $nIterations -i $ifname -o $ofname --nbPoints $to --integrandType $integrandType -g $nItegrandsPerIteration --limit $limit --writingInterval $writingInterval --base $base --previous $previous --innoctave $innoctave --nbTotal $nbTot";

    os << " --dimensionToOptimize";

    for (int dimToOptimize = 0; dimToOptimize < (int) configMB->ConfigVector.at(nbOfOptim).dimToOpt.size() ; dimToOptimize++) {
      os << " "+std::to_string(configMB->ConfigVector.at(nbOfOptim).dimToOpt.at(dimToOptimize));
    }

    os << "\ndone\n\n";

  }
  os.close();
  fs::permissions(scriptLocation,fs::perms::owner_all | fs::perms::others_exec,fs::perm_options::add);
}





void fullFillPointHolder(double** pointHolder, ConfigurationConfigFile* configurationConfigFile, ConfigurationMB* configurationMB){
  std::ifstream in(configurationConfigFile->input_matrices);

  std::string input_matrices = configurationConfigFile->input_matrices;
  int nDims = configurationMB->nDims;
  int m = configurationMB->m;

  if (in.fail()) {
    std::cerr << "Error: Could not open input file: " << input_matrices << std::endl;
    exit(-1);
  }

  std::vector<std::vector<int> > Bs(nDims, std::vector<int>(m*m));
  std::vector<std::vector<int> > Cs(nDims, std::vector<int>(m*m));

  readMatrices(in, m, nDims, Cs);

  std::minstd_rand gen(configurationConfigFile->seed);
  std::uniform_int_distribution<int> unif;

  int nbReal = configurationConfigFile->nbReal;
  for (int real = 0; real < nbReal; ++real) {
    int real_seed = unif(gen);
    for (int indpt = 0; indpt < configurationConfigFile->nbPtsToOptimize ; ++indpt) {
      for (int inddim = 0; inddim < nDims; ++inddim) {
        double pos;
        if (configurationConfigFile->owen_permut_flag){
          pos = getScrambledDouble(Cs[inddim], m, configurationMB->base, indpt, real_seed + inddim, configurationConfigFile->depth);
          pointHolder[indpt][inddim] = pos;
        } else {
          pos = getDouble(Cs[inddim], m, configurationMB->base, indpt);
          pointHolder[indpt][inddim] = pos;
        }
      }
    }
  }

}


void writePoints(double** pointHolder,std::string outputFileName,int nbPtsToOptimize, int nDims){
  std::ofstream out;
  out.open(outputFileName,std::ios_base::app);
  out << std::setprecision(17);
  for (int i = 0; i < nbPtsToOptimize; i++) {
    for (int j = 0; j < nDims; j++) {
     out << pointHolder[i][j] << '\t';
    }
    out << std::endl;
  }

  out.close();

}




int main(int argc, char const *argv[]) {
  ConfigurationConfigFile* configurationConfigFile = new ConfigurationConfigFile();
  ConfigurationMB* configurationMB = new ConfigurationMB();


  double **pointHolder = new double*[configurationConfigFile->nbPtsToOptimize];
  for (int i = 0; i < configurationConfigFile->nbPtsToOptimize; ++ i) {
  pointHolder[i] = new double[configurationMB->nDims];
  }
  std::string configFilePath = "../../data/config.txt";
  std::string MBProfilePath = "../../data/profile.txt";
  std::string scriptLocation = "../../Scripts/LaunchScript.sh";

  /* ------------- Configuration de CLI11 --------------- */

      CLI::App app { "ScriptCreator, a util that reads a MB profile and a configuration file in order to write and launch the scripts to realize the optimisations " };

      app.add_option("--configFilePath",configFilePath,"Path to the configuration file. Default: "+configFilePath)->check(CLI::ExistingFile);
      app.add_option("--MBProfilePath",MBProfilePath,"Path to the profile. Default: "+MBProfilePath)->check(CLI::ExistingFile);
      app.add_option("--ScriptLocation",scriptLocation,"Path for the script to output. Default: "+scriptLocation);
      CLI11_PARSE(app, argc, argv)

  /* ------------- Initisalization of the configuration objects --------- */

  if (!(fs::exists(configFilePath) && fs::exists(MBProfilePath)) ) { // If the default paths don't exist either
    exit(-1);
  }

  initializeConfigFiles(configurationConfigFile,configFilePath,configurationMB,MBProfilePath);

  /* ------------- Writing of the scripts ----------- */

  writeScript(configurationConfigFile,configurationMB,scriptLocation);

  /* ------------- Filling the point array by reading MB matrices ---------- */

  fullFillPointHolder(pointHolder,configurationConfigFile,configurationMB);

  /* ------------- Writing of the points that will be modified by the optimizing scripts --------- */

  writePoints(pointHolder, configurationConfigFile->output_points,  configurationConfigFile->nbPtsToOptimize, configurationMB->nDims);

  return 0;
}
