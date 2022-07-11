#include <iostream>
#include "ConfigurationConfigFile.hpp"
#include <fstream>
#include <sstream>

ConfigurationConfigFile::ConfigurationConfigFile(){

  noctave = 10;
  owen_permut_flag = true;
  depth = 19;
  seed = rand() %  65537;
  dbg_flag = false;
  output_fname = "tmp/pts.dat";
  input_matrices = "MatBuilder_matrices/2D_0m2net_000001.dat";
  nbPtsToOptimize = 729;
  thread = 8;
  nIter = 400;
  nbReal = 1;
  output_points = "output_points.dat";
  writingInterval = 1000;
}

ConfigurationConfigFile::~ConfigurationConfigFile(){

}

void ConfigurationConfigFile::printParameters(){
  std::cout << "noctave : "<< noctave << std::endl;
  std::cout << "owen_permut_flag : " << owen_permut_flag << std::endl;
  std::cout << "depth : " << depth << std::endl;
  std::cout << "seed : " << seed << std::endl;
  std::cout << "dbg_flag : "<< dbg_flag << std::endl;
  std::cout << "output_fname : " << output_fname << std::endl;
  std::cout << "input_matrices : " << input_matrices << std::endl;
  std::cout << "nbPtsToOptimize : " << nbPtsToOptimize << std::endl;
  std::cout << "thread : " << thread << std::endl;
  std::cout << "nIter : " << nIter << std::endl;
  std::cout << "nbReal : " << nbReal << std::endl;
  std::cout << "output_points : " << output_points << std::endl;
  std::cout << "writingInterval : " << writingInterval << std::endl;
}


void ConfigurationConfigFile::initByFile(std::string configFile){
  std::ifstream in(configFile,std::ios::in);
  std::string line;
  if (in.is_open()) {
    while (std::getline(in, line)) {
      std::istringstream sline(line);
      char letter;
      sline >> letter;
      if (letter == '#') {
        continue;
      }
      sline.putback(letter);
      std::string full_left_part;
      getline(sline, full_left_part, '=');
      if (!sline.eof()) {
        std::string left_part;
        std::istringstream left_ss(full_left_part);
        left_ss >> left_part;
        if (left_part == "noctave") {
          sline >> noctave;
        } else if (left_part == "depth") {
          sline >> depth;
        }else if (left_part == "nbPtsToOptimize") {
          sline >> nbPtsToOptimize;
        }else if (left_part == "thread") {
          sline >> thread;
        }else if (left_part == "seed") {
          sline >> seed;
        }else if (left_part == "nIter") {
          sline >> nIter;
        }else if (left_part == "owen_permut_flag") {
          std::string flag;
          sline >> flag;
          if (flag == "false") {
            owen_permut_flag = false;
          } else if (flag == "true"){
            owen_permut_flag = true;
          }
        }else if (left_part == "dbg_flag") {
          std::string flag;
          sline >> flag;
          if (flag == "false") {
            dbg_flag = false;
          } else if (flag == "true"){
            dbg_flag = true;
          }
        }else if (left_part == "output_fname") {
          sline >> output_fname;
        }else if (left_part == "input_matrices") {
          sline >> input_matrices;
        } else if (left_part == "nbReal" ){
        sline >> nbReal;
        } else if (left_part == "output_points") {
          sline >> output_points;
        } else if (left_part == "writingInterval") {
          sline >> writingInterval;
        }
      }
    }
  } else {
    std::cout << "File cannot be open : check if the path specified was correct." << std::endl;
  }
  in.close();
}
