#include <iostream>
#include "ConfigurationMB.hpp"
#include <fstream>
#include <sstream>

ConfigurationMB::ConfigurationMB(){
  OneConfig confDefault;
  confDefault.inteType = SoftEllipse;
  nDims = 2;
  base = 3;
  m = 19;
  confDefault.batch = 4096;
  confDefault.dimToOpt.push_back(0);
  confDefault.dimToOpt.push_back(1);
  confDefault.tileDim = confDefault.dimToOpt.size();
}

ConfigurationMB::~ConfigurationMB(){

}

void ConfigurationMB::printParameters(){
  for (int i = 0; i < (int) ConfigVector.size(); i++) {
    if (ConfigVector.at(i).inteType == InteType::SoftEllipse) {
      std::cout << "inteType : SoftEllipse" << std::endl;
    }else if (ConfigVector.at(i).inteType == InteType::Heaviside) {
      std::cout << "inteType : Heaviside" << std::endl;
    }
    std::cout << "nDims :" << nDims << std::endl;
    std::cout << "base :" << base << std::endl;
    std::cout << "m :" << m << std::endl;
    std::cout << "previous :" << ConfigVector.at(i).previous << std::endl;
    std::cout << "batch :" << ConfigVector.at(i).batch << std::endl;
    std::cout << "tileDim :" << ConfigVector.at(i).tileDim << std::endl;
    for (auto i : ConfigVector.at(i).dimToOpt){
      std::cout << i << " ";
    }
    std::cout << "\n" << '\n';
  }

}


void ConfigurationMB::initByFile(std::string configFile){
  std::ifstream in(configFile,std::ios::in);
  std::string mode;
  std::string type;
  std::string line;
  bool toPush;
  while (getline(in, line)) {
      OneConfig tempConf;
      toPush = false;
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
          if (left_part == "s") {
                sline >> nDims;
          } else if (left_part == "m") {
              sline >> m;
          } else if (left_part == "b" || left_part == "p") {
              sline >> base;
          }
      } else {
          sline = std::istringstream(line);
          sline >> type;
          bool loop = true;
          do {
              if (type == "from") {
                  sline >> type;
              } else if (type == "to") {
                  sline >> type;
              } else if (type == "weak") {
                  sline >> type;
              } else {
                  loop = false;
              }
          } while (loop);
           if (type == "optimize") {
              toPush = true;
              char c;
              std::string nbInteString;
              sline >> c;
              if (c == 'p') {
                tempConf.previous = true;
              }else if (c == 'c'){
                tempConf.previous = false;
              } else {
                sline.putback(c);
              }
              sline >> c;
              if (c == 'g') {
                tempConf.inteType = SoftEllipse;
                sline >> nbInteString;
                tempConf.batch = stoi(nbInteString);
              } else if (c == 'h') {
                tempConf.inteType = Heaviside;
                sline >> nbInteString;
                tempConf.batch = stoi(nbInteString);
              } else {
                std::cerr << "Error Parsing: Unknown integrand type " << c << std::endl;
                sline.putback(c);

              }
          }
          int d;

          while (sline >> d) {
              tempConf.dimToOpt.push_back(d);
          }
          tempConf.tileDim = tempConf.dimToOpt.size();
      }
      if (toPush) {
        ConfigVector.push_back(tempConf);
      }
  }
  in.close();
}
