#include <iostream>
#include <vector>
#include <fstream>
#include "Tiles.hpp"

Tiles::Tiles(){
  samplingPoint = PointTernaire();
  afterSamplingPoint = "";
}

Tiles::Tiles(PointTernaire samplingPointParam, std::string afterSamplingPointParam){
  samplingPoint = samplingPointParam;
  afterSamplingPoint = afterSamplingPointParam;
}

Tiles::~Tiles(){

}

std::vector<Tiles>* importTiles(std::string fileToRead,int octave){
  std::ifstream tilesFile(fileToRead,std::ios::in);
  std::vector<Tiles>* vectorOfTiles = new std::vector<Tiles>;
  if (tilesFile.is_open()) {
    std::string line;
    size_t pos = 0;
    std::string token;
    int indexLocal = 0;
    int xcoLocal = 0;
    int ycoLocal = 0;
    std::string afterSamplingPointLocal;
    while (getline(tilesFile,line)) {
      pos = 0;
      token = "";
      afterSamplingPointLocal = "";

      pos = line.find("\t");
      indexLocal = std::stoi(line.substr(0,pos));
      line.erase(0, pos + 1);

      pos = line.find("\t");
      xcoLocal = std::stoi(line.substr(0,pos));
      line.erase(0, pos + 1);

      pos = line.find("\t");
      ycoLocal = std::stoi(line.substr(0,pos));
      line.erase(0, pos + 1);

      afterSamplingPointLocal = line.substr(0,std::string::npos);
      line.erase(0, std::string::npos);

      vectorOfTiles->push_back(Tiles(PointTernaire(indexLocal,xcoLocal,ycoLocal,octave),afterSamplingPointLocal));
    }
  } else {
    std::cerr << "Unable to open file" << '\n';
  }
  return vectorOfTiles;
}

void exportTiles(std::vector<Tiles>* vectorOfTiles,std::string outputString){
  std::ofstream o;
  o.open(outputString);
  for (std::vector<Tiles>::iterator it = vectorOfTiles->begin();it != vectorOfTiles->end();  it++) {
    o << (*it) << '\n';
  }
  o.close();
}

std::vector<PointTernaire>* extractSamplingPoints(std::vector<Tiles>* vectorOfTiles){
  std::vector<PointTernaire>* vectorOfPoints = new std::vector<PointTernaire>;
  vectorOfPoints->reserve(vectorOfTiles->size());
  for (std::vector<Tiles>::iterator it = vectorOfTiles->begin();it != vectorOfTiles->end();  it++) {
        vectorOfPoints->push_back((*it).samplingPoint);
    };
  return vectorOfPoints;
}
