#include <iostream>
#include <vector>
#include "fstream"
#include <algorithm>
#include "Lib/Tile/Tiles.h"
#include <list>
#include <string.h>
#include <cmath>
#include <random>
#include <omp.h>
#include <chrono>
#include "CLI11.hpp"

/* ==================== Struct ==================== */

struct newPointHolder {
  Points_2D point;
  double apportOfNewPoint;
} ;

bool compareTwoNewPointHolder(const newPointHolder &a,const newPointHolder &b){
  if (a.apportOfNewPoint < b.apportOfNewPoint) {
    return true;
  }
  return false;
}

/* ==================== Utils ==================== */

float mygen (int i) {
  return rand()%i;
}

/* ==================== Evaluate the contribution of a specific point ==================== */


double evaluateFirstPartPointContribution(std::vector<Points_2D>* v,int oldpointindice){
  double somme2 = 0.0;
  double produit = 1.0;
  int nb_pts = v->size();
  // std::cout << oldpointindice << '\n';
  // std::cout << nb_pts << '\n';
  for (int i = 1; i <= nb_pts; i++) {
    produit = 1.0;
      for (int k = 1; k <= 2; k++) {
        produit *= (1 - (std::max(v->at(i-1).get_pos_dim(k),v->at(oldpointindice).get_pos_dim(k)) ) );
        // std::cout << "(1 -" << (std::max(v->at(i-1).get_pos_dim(k),v->at(oldpointindice).get_pos_dim(k)) ) << ")*" << '\n';
      }
      // std::cout << " == " << produit << '\n';
      somme2 += produit;
    }

    for (int i = 1; i <= nb_pts; i++) {
      produit = 1.0;
      if (i-1 != oldpointindice) {
        for (int k = 1; k <= 2; k++) {
          produit *= (1 - (std::max(v->at(i-1).get_pos_dim(k),v->at(oldpointindice).get_pos_dim(k)) ) );
          // std::cout << "(1 -" << (std::max(v->at(i-1).get_pos_dim(k),v->at(oldpointindice).get_pos_dim(k)) ) << ")*" << '\n';
        }
        // std::cout << " == " << produit << '\n';
        somme2 += produit;
      }
      }
      // std::cout << somme2 << '\n';
       // std::cout << "Voici l'évaluation de l'apport du " << oldpointindice << "er point en 1ere partie ==>" << ((1.0/double ((pow(nb_pts,2))))  *somme2)<< '\n';

  return ((1.0/double ((pow(nb_pts,2)))) * somme2);
}

double evaluateSecondPartPointContribution(std::vector<Points_2D>* v,int oldpointindice){
  double produit = 1.0;
  int nb_pts = v->size();
        for (int k = 1; k <= 2; k++) {
            produit *= (1 -  pow((v->at(oldpointindice).get_pos_dim(k)),2));
        }
    // std::cout << "Voici l'évaluation de l'apprt du " << oldpointindice << "er point ==>" << (   (  (  pow(2, -2+1 )  )   / double(nb_pts) ) * produit  ) << '\n';
  return  (   (  (  pow(2, -2+1 )  )   / double(nb_pts) ) * produit  );
}

double computel2StarDiscPointContribution(std::vector<Points_2D>* v,int index){
      double firstPart = evaluateFirstPartPointContribution(v,index);
      double secondPart = evaluateSecondPartPointContribution(v, index) ;
      double thirdPart = 1.0/9.0;
      return sqrt(firstPart - secondPart + thirdPart);
}

/* ==================== Evaluate the discrepancy of a PointSet ==================== */

double evaluateFirstPart(std::vector<Points_2D>* v, int nb_pts,int dimension){
  double somme1 = 0.0;
  double somme2 = 0.0;
  double produit = 0.0;
  for (int i = 1; i <= nb_pts; i++) {
    somme2 = 0.0;
    for (int j = 1; j <= nb_pts; j++) {
      produit = 1.0;
        for (int k = 1; k <= dimension; k++) {
            produit *= (1 - (std::max(v->at(i-1).get_pos_dim(k),v->at(j-1).get_pos_dim(k)) ) );
        }
        somme2 += produit;
    }
    somme1 += somme2;
    }
  return ((1.0/double ((pow(nb_pts,2)))) * somme1);
}

double evaluateSecondPart(std::vector<Points_2D>* v, int nb_pts,int dimension){
  double somme = 0.0;
  double produit = 0.0;
  for (int i = 1; i <= nb_pts; i++) {
      produit = 1.0;
        for (int k = 1; k <= dimension; k++) {
            produit *= (1 -  pow((v->at(i-1).get_pos_dim(k)),2));
        }
        somme += produit;
    }
  return  (   (  (  pow(2, -dimension+1 )  )   / nb_pts ) * somme  );
}

double computel2StarDisc(std::vector<Points_2D>* v){
    int dimension = 2;
    double firstPart = evaluateFirstPart(v, v->size(),2);
    double secondPart = evaluateSecondPart(v,  v->size(),2) ;
    double thirdPart = + pow(  3  , -dimension  );
    return sqrt(firstPart - secondPart + thirdPart);
}


/* ==================== Import the Tiles from a file and store them in a Vector of Tiles ==================== */

std::vector<Tiles>* importTiles(std::string fileToRead){
  std::ifstream tilesFile(fileToRead, std::ios::in);
  std::vector<Tiles>* vectorOfTiles = new std::vector<Tiles>;
  if (tilesFile.is_open()) {
    std::string bef = "";
    Points_2D samplingPoint;
    Points_2D previousRefPoint;
    Points_2D previousv1;
    Points_2D previousv2;
    std::string afterPreviousv2 = "";
    double val1;
    double val2;
    std::string line;
    size_t pos = 0;
    std::string token;

    while (getline(tilesFile,line)) {
      bef = "";
      afterPreviousv2 = "";
      pos = 0;
      token = "";
      // ==================== Parsing of a Tile ==================== //

      // ======= beforeSamplingPoint ======= //
      // === TileType === //
      pos = line.find("\t");
      bef += line.substr(0,pos)+"\t";
      line.erase(0, pos + 1);

      // === Structural ID === //
      pos = line.find("\t");
      bef += line.substr(0,pos);
      line.erase(0, pos + 1);

      // ======= samplingPoint ======= //
      // === x value (first dimension) === //
      pos = line.find("\t");
      val1 = std::stod(line.substr(0,pos));
      line.erase(0, pos + 1);

      // === y value (second dimension) === //
      pos = line.find("\t");
      val2 = std::stod(line.substr(0,pos));
      line.erase(0, pos + 1);

      // === Point Setting === //
      samplingPoint.set_pos_x(val1);
      samplingPoint.set_pos_y(val2);

      // ======= previousRefPoint ======= //
      // === x value (first dimension) === //
      pos = line.find("\t");
      val1 = std::stod(line.substr(0,pos));
      line.erase(0, pos + 1);

      // === y value (second dimension) === //
      pos = line.find("\t");
      val2 = std::stod(line.substr(0,pos));
      line.erase(0, pos + 1);

      // === Point Setting === //
      previousRefPoint.set_pos_x(val1);
      previousRefPoint.set_pos_y(val2);

      // ======= previousv1 ======= //
      // === x value (first dimension) === //
      pos = line.find("\t");
      val1 = std::stod(line.substr(0,pos));
      line.erase(0, pos + 1);

      // === y value (second dimension) === //
      pos = line.find("\t");
      val2 = std::stod(line.substr(0,pos));
      line.erase(0, pos + 1);

      // === Point Setting === //
      previousv1.set_pos_x(val1);
      previousv1.set_pos_y(val2);

      // ======= previousv2 ======= //
      // === x value (first dimension) === //
      pos = line.find("\t");
      val1 = std::stod(line.substr(0,pos));
      line.erase(0, pos + 1);

      // === y value (second dimension) === //
      pos = line.find("\t");
      val2 = std::stod(line.substr(0,pos));
      line.erase(0, pos + 1);

      // === Point Setting === //
      previousv2.set_pos_x(val1);
      previousv2.set_pos_y(val2);

      // ======= afterPreviousv2 ======= //
      afterPreviousv2 = line.substr(0,std::string::npos);
      line.erase(0, std::string::npos);

      // ==================== Push of a Tile in the Vector ==================== //
      vectorOfTiles->push_back(Tiles(bef,samplingPoint,previousRefPoint,previousv1,previousv2,afterPreviousv2));
    }
    tilesFile.close();
  }
  else {
    std::cerr << "Unable to open file\n";
  }
  return vectorOfTiles;
}

/* ==================== Export the Tiles from the Vector and write them in a file ==================== */

void exportTiles(std::vector<Tiles>* v,std::string outputString){
  std::ofstream o;
  o.open(outputString);
  for (std::vector<Tiles>::iterator it = v->begin();it != v->end();  it++) {
    o << (*it) << '\n';
  }
  o.close();
}

/* ==================== Import the Tiles from a file and store them in a Vector of Tiles ==================== */

void exportPoints(std::vector<Points_2D>* v,std::string outputString){
  std::ofstream o;
  o.open(outputString);
  for (std::vector<Points_2D>::iterator it = v->begin();it != v->end();  it++) {
    o << (*it) << '\n';
  }
  o.close();
}

/* ==================== Extract only the samplingPoints of the Tiles Vector ==================== */

std::vector<Points_2D> extractSP(std::vector<Tiles>* v){
  std::vector<Points_2D> v2;
  for (std::vector<Tiles>::iterator it = v->begin();it != v->end();  it++) {
      v2.push_back((*it).getSamplingPoint());
  };
  return v2;
}

/* ==================== Inject the samplingPoints back into the Tiles Vector ==================== */

void injectSP(std::vector<Tiles>* v,std::vector<Points_2D>* vtoinj){
  std::vector<Points_2D>::iterator itP = vtoinj->begin();
  for (std::vector<Tiles>::iterator itT = v->begin();itT != v->end();  itT++) {
    (*itT).setSamplingPoint(*itP);
    if (itP != vtoinj->end()) {
      itP++;
    }
  }
}

/* ==================== Creates a Vector in order to go all over the Tiles in a different order each time ==================== */

std::vector<int> randomAccessMatriceGenerator(int nbpts){
  std::vector<int> v;
  for (int i = 0; i < nbpts; i++) {
    v.push_back(i);
  }
  random_shuffle(v.begin(), v.end(),mygen);
  return v;
}

/* ==================== Main Function : Optimize the discrepancy of a PointSet contained within a Tile-based pattern  ==================== */

void optimTilesParallel(std::vector<Tiles>* v,int nbThrow,size_t niters,size_t writeEachNIterations,std::string outputString,bool debug) {

  // ==================== Variable Initializer ==================== //
  Points_2D oldPoint;
  double apportdupointavant;
  double newPointx;
  double newPointy;

  // ======= Randomizer ======= //
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  // ======= PointSet ======= //
  std::vector<Points_2D> pointSetToOptimize = extractSP(v);

  // ======= Matrix to access Tiles in a different order each time ======= //
  std::vector<int> rAM = randomAccessMatriceGenerator(pointSetToOptimize.size());

  // ======= NewPointHolder array to store values of the different threads ======= //
  newPointHolder discTab[nbThrow];

  // ======= Discrepancy of the current PointSet (Here, the initial one) ======= //
  double currentDisc = computel2StarDisc(&pointSetToOptimize);

  // ==================== Iteration over the Tiles N times ==================== //

  for (size_t iter = 0; iter < niters; iter++) {
    rAM = randomAccessMatriceGenerator(pointSetToOptimize.size());

    // ============= Beginning to go all over the PointSet ================== //
    for (long unsigned int i = 0; i < pointSetToOptimize.size(); i++) {

      oldPoint = pointSetToOptimize.at(rAM.at(i));
      apportdupointavant = computel2StarDiscPointContribution(&pointSetToOptimize,rAM.at(i));
      // ============= Parallel Beginning ================== //

      #pragma omp parallel for private(pointSetToOptimize,newPointx,newPointy)
        for (int j = 0; j < nbThrow; j++) {
          pointSetToOptimize = extractSP(v);
          generator.seed((omp_get_thread_num()*1234+5678)+std::chrono::system_clock::now().time_since_epoch().count());
          newPointx = v->at(rAM.at(i)).getPreviousRefPoint().get_pos_x() + v->at(rAM.at(i)).getPreviousv1().get_pos_x() * distribution(generator) +  v->at(rAM.at(i)).getPreviousv2().get_pos_x() * distribution(generator);
          newPointy = v->at(rAM.at(i)).getPreviousRefPoint().get_pos_y()+ v->at(rAM.at(i)).getPreviousv1().get_pos_y() * distribution(generator) + v->at(rAM.at(i)).getPreviousv2().get_pos_y() * distribution(generator);
          discTab[omp_get_thread_num()].point.set_pos_x(newPointx);
          discTab[omp_get_thread_num()].point.set_pos_y(newPointy);
          pointSetToOptimize.at(rAM.at(i)).set_pos_x(newPointx);
          pointSetToOptimize.at(rAM.at(i)).set_pos_y(newPointy);
          discTab[omp_get_thread_num()].apportOfNewPoint = computel2StarDiscPointContribution(&pointSetToOptimize,rAM.at(i));
        }

        // ============= Parallel Ending ================== //

        // ============= Decide to replace the point or not to replace it if necessary, along with the discrepancy of the PointSet ================== //
        newPointHolder theChosenOne = *std::min_element(discTab+0,discTab+nbThrow,compareTwoNewPointHolder);

        double newDisc = sqrt(pow(currentDisc,2) - pow(apportdupointavant,2) + pow(theChosenOne.apportOfNewPoint,2));

        if (currentDisc > newDisc ) {
          std::cout << iter << " : " << currentDisc << " becoming " << newDisc << std::endl;
          currentDisc = newDisc;
          pointSetToOptimize.at(rAM.at(i)).set_pos_x(theChosenOne.point.get_pos_x());
          pointSetToOptimize.at(rAM.at(i)).set_pos_y(theChosenOne.point.get_pos_y());
          injectSP(v,&pointSetToOptimize);
        }

      }

      // ============= Reinject the modified points into the Tiles ================== //

      // ============= Export the file at regular interval (Number of iterations) ================== //
      if (writeEachNIterations > 0) {
        if (iter % writeEachNIterations == 0) {
          if (debug) {
            std::cout << "[DEBUG] Iter  " << iter << " : " << " exporting into " << outputString << std::endl;
            exportPoints(&pointSetToOptimize,outputString);
          }else{
            std::cout << iter << " : " << currentDisc << " exporting into " << outputString << std::endl;
            exportTiles(v,outputString);
          }
        }
      }
  }
  if (debug) {
    std::cout << "J'exporte un pointset avec une discrépance de " << currentDisc << '\n';
    exportPoints(&pointSetToOptimize,outputString);
  }
}


int main(int argc, char const *argv[]) {
  // =========== Variables =========== //
  srand (time(NULL));
  int nbIterationsPerTile = omp_get_max_threads() == 64 ? 64 : omp_get_max_threads();
  size_t niters = 1024*1024;
  size_t writeEachNIterations = 1024;
  bool debug = false;
  std::string inputString ="pts.dat";
  std::string outputString ="OptimizedPts.dat";

  // =========== CLI11 Configuration =========== //

  CLI::App app { "OptimTiles2D" };

  app.add_option("-t,--nbIterationsPerTile",nbIterationsPerTile,"Number of thread used (is also the number of throw on a tile), default: "+std::to_string(nbIterationsPerTile))->check(CLI::Range(1,nbIterationsPerTile));
  app.add_option("-n,--iterationNumber",niters,"Number of iterations over the pointset, default: "+std::to_string(niters))->check(CLI::PositiveNumber);
  app.add_option("-i,--input",inputString,"Path to input file, default: "+inputString)->check(CLI::ExistingFile)->required();
  app.add_option("-o,--output",outputString,"Path to output file, default: "+outputString);
  app.add_option("-w,--writeEachNIterations",writeEachNIterations,"Will output the result at each N iterations, default: "+std::to_string(writeEachNIterations))->check(CLI::Range(100,100000));
  app.add_option("-d,--debug",debug,"Will print debug logs and output only Points within the file");

  CLI11_PARSE(app, argc, argv)
  // =========== OpenMP Configuration =========== //
  omp_set_dynamic(0);
  omp_set_num_threads(nbIterationsPerTile);

  // =========== Tiles Import =========== //

  std::vector<Tiles>* v = importTiles(inputString);

  // =========== Tiles Optimisation =========== //

  optimTilesParallel(v,nbIterationsPerTile,niters,writeEachNIterations,outputString,debug);

  // =========== Output Writing =========== //
  if (!debug) {
    exportTiles(v,outputString);
  }


  return 0;
}
