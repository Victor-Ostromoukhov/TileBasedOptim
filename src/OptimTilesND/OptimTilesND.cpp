#include <iostream>
#include <random>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include "CLI11.hpp"
#include <chrono>
#include "Lib/Tiles/Tiles.hpp"
#ifndef DIM
  #define  DIM 2
#endif
template<int dimension>
struct newPointHolder {
  Points<dimension> point;
  double apportOfNewPoint;
} ;

template <int dimension>
bool compareTwoNewPointHolder(const newPointHolder<dimension> &a,const newPointHolder<dimension> &b){
  if (a.apportOfNewPoint < b.apportOfNewPoint) {
    return true;
  }
  return false;
}

template<int dimension>
double evaluateFirstPartPointContribution(std::vector<Points<dimension>>* pointSet,int oldPointIndex){
  double sum = 0.0;
  double produit = 1.0;
  int nb_pts = pointSet->size();
  for (int i = 1; i <= nb_pts; i++) {
    produit = 1.0;
    for (int k = 1; k <= dimension; k++) {
      produit *= (1 - (std::max(pointSet->at(i-1).dimensionnalArray[k-1],pointSet->at(oldPointIndex).dimensionnalArray[k-1]) )  );
    }
    sum += produit;
  }
  for (int i = 1; i <= nb_pts; i++) {
    produit = 1.0;
    if (i-1 != oldPointIndex) {
      for (int k = 1; k <= dimension; k++) {
          produit *= ( 1 - (std::max(pointSet->at(i-1).dimensionnalArray[k-1],pointSet->at(oldPointIndex).dimensionnalArray[k-1])));
      }
      sum+= produit;
    }
  }
  return ((1.0/double ((pow(nb_pts,2)))) * sum);
}

template<int dimension>
double evaluateSecondPartPointContribution(std::vector<Points<dimension>>* pointSet,int oldPointIndex){
  double produit = 1.0;
  int nb_pts = pointSet->size();
  for (int k = 1; k <= dimension; k++) {
    produit *= ( 1  -  pow((pointSet->at(oldPointIndex).dimensionnalArray[k-1]),2));
  }
  return  (((pow(2,-dimension+1))/double(nb_pts))*produit);
}

template<int dimension>
double evaluateL2DiscrepancyPointContribution(std::vector<Points<dimension>>* pointSet,int oldPointIndex){
  double firstPartPointContribution = evaluateFirstPartPointContribution(pointSet,oldPointIndex);
  double secondPartPointContribution = evaluateSecondPartPointContribution(pointSet,oldPointIndex);
  double thirdPartPointContribution = pow(  3  , -dimension);
  return sqrt(firstPartPointContribution - secondPartPointContribution + thirdPartPointContribution);
}

template<int dimension>
double evaluateFirstPartPointSet(std::vector<Points<dimension>>* pointSet){
  double firstSum = 0.0;
  double secondSum = 0.0;
  double product = 0.0;
  int nb_pts = pointSet->size();
  for (int i = 1; i <= nb_pts; i++) {
    secondSum = 0.0;
    for (int j = 1; j <= nb_pts; j++) {
      product = 1.0;
      for (int k = 1; k <= dimension; k++) {
        product *= (1 - (std::max(pointSet->at(i-1).dimensionnalArray[k-1],pointSet->at(j-1).dimensionnalArray[k-1]) ) );
      }
      secondSum += product;
    }
    firstSum += secondSum;
  }
  return ((1.0/double ((pow(nb_pts,2)))) * firstSum);
}

template<int dimension>
double evaluateSecondPartPointSet(std::vector<Points<dimension>>* pointSet){
  double sum = 0.0;
  double product = 1.0;
  int nb_pts = pointSet->size();
  for (int i = 1; i <= nb_pts; i++) {
    product = 1.0;
    for (int k = 1; k <= dimension; k++) {
        product *= (1 -  pow((pointSet->at(i-1).dimensionnalArray[k-1]),2));
    }
    sum += product;
  }
  return  (   (  (  pow(2, -dimension+1 )  )   / nb_pts ) * sum  );
}

template<int dimension>
double evaluateL2DiscrepancyPointSet(std::vector<Points<dimension>>* pointSet){
  double firstPartPointSet = evaluateFirstPartPointSet(pointSet);
  double secondPartPointSet = evaluateSecondPartPointSet(pointSet);
  double thirdPartPointSet = pow(  3  , -dimension);
  return sqrt(firstPartPointSet - secondPartPointSet + thirdPartPointSet);
}

template<int dimension>
std::vector<Points<dimension>> extractSP(std::vector<Tiles<dimension>>* v){
  std::vector<Points<dimension>> v2;
  for (typename std::vector<Tiles<dimension>>::iterator it = v->begin();it != v->end();  it++) {
      v2.push_back(*(*it).getSamplingPoint());
  };
  return v2;
}

template<int dimension>
void injectSP(std::vector<Tiles<dimension>>* v,std::vector<Points<dimension>>* vtoinj){
  typename std::vector<Points<dimension>>::iterator itP = vtoinj->begin();
  for (typename std::vector<Tiles<dimension>>::iterator itT = v->begin();itT != v->end();  itT++) {
    (*itT).setSamplingPoint(*itP);
    if (itP != vtoinj->end()) {
      itP++;
    }
  }
}

template<int dimension>
void exportTiles(std::vector<Tiles<dimension>>* v,std::string outputString){
  std::ofstream o;
  o.open(outputString);
  for (typename std::vector<Tiles<dimension>>::iterator it = v->begin();it != v->end();  it++) {
    // std::cout << (*it) << '\n';
    o << (*it);
  }
  o.close();
}

template<int dimension>
void exportPoints(std::vector<Points<dimension>>* v,std::string outputString){
  std::ofstream o;
  o.open(outputString);
  for (typename std::vector<Points<dimension>>::iterator it = v->begin();it != v->end();  it++) {
    // std::cout << (*it) << '\n';
    o << (*it) << '\n';
  }
  o.close();
}

std::vector<int> randomAccessMatriceGenerator(int nbpts){
  std::vector<int> v;
  for (int i = 0; i < nbpts; i++) {
    v.push_back(i);
  }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(v.begin(), v.end(),g);
  return v;
}

template<int dimension> // c plus robssute car en le sortant DIM existe plus
void optimTilesParallel(std::vector<Tiles<dimension>>* v,int nbThrow,size_t niters,size_t writeEachNIterations,std::string outputString,bool debug){
  double previousPointContribution;

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
/* la je peux laisser un vector */
  std::vector<Points<dimension>> pointSetToOptimize = extractSP(v);

  std::vector<int> rAM = randomAccessMatriceGenerator(pointSetToOptimize.size());

  newPointHolder<dimension> discTab[nbThrow];


  double initialDiscrepency = evaluateL2DiscrepancyPointSet(&pointSetToOptimize);
  double currentDisc = initialDiscrepency;

  for (size_t iter_over_pointset = 0; iter_over_pointset < niters; iter_over_pointset++) {
        rAM = randomAccessMatriceGenerator(pointSetToOptimize.size());
        for (long unsigned int i_tile = 0; i_tile < pointSetToOptimize.size(); i_tile++) {
            previousPointContribution = evaluateL2DiscrepancyPointContribution(&pointSetToOptimize,rAM.at(i_tile));
            #pragma omp parallel for private(pointSetToOptimize)
            for (int i_pt_in_tile = 0; i_pt_in_tile < nbThrow; i_pt_in_tile++) {
                pointSetToOptimize = extractSP(v);
                generator.seed((i_pt_in_tile*1234+5678)+std::chrono::system_clock::now().time_since_epoch().count());
                for (int currentDim = 1; currentDim <= dimension; currentDim++) {
                    double rand = distribution(generator);
                    v->at(rAM.at(i_tile)).vectorPointSumRand(&(discTab[i_pt_in_tile].point),*(v->at(rAM.at(i_tile)).getPreviousRefPoint()),rand,iter_over_pointset,(v->at(rAM.at(i_tile)).getPreviousDirectionnalVector(currentDim)));
                    v->at(rAM.at(i_tile)).vectorPointSumRand(&(pointSetToOptimize.at(rAM.at(i_tile))),*(v->at(rAM.at(i_tile)).getPreviousRefPoint()),rand,iter_over_pointset,v->at(rAM.at(i_tile)).getPreviousDirectionnalVector(currentDim));
                }
                discTab[i_pt_in_tile].apportOfNewPoint = evaluateL2DiscrepancyPointContribution(&pointSetToOptimize,rAM.at(i_tile));
            }
            newPointHolder<dimension> theChosenOne = *std::min_element(discTab+0,discTab+nbThrow,compareTwoNewPointHolder<dimension>);

            double newDisc = sqrt(pow(currentDisc,2) - pow(previousPointContribution,2) + pow(theChosenOne.apportOfNewPoint,2));
            if (currentDisc > newDisc ) {
              std::cout <<" Iteration " << iter_over_pointset << " : " << "Initial discrepancy : " << initialDiscrepency << ",  current discrepancy : "<< currentDisc << " becoming " << newDisc << std::endl;
              currentDisc = newDisc;
              for (int currentDim = 1; currentDim <= pointSetToOptimize.at(0).getDim(); currentDim++) {
                pointSetToOptimize.at(rAM.at(i_tile)).dimensionnalArray[currentDim-1] = theChosenOne.point.dimensionnalArray[currentDim-1];

              }
              injectSP(v,&pointSetToOptimize);
            }
        }
        if (writeEachNIterations > 0) {
          if (iter_over_pointset % writeEachNIterations == 0) {
            if (debug) {
                std::cout << "[DEBUG] Iter  " << iter_over_pointset << " : " << " exporting into " << outputString << std::endl;
                exportPoints(&pointSetToOptimize,outputString);
            }else{
                std::cout << iter_over_pointset << " : " << currentDisc << " exporting into " << outputString << std::endl;
                exportTiles(v,outputString);
            }
         }
        }
  }
  if (debug) {
    std::cout << "Exporting Pointset with a computed discrepancy of " << currentDisc  << " and real of " << evaluateL2DiscrepancyPointSet(&pointSetToOptimize)<< '\n';
    exportPoints(&pointSetToOptimize,outputString);
  }
}

template<int dimension>
int submain(int argc, char const *argv[]) {
    // =========== Variables =========== //
  srand (time(NULL));
  int nbIterationsPerTile = 64;
  // int dimension = 2;
  int nbThreads = omp_get_max_threads() == 64 ? 64 : omp_get_max_threads();
  size_t niters = 1024*1024;
  size_t writeEachNIterations = 1024;
  bool debug = false;
  std::string inputString ="pts.dat";
  std::string outputString ="OptimizedPts.dat";

  // =========== CLI11 Configuration =========== //

  CLI::App app { "OptimTilesND" };

  app.add_option("-t,--nbThreads",nbThreads,"Number of threads used , default: "+std::to_string(nbThreads) );
  app.add_option("-n,--iterationNumber",niters,"Number of iterations over the pointset, default: "+std::to_string(niters))->check(CLI::PositiveNumber);
  app.add_option("-i,--input",inputString,"Path to input file, default: "+inputString)->check(CLI::ExistingFile)->required();
  app.add_option("-o,--output",outputString,"Path to output file, default: "+outputString);
  app.add_option("-w,--writeEachNIterations",writeEachNIterations,"Will output the result at each N iterations, default: "+std::to_string(writeEachNIterations))->check(CLI::Range(100,100000));
  app.add_option("-d,--debug",debug,"Will print debug logs and output only Points within the file");
  // app.add_option("-D,--dimension",dimension,"Dimension of the tiles from the input, default: "+std::to_string(dimension))->required();

  CLI11_PARSE(app, argc, argv)

  // =========== OpenMP Configuration =========== //

    omp_set_dynamic(0);
    omp_set_num_threads(nbThreads);

    // =========== Tiles Import =========== //

    std::vector<Tiles<dimension>>* v = importTiles<dimension>(inputString);

      // =========== Tiles Optimisation =========== //

    optimTilesParallel(v,nbIterationsPerTile,niters,writeEachNIterations,outputString,debug);

      // =========== Output Writing =========== //

    if (!debug) {
      exportTiles(v,outputString);
    }

  return 0;
}

int main(int argc, char const *argv[]) {
  std::cout <<  "Dimension ==>  " <<DIM << '\n';
  return submain<DIM>(argc,argv);
}
