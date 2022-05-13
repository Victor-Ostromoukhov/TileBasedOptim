#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstring>
#include <random>
#include <omp.h>
#include <chrono>
#include "Lib/Tiles/Tiles.hpp"
#include "Lib/MultivariateGaussian/multivariateGaussian.h"
#include "Lib/MultivariateGaussian/Integration.h"
#include "Data/Intégrandes/Ellipses/HardEllipses/Ellipses2D_nIntegrands16384_optimSet.cpp"
#include "CLI11.hpp"

/* ----------- Déclaration des constantes ----------- */

#ifndef DIM
  #define DIM 2
#endif

#ifndef NBGAUSS
  #define NBGAUSS 16384
#endif
/* ----------- Déclaration des structures et leur "méthodes" ----------- */

template<int dimension>
struct newPointHolder {
  int index;
  VecX<dimension> point;
  double apportOfNewPoint;
} ;

template<int dimension>
bool compareTwoNewPointHolder(const newPointHolder<dimension> &a,const newPointHolder<dimension> &b){
  if (a.apportOfNewPoint < b.apportOfNewPoint) {
    return true;
  }
  return false;
}

/* ----------- Fonctions import ----------- */
template<int dimension>
void importPoints(double pointsTab[][dimension], std::string inputString){
  std::ifstream inputPointsFile(inputString, std::ios::in);
  if (inputPointsFile.is_open()) {
    std::string line;
    int i =0;
    size_t pos = 0;
    std::string token;
    while (getline(inputPointsFile,line)) {
        pos = 0;
        token = "";
        pos = line.find("\t");
        pointsTab[i][0] = std::stod(line.substr(0,pos));
        line.erase(0, pos + 1);
        pointsTab[i][1] = std::stod(line.substr(0,std::string::npos));
        line.erase(0, std::string::npos);
        i++;
      }
      inputPointsFile.close();
    }else{
      std::cerr << "Unable to open file\n";
    }
} // Ne sert pas pour la version tuile

template<int dimension>
std::vector<VecX<DIM>>* initializePointsVector(double pointsTab[][dimension],int nbpts){
  std::vector<VecX<DIM>>* points = new std::vector<VecX<DIM>>;
  points->resize(nbpts);
  for (int i = 0; i < nbpts; ++i) {
      for (int d = 0; d < dimension; ++d) {
          (*points)[i][d] = pointsTab[i][d];
      }
  }
  return points;
} // Ne sert pas pour la version tuile

template<int dimension>
std::vector<VecX<DIM>>* extractSP(std::vector<Tiles<dimension>>* v){
  std::vector<VecX<DIM>>* v2 = new  std::vector<VecX<DIM>>;
  v2->resize(v->size());
  v2->clear();
  for (typename std::vector<Tiles<dimension>>::iterator it = v->begin();it != v->end();  it++) {
      v2->push_back(*(*it).getSamplingPoint());
  };
  return v2;
}

       /* importTiles dans le Tiles.hpp */

/* ----------- Fonctions export ----------- */

void exportValue(int nbpts,double value,std::string outputString){
  std::ofstream o;
  o.open(outputString,std::ios_base::app);
  o << nbpts << '\t' << value << '\n';
  o.close();
} // Exports the number of points and the mse of it on a file

void exportPoints(std::vector<VecX<DIM>>* v,std::string outputString){
  std::ofstream o;
  o.open(outputString);
  for (std::vector<VecX<DIM>>::iterator it = v->begin();it != v->end();  it++) {
    o << (*it)[0] << '\t' << (*it)[1] << '\n';
  }
  o.close();
} // Exports the points in a file, not useful for tile version

template<int dimension>
void injectSP(std::vector<Tiles<dimension>>* v,std::vector<VecX<dimension>>* vtoinj){
  typename std::vector<VecX<dimension>>::iterator itP = vtoinj->begin();
  for (typename std::vector<Tiles<dimension>>::iterator itT = v->begin();itT != v->end();  itT++) {
    (*itT).setSamplingPoint(*itP);
    if (itP != vtoinj->end()) {
      itP++;
    }
  }
} // Reinjects the samplingPoints in their tiles

template<int dimension>
void exportTiles(std::vector<Tiles<dimension>>* v,std::string outputString){
  std::ofstream o;
  o.open(outputString);
  for (typename std::vector<Tiles<dimension>>::iterator it = v->begin();it != v->end();  it++) {
    o << (*it);
  }
  o.close();
}

/* ----------- Fonctions Calcul MSE ----------- */

template<int dimension>
void initializeGaussianVectors(std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,std::vector<double>* anal,int offset,t_GaussianStruct2D gaussTab[],int gaussianSubSetSize){
  for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; ++z) {

    for (int i = 0; i < dimension; ++i) {
        (*shift)[z][i] = gaussTab[z].mu[i];  // Déplacement
    }
    for (int j = 0; j < dimension; ++j) {
        for (int i = 0; i < dimension; ++i) {
            (*sigma)[z](i, j) = gaussTab[z].mxCInv[i + j * dimension]; // Matrice SR
        }
    }
    (*anal)[z] = gaussTab[z].integral; // Valeur analytique
  }
} // Fulfill the gaussian array with the parameters

template<int dimension>
double recalculateGaussianValue(double oldVal, VecX<dimension> oldPoint,VecX<dimension> newPoint, MatXDynamic sigma,VecXDynamic shift,int nbpts){
  return multivariateGaussianIntegrationPointModif(oldPoint,newPoint,shift,sigma,oldVal,nbpts);
} // Recalculate the value of the integration with a point shift

template<int dimension>
void changeAllValueGaussTab(VecX<dimension> oldPoint,VecX<dimension> newPoint, std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,double tabPtsValGauss[],int nbGauss,int nbpts){
  for (int i = 0; i < nbGauss; i++) {
    tabPtsValGauss[i] = recalculateGaussianValue(tabPtsValGauss[i],oldPoint,newPoint,(*sigma)[i],(*shift)[i],nbpts);
    }
} // Recalculate the value of the integration with a point shift for all gaussian

template<int dimension>
double recalculateGaussianValueAllGauss(VecX<dimension> oldPoint,VecX<dimension> newPoint, std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,double tabPtsValGauss[],int nbGauss,int nbpts,std::vector<double>* anal){
  double val = 0.;
  for (int i = 0; i < nbGauss; i++) {
	  double tmp = (recalculateGaussianValue(tabPtsValGauss[i],oldPoint,newPoint,(*sigma)[i],(*shift)[i],nbpts) - (*anal)[i]);

	  val+= tmp*tmp;
    }
  val /= nbGauss;
  return val;
}

template<int dimension>
double seOfAPointsetOnOneGaussian(std::vector<VecX<dimension>>* points,MatXDynamic sigma,VecXDynamic shift,double anal,int nbGauss,double tabPtsValGauss[],int indice){
    double val = 0.;
    tabPtsValGauss[indice] = multivariateGaussianIntegration((*points), shift, sigma);
    val = pow(( anal  -tabPtsValGauss[indice]),2);
    return val;
}

template<int dimension>
double mseOfAPointsetOnAllGaussian(std::vector<VecX<dimension>>* points,std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,std::vector<double>* anal,int nbGauss,double tabPtsValGauss[]){
  double val = 0.;
  for (int z = 0; z < nbGauss; ++z) {
    val += seOfAPointsetOnOneGaussian(points,(*sigma)[z],(*shift)[z],(*anal)[z],nbGauss,tabPtsValGauss,z);
  }
  val /=  nbGauss;
  tabPtsValGauss[nbGauss] = val;

  return val;
}

/* ----------- Fonctions autre ----------- */

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

/* ----------- Fonction Principale ----------- */
template<int dimension>
double optimPointME(std::vector<Tiles<DIM>>* v,int nbpts,std::string inputString,int niters,int nbThrow,std::string outputString,int gaussianSubSetSize){

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  // =========== Initialisation des variables post-command =========== //

  std::vector<MatXDynamic> sigma(NBGAUSS, MatXDynamic(dimension, dimension)); // Initialisation des matrices SR

  std::vector<VecXDynamic> shift(NBGAUSS, VecXDynamic(dimension)); // Initialisation des vecteurs de déplacements

  std::vector<double> anal(NBGAUSS, 0.); // Initialisation des valeurs calculées (analytiques) des gaussiennes

  // =========== Remplissage des variables précédentes  =========== //

  initializeGaussianVectors<dimension>(&sigma,&shift,&anal,0,tab_Ellipses2D,gaussianSubSetSize);

  // =========== Fin Remplissage des variables précédentes  =========== //

  // =========== Lecture des points  =========== //

    double pointsTab[nbpts][dimension];
    importPoints(pointsTab,inputString);

    // =========== Fin lecture  =========== //

    // =========== Début Remplissage du vecteur après lecture =========== //

    std::vector<VecX<DIM>> points = (*extractSP(v));

    // =========== Fin Remplissage du vecteur après lecture  =========== //

    // =========== Calcul de l'erreur d'intégration sur les 1024 gaussiennes pour un pointset  =========== //
    double tabPtsValGauss[gaussianSubSetSize+1];
    double initialSE = mseOfAPointsetOnAllGaussian(&points,&sigma,&shift,&anal,gaussianSubSetSize,tabPtsValGauss);
    // std::cout << initialSE << '\n';
    newPointHolder<dimension> mseTab[nbThrow];

    // =========== Initialisation matrice pour parcours aléatoire des points et des sous-sets d gaussiennes  =========== //

    std::vector<int> rAM = randomAccessMatriceGenerator(nbpts);

    std::vector<int> rAMGaussiennes;



    for (int  iter_over_pointset = 0;  iter_over_pointset < niters;  iter_over_pointset++) {
      rAM = randomAccessMatriceGenerator(nbpts);
      if (iter_over_pointset % (( NBGAUSS / gaussianSubSetSize ) - 1) == 0) {
        rAMGaussiennes = randomAccessMatriceGenerator(( NBGAUSS / gaussianSubSetSize ) - 1);
      }
      initializeGaussianVectors<dimension>(&sigma,&shift,&anal,(iter_over_pointset % (( NBGAUSS / gaussianSubSetSize ) - 1)),tab_Ellipses2D,gaussianSubSetSize);
      for (int i_pts = 0; i_pts < nbpts; i_pts++) {
          #pragma omp parallel for
          for (int i_pt_in_tile = 0; i_pt_in_tile < nbThrow; i_pt_in_tile++) {
            generator.seed((i_pt_in_tile*1234+5678)+std::chrono::system_clock::now().time_since_epoch().count());
            mseTab[i_pt_in_tile].index = rAM.at(i_pts);
            for (int currentDim = 1; currentDim <= dimension; currentDim++) {
                  mseTab[i_pt_in_tile].point.coefs[currentDim-1] = points[rAM.at(i_pts)][currentDim-1];
            }
            for (int currentDim = 1; currentDim <= dimension; currentDim++) {
                  double rand = distribution(generator);
                  v->at(rAM.at(i_pts)).vectorPointSumRand(&(mseTab[i_pt_in_tile].point),*(v->at(rAM.at(i_pts)).getPreviousRefPoint()),rand,iter_over_pointset,(v->at(rAM.at(i_pts)).getPreviousDirectionnalVector(currentDim)));
            }
            mseTab[i_pt_in_tile].apportOfNewPoint = recalculateGaussianValueAllGauss(points[rAM.at(i_pts)],mseTab[i_pt_in_tile].point, &sigma,&shift,tabPtsValGauss,gaussianSubSetSize,nbpts,&anal);
          }
          newPointHolder<dimension> theChosenOne = *std::min_element(mseTab+0,mseTab+nbThrow,compareTwoNewPointHolder<dimension>);
          std::cout << "Voila ce que j'avais " <<  tabPtsValGauss[gaussianSubSetSize] << "Et voilà ce qu'on me propose" << theChosenOne.apportOfNewPoint << '\n';
          if (theChosenOne.apportOfNewPoint < tabPtsValGauss[gaussianSubSetSize]) {
            std::cout <<" Iteration " << iter_over_pointset << " : " << "Initial SE : " << initialSE << ",  current SE : "<< tabPtsValGauss[gaussianSubSetSize] << " becoming " << theChosenOne.apportOfNewPoint << std::endl;
            tabPtsValGauss[gaussianSubSetSize] = theChosenOne.apportOfNewPoint;
            changeAllValueGaussTab(points[theChosenOne.index],theChosenOne.point, &sigma,&shift, tabPtsValGauss,gaussianSubSetSize,nbpts);
            points[theChosenOne.index][0] = theChosenOne.point[0];
            points[theChosenOne.index][1] = theChosenOne.point[1];
            injectSP(v,&points);
          }
      }
    }

    exportTiles(v,outputString);
    // exportPoints(&points,outputString);
    return tabPtsValGauss[gaussianSubSetSize];
}

int main(int argc, char const *argv[]) {

  /* ----------- Début Initialisation des variables ----------- */
  int nbpts = 2;
  int nbThreads = omp_get_max_threads() >= 64 ? 64 : omp_get_max_threads();
  std::string inputString ="pts.dat";
  std::string outputString ="OptimizedPts.dat";
  size_t niters = 1024*1024;
  int gaussianSubSetSize = 1024;
  int nbThrow = 64;
  std::string outputStringMSE ="MSE.dat";
  /* ----------- Fin Initialisation des variables ----------- */

  // =========== Début CLI11 Configuration =========== //

      CLI::App app { "OptimME" };

      app.add_option("--nbPoints",nbpts,"Number of Points, default: "+nbpts)->required();
      app.add_option("-t,--nbThreads",nbThreads,"Number of threads used , default: "+std::to_string(nbThreads))->check(CLI::Range(1,omp_get_max_threads()));
      app.add_option("-n,--iterationNumber",niters,"Number of iterations over the pointset, default: "+std::to_string(niters))->check(CLI::PositiveNumber);
      app.add_option("-i,--input",inputString,"Path to input file, default: "+inputString)->check(CLI::ExistingFile)->required();
      app.add_option("-o,--output",outputString,"Path to output file, default: "+outputString);
      app.add_option("-g",gaussianSubSetSize,"Number of gaussian to iterate over in one iteration, default: "+std::to_string(gaussianSubSetSize));
      app.add_option("-m,--writeMSE",outputStringMSE,"If precised, will write the MSE of the pointset by appending it to this file, default: "+outputStringMSE);


      CLI11_PARSE(app, argc, argv)

      // =========== Fin CLI11 Configuration =========== //
                                          /*****/
      // =========== Début OpenMP Configuration =========== //

      omp_set_dynamic(0);
      omp_set_num_threads(nbThreads);

      // =========== Fin OpenMP Configuration =========== //

        std::vector<Tiles<DIM>>* v = importTiles<DIM>(inputString);
        double val = optimPointME<DIM>(v,nbpts,inputString,niters,nbThrow,outputString,gaussianSubSetSize);
        // =========== Ecriture de l'erreur associée à un pointset  =========== //
        // std::cout << val << '\n';
        if (outputStringMSE.compare("MSE.dat") != 0) {
          exportValue(nbpts,val,outputStringMSE);
        }

        // =========== Fin du main  =========== //
  return 0;
}
