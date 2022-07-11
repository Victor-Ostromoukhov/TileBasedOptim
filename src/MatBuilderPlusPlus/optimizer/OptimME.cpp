#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstring>
#include <random>
#include <omp.h>
#include <chrono>
#include <iomanip>
#include "Lib/Tiles/Tiles.hpp"
#include "Lib/Integrand/Integrand.hpp"
#include "../include/CLI11.hpp"
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
/* ----------- Declaration of constants  ----------- */

 /* -------- Number of dimension, meant to be changed using compiling options -------- */

#ifndef DIM
  #define DIM 2
#endif

/* -------- Integrands type -------- */

#define INTEGRAND_TYPE_HEAVISIDE 1
#define INTEGRAND_TYPE_SOFTELLIPSES 2
#define INTEGRAND_TYPE_HARDELLIPSES 3
#define INTEGRAND_TYPE_SOFTRECTANGLES 4
#define INTEGRAND_TYPE_HARDRECTANGLES 5

/* -------- Number of integrands for each integrands, for each dimnsions -------- */

  /* -------- Dimension 2 -------- */

#define INTEGRAND_TYPE_HEAVISIDE_2D_NUMBER 262144
#define INTEGRAND_TYPE_SOFTELLIPSES_2D_NUMBER 524288
#define INTEGRAND_TYPE_HARDELLIPSES_2D_NUMBER 16384
#define INTEGRAND_TYPE_SOFTRECTANGLES_2D_NUMBER 16384
#define INTEGRAND_TYPE_HARDRECTANGLES_2D_NUMBER 16384

  /* -------- Dimension 3 //TODO -------- */

#define INTEGRAND_TYPE_HEAVISIDE_3D_NUMBER 16384
#define INTEGRAND_TYPE_SOFTELLIPSES_3D_NUMBER 16384

  /* -------- Dimension 4 //TODO -------- */


/* ----------- Declaration of global variables  ----------- */

int total_N_integrands;

/* ----------- Declaration of integrands structures  ----------- */

  /* -------- Dimension 2 -------- */

  struct t_GaussianStruct2D {
      double integral;
      double mu[2];
      double mxCInv[2 * 2];
  };

  struct t_SoftRectanglesStruct2D{
    double integral;
    double mu[2];
    double sigma[2];
  };

  struct t_Heaviside2D {
    double integral;
    double muDiscotinuity[2];
    double normal[2];
  };

  struct t_RectanglesStruct2D{
    double integral;
    double mu[2];
    double sigma[2];
  };

  /* -------- Dimension 3 //TODO -------- */

  struct t_Gauss3D {
    double integral;
    double mu[3] ;
    double mxCInv[3 * 3] ;
  };

  struct t_Heaviside3D {
    double integral;
    double muDiscotinuity[3];
    double normal[3];
  };

  /* -------- Dimension 4 //TODO -------- */


/* ----------- Declaration of integrands arrays  ----------- */

  /* -------- Dimension 2 -------- */

  extern t_Heaviside2D tab_Heaviside2D[INTEGRAND_TYPE_HEAVISIDE_2D_NUMBER];
  extern t_GaussianStruct2D tab_SoftEllipses2D[INTEGRAND_TYPE_SOFTELLIPSES_2D_NUMBER];
  extern t_GaussianStruct2D tab_Ellipses2D[INTEGRAND_TYPE_HARDELLIPSES_2D_NUMBER];
  extern t_SoftRectanglesStruct2D tab_SoftRectangles2D[INTEGRAND_TYPE_SOFTRECTANGLES_2D_NUMBER];
  extern t_RectanglesStruct2D tab_Rectangles2D[INTEGRAND_TYPE_HARDRECTANGLES_2D_NUMBER];
  extern t_RectanglesStruct2D tab_Rectangles2D[INTEGRAND_TYPE_HARDRECTANGLES_2D_NUMBER];

  /* -------- Dimension 3 -------- */

  extern t_Gauss3D tab_Gauss3D[INTEGRAND_TYPE_SOFTELLIPSES_3D_NUMBER];
  extern t_Heaviside3D tab_Heaviside3D[INTEGRAND_TYPE_HEAVISIDE_3D_NUMBER];

  /* -------- Dimension 4 //TODO -------- */

/* ----------- Other structures and function related to it ----------- */

  /* -------- Structure meant to hold the points and their contribution before comparing them within the optimization process -------- */

template<int dimension>
struct newPointHolder {
  int index;
  VecX<dimension> point;
  double apportOfNewPoint;
};

  /* -------- Compares two pointHolder structures by comparing their contribution to the MSE -------- */
template<int dimension>
bool compareTwoNewPointHolder(const newPointHolder<dimension> &a,const newPointHolder<dimension> &b){
  if (a.apportOfNewPoint < b.apportOfNewPoint) {
    return true;
  }
  return false;
}

/* ----------- Import Functions ----------- */

  /* -------- Creates the tiles using the points from pointHolder -------- */
template<int dimension>
std::vector<Tiles<dimension>>* tilesOnTheFly(int nbpts,int base,bool previous,std::vector<int>* dimToOpt,int innoctave,double** pointHolder){

  std::vector<Tiles<dimension>>* vectorOfTiles = new std::vector<Tiles<dimension>>;

  switch (dimension) {
    case 2:
      {
      /* -------- Creation on the first, always the same and independant of the options specified -------- */

      Tiles<dimension> firstTile;
      firstTile.beforeSamplingPoint = "0  ";
      firstTile.samplingPoint[0] = pointHolder[0][dimToOpt->at(0)];
      firstTile.samplingPoint[1] = pointHolder[0][dimToOpt->at(1)];
      firstTile.previousRefPoint[0] = 0;
      firstTile.previousRefPoint[1] = 0;
      firstTile.previousDirectionnalVectors[0][0] = 1;
      firstTile.previousDirectionnalVectors[0][1] = 0;
      firstTile.previousDirectionnalVectors[1][0] = 0;
      firstTile.previousDirectionnalVectors[1][1] = 1;
      vectorOfTiles->push_back(firstTile);


      for (int ioctave = 1; ioctave <= innoctave; ioctave++) {
        int iOrdinalAbsoluteFrom = std::pow(base, ioctave-1)+1;
        int iOrdinalAbsoluteTo = pow(base,ioctave);
        for (int iOrdinalAbsolute = iOrdinalAbsoluteFrom; iOrdinalAbsolute <= iOrdinalAbsoluteTo; iOrdinalAbsolute ++) {

          int npts = pow(base,ioctave);
          int prevoctave = ioctave-1;
          double x = npts * pointHolder[iOrdinalAbsolute-1][dimToOpt->at(0)];
          double y = npts * pointHolder[iOrdinalAbsolute-1][dimToOpt->at(1)];
          std::vector<double> v1;
          std::vector<double> v2;
          std::vector<double> prevv1;
          std::vector<double> prevv2;
          double prevrefx;
          double prevrefy;
          double refx;
          double refy;

          if (prevoctave % 2 == 0) {
            int ix = x / (pow(base, ((ioctave+1)/2)));
            int iy = y / (pow(base, ((ioctave+1)/2)));
            int k1,k2;
            if ((int)(ix+iy) % 2 == 0 ) {
              k1 = pow(base, ((ioctave-1)/2));
              k2 = pow(base, ((ioctave+1)/2));
            } else {
              k1 = pow(base, ((ioctave+1)/2));
              k2 = pow(base, ((ioctave-1)/2));
            }
            std::vector<double> v1;
            v1.push_back(1.0/k1);
            v1.push_back(0);
            std::vector<double> v2;
            v2.push_back(0);
            v2.push_back(1.0/k2);
            refx =(double) (((int) x / k2)) / k1;
            refy =(double) (((int) y / k1)) / k2;
            prevv1.push_back((1.0)/(pow(base,((ioctave-1)/2) )));
            prevv1.push_back(0);
            prevv2.push_back(0);
            prevv2.push_back((1.0)/(pow(base,((ioctave-1)/2) )));
            prevrefx = (ix) / pow(base, ((ioctave-1)/2));
            prevrefy = (iy) / pow(base, ((ioctave-1)/2));
          }else {
            int dx = pow(base, (ioctave/2));
            int dy = pow(base, (ioctave/2));
            refx = (double) ((int) x / dx ) / dx;
            refy =(double) ((int)y / dy ) / dy;
            v1.push_back(1.0 / dx);
            v1.push_back(0.0);
            v2.push_back(0.0);
            v2.push_back(1.0 / dy);
            int ix = x / pow(base, ((ioctave+2)/2));
            int iy = y / pow(base, ((ioctave+2)/2));
            int k1,k2;
            if ((int)(ix+iy) % 2 == 0 ) {
              k1 = pow(base, (ioctave)/2);
              k2 = pow(base, ((ioctave+2)/2));
            } else {
              k1 = pow(base, ((ioctave+2)/2));
              k2 = pow(base, (ioctave)/2);
            }
            prevv1.push_back(3 * (1.0 /(double) k1));
            prevv1.push_back(0);
            prevv2.push_back(0);
            prevv2.push_back(3 * (1.0 / k2));
            prevrefx = 3 * ((double) ((int)x / k2) / k1);
            prevrefy = 3 * ((double)((int)y / k1) / k2);
          }
          if (previous) {
            Tiles<dimension> newTile;
            newTile.beforeSamplingPoint = std::to_string(iOrdinalAbsolute-1) + '\t';
            newTile.samplingPoint[0] = x/npts;
            newTile.samplingPoint[1] = y/npts;
            newTile.previousRefPoint[0] = prevrefx;
            newTile.previousRefPoint[1] = prevrefy;
            for (int i = 0; i < dimension; i++) {
              newTile.previousDirectionnalVectors[0][i] = prevv1.at(i);
            }
            for (int i = 0; i < dimension; i++) {
              newTile.previousDirectionnalVectors[1][i] = prevv2.at(i);
            }
            vectorOfTiles->push_back(newTile);
            if ((int)vectorOfTiles->size() == nbpts) {
              return vectorOfTiles;
            }
           } else {
             Tiles<dimension> newTile;
             newTile.beforeSamplingPoint = std::to_string(iOrdinalAbsolute-1) + '\t';
             newTile.samplingPoint[0] = x/npts;
             newTile.samplingPoint[1] = y/npts;
             newTile.previousRefPoint[0] = refx;
             newTile.previousRefPoint[1] = refy;
             for (int i = 0; i < dimension; i++) {
               newTile.previousDirectionnalVectors[0][i] = v1.at(i);
             }
             for (int i = 0; i < dimension; i++) {
               newTile.previousDirectionnalVectors[1][i] = v2.at(i);
             }
             vectorOfTiles->push_back(newTile);
             if ((int)vectorOfTiles->size() == nbpts) {
               return vectorOfTiles;
             }
           }
        }
      }
      break;
      }
    /* -------- TODO Implement the creation of the tiles in the other dimensions
    case 3:

      break;
    case 4;

      break;
    -------- */
    default:
      std::cerr << "The creation of the tiles in this dimension is not supported ( yet ? )." << '\n';
      exit(-1);
  }

  return vectorOfTiles;
}

  /* -------- Initilazes pointHolder by reading the file given as a parameter -------- */
template<int dimension>
void importPoints(double** pointHolder, std::string inputString,int nbTotal){
  std::ifstream inputPointsFile(inputString, std::ios::in);
  std::string line;
  if (inputPointsFile.is_open()) {
    for (int nbPts = 0; nbPts < nbTotal; nbPts++) {
      std::getline(inputPointsFile,line);
      std::istringstream sline(line);
      for (int j = 0; j < dimension; j++) {
        sline >> pointHolder[nbPts][j];
      }
    }
    inputPointsFile.close();
    }else{
      std::cerr << "Unable to open file\n";
    }
}

  /* -------- Creates a vector of sampling points by extracting it from the vector of tiles -------- */
template<int dimension>
std::vector<VecX<DIM>>* extractSP(std::vector<Tiles<dimension>>* vectorTiles){
  std::vector<VecX<DIM>>* vectorSP = new  std::vector<VecX<DIM>>;
  vectorSP->resize(vectorTiles->size());
  vectorSP->clear();
  for (typename std::vector<Tiles<dimension>>::iterator it = vectorTiles->begin();it != vectorTiles->end();  it++) {
      vectorSP->push_back(*(*it).getSamplingPoint());
  };
  return vectorSP;
}

  /* -------- Injects the sampling points contained within the tiles back into pointHolder, at the right dimension -------- */
template<int dimension>
void injectSPInPointHolder(std::vector<Tiles<dimension>>* v, double** pointHolder,int nbpts, std::vector<int>* dimToOpt){
  for (int i = 0; i < nbpts; i++) {
    for (int j = 0; j < (int)dimToOpt->size(); j++) {
      pointHolder[i][dimToOpt->at(j)] = v->at(i).samplingPoint[j];
    }
  }
}


/* ----------- Export functions ----------- */

  /* -------- Writes all of the points of pointHolder in the file specified as a parameter -------- */
template<int dimension>
void exportPoints(double** pointHolder,int nbpts, std::string outputString){
  std::ofstream out;
  out.open(outputString);
  out << std::setprecision(17);
  for (int i = 0; i < nbpts; i++) {
    for (int j = 0; j < dimension; j++) {
      out << pointHolder[i][j] << '\t';
    }
    out << std::endl;
  }
  out.close();
}

  /* -------- Injects the sampling points back into their tiles -------- */
template<int dimension>
void injectSP(std::vector<Tiles<dimension>>* vectorTiles,std::vector<VecX<dimension>>* vectorToInject){
  typename std::vector<VecX<dimension>>::iterator iteratorToInject = vectorToInject->begin();
  for (typename std::vector<Tiles<dimension>>::iterator iteratorTiles = vectorTiles->begin(); iteratorTiles != vectorTiles->end();  iteratorTiles++) {
    (*iteratorTiles).setSamplingPoint(*iteratorToInject);
    if (iteratorToInject != vectorToInject->end()) {
      iteratorToInject++;
    }
  }
}

/* ----------- Functions for the Mean Squared Error ----------- */

  /* -------- Initilazes the batch of integrands, depending on the size of the batch and the integrand chosen -------- */
template<int dimension>
void initializeGaussianVectors(std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,std::vector<double>* anal,int offset,int integrandType,int gaussianSubSetSize){
  switch (integrandType) {
    case 1: // HeaviSide
      for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
        for (int i = 0; i < dimension; ++i) {
            (*shift)[z-(offset*gaussianSubSetSize)][i] = tab_Heaviside2D[z].muDiscotinuity[i];  // Déplacement
        }
        for (int j = 0; j < dimension; ++j) {
            (*sigma)[z-(offset*gaussianSubSetSize)](0, j) = tab_Heaviside2D[z].normal[j]; // Matrice SR
        }
        (*anal)[z-(offset*gaussianSubSetSize)] = tab_Heaviside2D[z].integral; // Valeur analytique
      }
    break;
    case 2: // SoftEllipses
        for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {

          for (int i = 0; i < dimension; ++i) {
              (*shift)[z-(offset*gaussianSubSetSize)][i] = tab_SoftEllipses2D[z].mu[i];  // Déplacement
          }
          for (int j = 0; j < dimension; ++j) {
              for (int i = 0; i < dimension; ++i) {
                  (*sigma)[z-(offset*gaussianSubSetSize)](i, j) = tab_SoftEllipses2D[z].mxCInv[i + j * dimension]; // Matrice SR
              }
          }
          (*anal)[z-(offset*gaussianSubSetSize)] = tab_SoftEllipses2D[z].integral; // Valeur analytique
        }
    break;
    case 3: // HardEllipses

      for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
        for (int i = 0; i < dimension; ++i) {
            (*shift)[z-(offset*gaussianSubSetSize)][i] = tab_Ellipses2D[z].mu[i];  // Déplacement
        }
        for (int j = 0; j < dimension; ++j) {
            for (int i = 0; i < dimension; ++i) {
                (*sigma)[z-(offset*gaussianSubSetSize)](i, j) = tab_Ellipses2D[z].mxCInv[i + j * dimension]; // Matrice SR
            }
        }
        (*anal)[z-(offset*gaussianSubSetSize)] = tab_Ellipses2D[z].integral; // Valeur analytique
      }
    break;
    case 4: // SoftRectangles
        for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
          for (int i = 0; i < dimension; ++i) {
              (*shift)[z-(offset*gaussianSubSetSize)][i] = tab_SoftRectangles2D[z].mu[i];  // Déplacement
          }
          for (int j = 0; j < dimension; ++j) {
              (*sigma)[z-(offset*gaussianSubSetSize)](0, j) = tab_SoftRectangles2D[z].sigma[j]; // Matrice SR
          }
          (*anal)[z-(offset*gaussianSubSetSize)] = tab_SoftRectangles2D[z].integral; // Valeur analytique
        }
    break;
    case 5: // HardRectangles
        for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
          for (int i = 0; i < dimension; ++i) {
              (*shift)[z-(offset*gaussianSubSetSize)][i] = tab_Rectangles2D[z].mu[i];  // Déplacement
          }
          for (int j = 0; j < dimension; ++j) {
              (*sigma)[z-(offset*gaussianSubSetSize)](0, j) = tab_Rectangles2D[z].sigma[j]; // Matrice SR
          }
          (*anal)[z-(offset*gaussianSubSetSize)] = tab_Rectangles2D[z].integral; // Valeur analytique
        }
    break;
    default:
      std::cerr << "The chosen integrand type hasn't been implemented yet, so stay tuned." << '\n';
    break;
  }
}

 /* -------- Recalculate the value of the integration when only one point has been moved -------- */
template<int dimension>
double recalculateGaussianValue(double oldVal, VecX<dimension> oldPoint,VecX<dimension> newPoint, MatXDynamic sigma,VecXDynamic shift,int nbpts,int integrandType){
  return multivariateGaussianIntegrationPointModif(oldPoint,newPoint,shift,sigma,oldVal,nbpts,integrandType);
}

  /* -------- Recalculate the value of the integration when only one point has been moved for all the integrands and stores it in the array  -------- */
template<int dimension>
void changeAllValueGaussTab(VecX<dimension> oldPoint,VecX<dimension> newPoint, std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,double tabPtsValGauss[],int nbGauss,int nbpts,int integrandType){
  for (int i = 0; i < nbGauss; i++) {
    tabPtsValGauss[i] = recalculateGaussianValue(tabPtsValGauss[i],oldPoint,newPoint,(*sigma)[i],(*shift)[i],nbpts,integrandType);
    }
}



template<int dimension>
double recalculateGaussianValueAllGauss(VecX<dimension> oldPoint,VecX<dimension> newPoint, std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,double tabPtsValGauss[],int nbGauss,int nbpts,std::vector<double>* anal,int integrandType){
  double val = 0.;
  for (int i = 0; i < nbGauss; i++) {
    val+= pow(((*anal)[i] - recalculateGaussianValue(tabPtsValGauss[i],oldPoint,newPoint,(*sigma)[i],(*shift)[i],nbpts,integrandType)),2);
    }
  val /= nbGauss;
  return val;
}

  /* -------- Calculate the squared error of the pointset for one integrand  -------- */
template<int dimension>
double seOfAPointsetOnOneGaussian(std::vector<VecX<dimension>>* points,MatXDynamic sigma,VecXDynamic shift,double anal,int nbGauss,double tabPtsValGauss[],int indice,int integrandType){
    double val = 0.;
    tabPtsValGauss[indice] = multivariateGaussianIntegration((*points), shift, sigma,integrandType);
    val = pow(( anal - tabPtsValGauss[indice]),2);
    return val;
}

  /* -------- Calculate the Mean Squared Error of the pointset for all the integrands  -------- */
template<int dimension>
double mseOfAPointsetOnAllGaussian(std::vector<VecX<dimension>>* points,std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,std::vector<double>* anal,int nbGauss,double tabPtsValGauss[],int integrandType){
  double val = 0.;
  for (int z = 0; z < nbGauss; ++z) {
    val += seOfAPointsetOnOneGaussian(points,(*sigma)[z],(*shift)[z],(*anal)[z],nbGauss,tabPtsValGauss,z,integrandType);
  }
  val /=  nbGauss;
  tabPtsValGauss[nbGauss] = val;
  return val;
}

/* ----------- Other functions ----------- */

  /* -------- Generates a vector filled with integer from limit to nbpts then shuffles it randomly  -------- */
std::vector<int> randomAccessMatriceGenerator(int nbpts,int limit){
  std::vector<int> v;
  for (int i = limit; i < nbpts; i++) {
    v.push_back(i);
  }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(v.begin(), v.end(),g);
  return v;
}

/* ----------- Function that realizes the optimization ----------- */
template<int dimension>
void optimPointME(std::vector<Tiles<DIM>>* v,int nbpts,int niters,int nbThrow,std::string outputString,int gaussianSubSetSize,int integrandType,int limit, int intervalToWrite,int nbPtsTotal,double **pointHolder, std::vector<int>* dimToOpt){

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  /* -------- Declaration of the integrand related variables -------- */

  std::vector<MatXDynamic> sigma(gaussianSubSetSize, MatXDynamic(dimension, dimension));

  std::vector<VecXDynamic> shift(gaussianSubSetSize, VecXDynamic(dimension));

  std::vector<double> anal(gaussianSubSetSize, 0.); // Ground-truth of the integration of the integral

  /* -------- Initialization of the integrand related variables -------- */

  initializeGaussianVectors<dimension>(&sigma,&shift,&anal,0,integrandType,gaussianSubSetSize);

  /* -------- Extracting the sampling points of the tile into a vector of VecX -------- */

  std::vector<VecX<DIM>> points = (*extractSP(v));

  /* -------- Declaration of the array that will contain the SE of the pointset for one integrand and the mean of those values in its last index, i.e. its size -------- */

  double tabPtsValGauss[gaussianSubSetSize+1];

  /* -------- Initialization of the initial MSE and of the array, using the same function -------- */

  double initialSE = mseOfAPointsetOnAllGaussian(&points,&sigma,&shift,&anal,gaussianSubSetSize,tabPtsValGauss,integrandType);

  /* -------- Declaration of the array that will contain the different tries on one tile, using the strucuture declared at the beginning -------- */

  newPointHolder<dimension> mseTab[nbThrow];

  /* -------- Initialization of the vector that will be used to iterate over the points, limit being there for the sequence -------- */

  std::vector<int> rAM = randomAccessMatriceGenerator(nbpts,limit);

  /* -------- Same as rAM, but for the integrands -------- */

  std::vector<int> rAMGaussiennes;

  double prevMSE = 0.;

  /* -------- The iterations over the pointset -------- */
  for (int  iter_over_pointset = 0;  iter_over_pointset < niters;  iter_over_pointset++) {

     /* -------- At each time that we change of iteration, we change the way we iterate over the points and take a different batch of integrands -------- */

    rAM = randomAccessMatriceGenerator(nbpts,limit);
    if (iter_over_pointset % ( total_N_integrands / gaussianSubSetSize ) == 0) {
      rAMGaussiennes = randomAccessMatriceGenerator(( total_N_integrands / gaussianSubSetSize ),0);
    }
    initializeGaussianVectors<dimension>(&sigma,&shift,&anal,rAMGaussiennes.at((iter_over_pointset % (( total_N_integrands / gaussianSubSetSize )))),integrandType,gaussianSubSetSize);
    mseOfAPointsetOnAllGaussian(&points,&sigma,&shift,&anal,gaussianSubSetSize,tabPtsValGauss,integrandType);

      /* -------- We then iterate over each points -------- */

    for (int i_pts = 0; i_pts < nbpts-limit; i_pts++) {
       /* -------- We do a certain number of attempt in parallel, using multiple threads -------- */
      #pragma omp parallel for
      for (int i_pt_in_tile = 0; i_pt_in_tile < nbThrow; i_pt_in_tile++) {
         /* -------- We change the seed of the generator -------- */
        generator.seed((i_pt_in_tile*1234+5678)+std::chrono::system_clock::now().time_since_epoch().count());
          /* -------- For each of the attempt, we store the index and the point coordonates in its emplacement in the mseTab -------- */
        mseTab[i_pt_in_tile].index = rAM.at(i_pts);
        for (int currentDim = 1; currentDim <= dimension; currentDim++) {
          mseTab[i_pt_in_tile].point[currentDim-1] = points[rAM.at(i_pts)][currentDim-1];
        }
          double rand = distribution(generator);
          /* -------- For each dimension of the tile, we multiply its vector by a random number between 0 and 1
                      However, we stratify the attempts, here 64 attempts, which explains the division by eight
                      TODO Stratification in other dimension
          -------- */
        for (int i = 0; i < dimension; i++) {
          double shiftInOneDimension = 0.0;
          for (int j = 1; j <= dimension; j++) {
            rand =  distribution(generator);
            shiftInOneDimension +=
            ( ( ( (double)(i_pt_in_tile/8) * v->at(rAM.at(i_pts)).getPreviousDirectionnalVector(j)[i] ) ) / 8.0 )
            + rand * ( (v->at(rAM.at(i_pts)).getPreviousDirectionnalVector(j)[i] / 8.0));
          }
          (mseTab[i_pt_in_tile].point)[i] = (*(v->at(rAM.at(i_pts)).getPreviousRefPoint()))[i] + shiftInOneDimension;
        }
          /* -------- We calculate and store the MSE with this new point -------- */
        mseTab[i_pt_in_tile].apportOfNewPoint = recalculateGaussianValueAllGauss(points[rAM.at(i_pts)],mseTab[i_pt_in_tile].point, &sigma,&shift,tabPtsValGauss,gaussianSubSetSize,nbpts,&anal,integrandType);
      }
        /* -------- We keep the one with the lowest MSE -------- */
      newPointHolder<dimension> theChosenOne = *std::min_element(mseTab+0,mseTab+nbThrow,compareTwoNewPointHolder<dimension>);
        /* -------- We verify wether the new MSE is lower than the previous one stored at the end of the array -------- */
      if (theChosenOne.apportOfNewPoint < tabPtsValGauss[gaussianSubSetSize]) {
        double delta = fabs(theChosenOne.apportOfNewPoint - prevMSE);
        std::cout << outputString << " iter=" << iter_over_pointset << " npts=" << nbpts  << " pts " << theChosenOne.index << " MSE : " << initialSE << " -> " << theChosenOne.apportOfNewPoint << " \t delta : "<< delta << std::endl;
        prevMSE = theChosenOne.apportOfNewPoint;
        /* -------- Since we accepted the change, we modify the point in the tiles and the points array -------- */
        tabPtsValGauss[gaussianSubSetSize] = theChosenOne.apportOfNewPoint;
        changeAllValueGaussTab(points[theChosenOne.index],theChosenOne.point, &sigma,&shift, tabPtsValGauss,gaussianSubSetSize,nbpts,integrandType);
        points[theChosenOne.index][0] = theChosenOne.point[0];
        points[theChosenOne.index][1] = theChosenOne.point[1];
        injectSP(v,&points);
      }
    }
      /* -------- Writes the points in order to save the optimization at regular intervals -------- */
    if (iter_over_pointset % intervalToWrite == 0 ) {
      injectSPInPointHolder(v,pointHolder,nbpts,dimToOpt);
      exportPoints<dimension>(pointHolder,nbPtsTotal,outputString);
    }
  }

    /* -------- Writes the points -------- */
  exportPoints<dimension>(pointHolder,nbPtsTotal,outputString);
  return;
}

int main(int argc, char const *argv[]) {

  /* ----------- Initialization of the variables for CLI11 ----------- */

  double **pointHolder;
  int nbpts = 2;
  int nbThreads = omp_get_max_threads() >= 64 ? 64 : omp_get_max_threads();
  std::string inputString ="pts.dat";
  std::string outputString ="OptimizedPts.dat";
  size_t niters = 1024*1024;
  int gaussianSubSetSize = 4096;
  int nbThrow = 64;
  int integrandType = 2;
  int limit = 1;
  int base = 3;
  bool previous = true;
  std::vector<int>* dimToOpt = new std::vector<int>();
  int innoctave = 6;
  int intervalToWrite = 10;
  int nbTotal = 136;

  /* ----------- CLI11 configuration ----------- */

  CLI::App app { "OptimME" };

  app.add_option("--nbPoints",nbpts,"Number of Points, default: "+std::to_string(nbpts))->required();
  app.add_option("-t,--nbThreads",nbThreads,"Number of threads used , default: "+std::to_string(nbThreads));	//->check(CLI::Range(1,omp_get_max_threads()));
  app.add_option("-n,--iterationNumber",niters,"Number of iterations over the pointset, default: "+std::to_string(niters))->check(CLI::PositiveNumber);
  app.add_option("-i,--input",inputString,"Path to input file, default: "+inputString)->check(CLI::ExistingFile)->required();
  app.add_option("-o,--output",outputString,"Path to output file, default: "+outputString);
  app.add_option("-g",gaussianSubSetSize,"Number of (gaussian) integrands to iterate over in one iteration, default: "+std::to_string(gaussianSubSetSize));
  app.add_option("--integrandType",integrandType,"Type of the integrand to compute MSE : 1|-> HeaviSide, 2|-> SoftEllipses, 3|-> HardEllipses, 4|-> SoftRectangles, 5|-> HardRectangles. Default: "+std::to_string(integrandType));
  app.add_option("-l,--limit",limit,"The limit number behind which we don't modify the points. If left by default, will take every point. Default: "+std::to_string(limit));
  app.add_option("--writingInterval",intervalToWrite,"Defines at which interval should the pointset be written. Default: "+std::to_string(intervalToWrite));
  app.add_option("--base",base,"Base of the pointset. Default: "+std::to_string(base));
  app.add_option("--previous",previous,"Whether to take the previous tile or the current one. Default: "+previous);
  app.add_option("--dimensionToOptimize",(*dimToOpt),"Specify which dimension to Optimize, default: the first two dimensions");
  app.add_option("--innoctave",innoctave,"innoctave of the pointset. Default: "+std::to_string(innoctave));
  app.add_option("--nbTotal",nbTotal,"Number of point in the input file. Default: "+std::to_string(nbTotal));
  CLI11_PARSE(app, argc, argv)


  /* ----------------- Initialization of the variables after the parsing of the arguments --------------- */

  if (dimToOpt->size() == 0) {
    dimToOpt->push_back(0);
    dimToOpt->push_back(1);
  }

  pointHolder = new double*[nbTotal];
  for (int i = 0; i < nbTotal; ++ i) {
    pointHolder[i] = new double[DIM];
  }
  /* ----------------- Initialization of the number of integrands depending of the type and the dimension --------------- */
  switch (integrandType) {
        case INTEGRAND_TYPE_HEAVISIDE:
          switch (DIM) {
            case 2:
              total_N_integrands = INTEGRAND_TYPE_HEAVISIDE_2D_NUMBER;
              break;
            case 3:
              total_N_integrands = INTEGRAND_TYPE_HEAVISIDE_3D_NUMBER;
              break;
            default:
              std::cerr << " L'integrande specifiee n'est pas disponible en dimension " << std::to_string(DIM) << std::endl;
              break;
          }
        case INTEGRAND_TYPE_SOFTELLIPSES:
          switch (DIM) {
            case 2:
              total_N_integrands = INTEGRAND_TYPE_SOFTELLIPSES_2D_NUMBER;
              break;
            case 3:
              total_N_integrands = INTEGRAND_TYPE_SOFTELLIPSES_3D_NUMBER;
              break;
            default:
              std::cerr << " L'integrande specifiee n'est pas disponible en dimension " << std::to_string(DIM) << std::endl;
              break;
          }
          break;
        case INTEGRAND_TYPE_HARDELLIPSES:
          total_N_integrands = INTEGRAND_TYPE_HARDELLIPSES_2D_NUMBER;
          break;
        case INTEGRAND_TYPE_SOFTRECTANGLES:
          total_N_integrands = INTEGRAND_TYPE_SOFTRECTANGLES_2D_NUMBER;
          break;
        case INTEGRAND_TYPE_HARDRECTANGLES:
          total_N_integrands = INTEGRAND_TYPE_HARDRECTANGLES_2D_NUMBER;
          break;
        default:
        std::cerr << "L'integrande specifiee n'est pas disponible." << '\n';
        return -1;
      }
  /* ----------------- OpenMP configuration --------------- */

  omp_set_dynamic(0);
  omp_set_num_threads(nbThreads);

  /* ----------------- We read the points --------------- */
  importPoints<DIM>(pointHolder,inputString,nbTotal);

  /* ----------------- We create the tiles --------------- */

  std::vector<Tiles<DIM>>* v = tilesOnTheFly<DIM>(nbpts, base, previous, dimToOpt, innoctave, pointHolder);

  /* ----------------- We optimize --------------- */
  optimPointME<DIM>(v,nbpts,niters,nbThrow,outputString,gaussianSubSetSize,integrandType,limit-1,intervalToWrite,nbTotal,pointHolder, dimToOpt);

  return 0;
}
