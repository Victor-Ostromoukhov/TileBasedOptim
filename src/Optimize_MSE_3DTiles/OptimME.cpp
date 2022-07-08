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
#include "CLI11.hpp"

/* ----------- Déclaration des constantes ----------- */

#ifndef DIM
  #define DIM 2
#endif

#define INTEGRAND_TYPE_HEAVISIDE 1
#define INTEGRAND_TYPE_SOFTELLIPSES 2
#define INTEGRAND_TYPE_HARDELLIPSES 3
#define INTEGRAND_TYPE_SOFTRECTANGLES 4
#define INTEGRAND_TYPE_HARDRECTANGLES 5

#define INTEGRAND_TYPE_HEAVISIDE_2D_NUMBER 1048576
#define INTEGRAND_TYPE_SOFTELLIPSES_2D_NUMBER 524288
#define INTEGRAND_TYPE_HARDELLIPSES_2D_NUMBER 16384
#define INTEGRAND_TYPE_SOFTRECTANGLES_2D_NUMBER 16384
#define INTEGRAND_TYPE_HARDRECTANGLES_2D_NUMBER 16384

#define INTEGRAND_TYPE_HEAVISIDE_3D_NUMBER 16384
#define INTEGRAND_TYPE_SOFTELLIPSES_3D_NUMBER 16384

// Declaration des variables globales

int total_N_integrands;


// Déclaration des structures
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

// Déclaration des tableaux d'intégrandes

extern t_Heaviside2D tab_Heaviside2D[INTEGRAND_TYPE_HEAVISIDE_2D_NUMBER];
extern t_GaussianStruct2D tab_SoftEllipses2D[INTEGRAND_TYPE_SOFTELLIPSES_2D_NUMBER];
extern t_GaussianStruct2D tab_Ellipses2D[INTEGRAND_TYPE_HARDELLIPSES_2D_NUMBER];
extern t_SoftRectanglesStruct2D tab_SoftRectangles2D[INTEGRAND_TYPE_SOFTRECTANGLES_2D_NUMBER];
extern t_RectanglesStruct2D tab_Rectangles2D[INTEGRAND_TYPE_HARDRECTANGLES_2D_NUMBER];
extern t_RectanglesStruct2D tab_Rectangles2D[INTEGRAND_TYPE_HARDRECTANGLES_2D_NUMBER];
extern t_Gauss3D tab_Gauss3D[INTEGRAND_TYPE_SOFTELLIPSES_3D_NUMBER];
extern t_Heaviside3D tab_Heaviside3D[INTEGRAND_TYPE_HEAVISIDE_3D_NUMBER];

/* ----------- Déclaration des structures et leur "méthodes" ----------- */

template<int dimension>
struct newPointHolder {
  int index;
  VecX<dimension> point;
  double apportOfNewPoint;
};

template<int dimension>
bool compareTwoNewPointHolder(const newPointHolder<dimension> &a,const newPointHolder<dimension> &b){
  if (a.apportOfNewPoint < b.apportOfNewPoint) {
    return true;
  }
  return false;
}

/* ----------- Fonctions import ----------- */

template<int dimension>
void importPoints(double pointsTab[][dimension], std::string inputString, int nbpts){
  std::ifstream inputPointsFile(inputString, std::ios::in);
  if (inputPointsFile.is_open()) {
    std::string line;
    int i =0;
    size_t pos = 0;
    std::string token;
    for (int currentPoint = 0; currentPoint < nbpts; currentPoint++) {
      getline(inputPointsFile,line);
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
}

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

template<int dimension>
void exportTilesNextStep(std::vector<Tiles<dimension>>* v,std::string outputStringNextStep,std::vector<std::string>* restOfTheDocument){
  std::ofstream o;
  o.open(outputStringNextStep);
  for (typename std::vector<Tiles<dimension>>::iterator it = v->begin();it != v->end();  it++) {
    o << (*it);
  }
  for (auto i = restOfTheDocument->begin(); i != restOfTheDocument->end(); ++i) {
    o << (*i) << '\n';
  }
  o.close();
}
/* ----------- Fonctions Calcul MSE ----------- */

template<int dimension>
void initializeGaussianVectors(std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,std::vector<double>* anal,int offset,int integrandType,int gaussianSubSetSize){
  switch (integrandType) {
    case 1: // HeaviSide
      if (dimension == 2) {
        for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
          for (int i = 0; i < dimension; ++i) {
              (*shift)[z-(offset*gaussianSubSetSize)][i] = tab_Heaviside2D[z].muDiscotinuity[i];  // Déplacement
          }
          for (int j = 0; j < dimension; ++j) {
            (*sigma)[z-(offset*gaussianSubSetSize)](0, j) = tab_Heaviside2D[z].normal[j]; // Matrice SR
          }
          (*anal)[z-(offset*gaussianSubSetSize)] = tab_Heaviside2D[z].integral; // Valeur analytique
        }
      }else {
        for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
          for (int i = 0; i < dimension; ++i) {
              (*shift)[z-(offset*gaussianSubSetSize)][i] = tab_Heaviside3D[z].muDiscotinuity[i];  // Déplacement
          }
          for (int j = 0; j < dimension; ++j) {
            (*sigma)[z-(offset*gaussianSubSetSize)](0, j) = tab_Heaviside3D[z].normal[j]; // Matrice SR
          }
          (*anal)[z-(offset*gaussianSubSetSize)] = tab_Heaviside3D[z].integral; // Valeur analytique
        }
      }
    break;
    case 2: // SoftEllipses
          if (dimension == 2) {
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
          }else{
            for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
              for (int i = 0; i < dimension; ++i) {
                (*shift)[z-(offset*gaussianSubSetSize)][i] = tab_Gauss3D[z].mu[i];  // Déplacement
              }
            for (int j = 0; j < dimension; ++j) {
              for (int i = 0; i < dimension; ++i) {
                (*sigma)[z-(offset*gaussianSubSetSize)](i, j) = tab_Gauss3D[z].mxCInv[i + j * dimension]; // Matrice SR
              }
            }
            (*anal)[z-(offset*gaussianSubSetSize)] = tab_Gauss3D[z].integral; // Valeur analytique
            }
          }
    break;
    case 3: // HardEllipses

      for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
        for (int i = 0; i < dimension; ++i) {
            (*shift)[z][i] = tab_Ellipses2D[z].mu[i];  // Déplacement
        }
        for (int j = 0; j < dimension; ++j) {
            for (int i = 0; i < dimension; ++i) {
                (*sigma)[z](i, j) = tab_Ellipses2D[z].mxCInv[i + j * dimension]; // Matrice SR
            }
        }
        (*anal)[z] = tab_Ellipses2D[z].integral; // Valeur analytique
      }
    break;
    case 4: // SoftRectangles
        for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
          for (int i = 0; i < dimension; ++i) {
              (*shift)[z][i] = tab_SoftRectangles2D[z].mu[i];  // Déplacement
          }
          for (int j = 0; j < dimension; ++j) {
              (*sigma)[z](0, j) = tab_SoftRectangles2D[z].sigma[j]; // Matrice SR
          }
          (*anal)[z] = tab_SoftRectangles2D[z].integral; // Valeur analytique
        }
    break;
    case 5: // HardRectangles
        for (int z = 0 + offset*gaussianSubSetSize ; z < (offset + 1)*gaussianSubSetSize; z++) {
          for (int i = 0; i < dimension; ++i) {
              (*shift)[z][i] = tab_Rectangles2D[z].mu[i];  // Déplacement
          }
          for (int j = 0; j < dimension; ++j) {
              (*sigma)[z](0, j) = tab_Rectangles2D[z].sigma[j]; // Matrice SR
          }
          (*anal)[z] = tab_Rectangles2D[z].integral; // Valeur analytique
        }
    break;
    default:
      std::cerr << "The chosen integrand type hasn't been implemented yet, so stay tuned." << '\n';
    break;
  }
} // Fulfill the gaussian array with the parameters

template<int dimension>
double recalculateGaussianValue(double oldVal, VecX<dimension> oldPoint,VecX<dimension> newPoint, MatXDynamic sigma,VecXDynamic shift,int nbpts,int integrandType){
  return multivariateGaussianIntegrationPointModif(oldPoint,newPoint,shift,sigma,oldVal,nbpts,integrandType);
} // Recalculate the value of the integration with a point shift

template<int dimension>
void changeAllValueGaussTab(VecX<dimension> oldPoint,VecX<dimension> newPoint, std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,double tabPtsValGauss[],int nbGauss,int nbpts,int integrandType){
  for (int i = 0; i < nbGauss; i++) {
    tabPtsValGauss[i] = recalculateGaussianValue(tabPtsValGauss[i],oldPoint,newPoint,(*sigma)[i],(*shift)[i],nbpts,integrandType);
    }
} // Recalculate the value of the integration with a point shift for all gaussian

template<int dimension>
double recalculateGaussianValueAllGauss(VecX<dimension> oldPoint,VecX<dimension> newPoint, std::vector<MatXDynamic>* sigma,std::vector<VecXDynamic>* shift,double tabPtsValGauss[],int nbGauss,int nbpts,std::vector<double>* anal,int integrandType){
  double val = 0.;
  for (int i = 0; i < nbGauss; i++) {
    val+= pow(((*anal)[i] - recalculateGaussianValue(tabPtsValGauss[i],oldPoint,newPoint,(*sigma)[i],(*shift)[i],nbpts,integrandType)),2);
    }
  val /= nbGauss;
  return val;
}

template<int dimension>
double seOfAPointsetOnOneGaussian(std::vector<VecX<dimension>>* points,MatXDynamic sigma,VecXDynamic shift,double anal,int nbGauss,double tabPtsValGauss[],int indice,int integrandType){
    double val = 0.;
    tabPtsValGauss[indice] = multivariateGaussianIntegration((*points), shift, sigma,integrandType);
    val = pow(( anal - tabPtsValGauss[indice]),2);
    return val;
}

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

/* ----------- Fonctions autre ----------- */

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

/* ----------- Fonction Principale ----------- */
template<int dimension>
double optimPointME(std::vector<Tiles<DIM>>* v,int nbpts,std::string inputString,int niters,int nbThrow,std::string outputString,int gaussianSubSetSize,int integrandType,int limit,std::string outputStringNextStep, std::vector<std::string>* restOfTheDocument){

  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  // =========== Initialisation des variables post-command =========== //

  std::vector<MatXDynamic> sigma(gaussianSubSetSize, MatXDynamic(dimension, dimension)); // Initialisation des matrices SR

  std::vector<VecXDynamic> shift(gaussianSubSetSize, VecXDynamic(dimension)); // Initialisation des vecteurs de déplacements

  std::vector<double> anal(gaussianSubSetSize, 0.); // Initialisation des valeurs calculées (analytiques) des gaussiennes

  // =========== Remplissage des variables précédentes  =========== //
  initializeGaussianVectors<dimension>(&sigma,&shift,&anal,0,integrandType,gaussianSubSetSize);
  // =========== Fin Remplissage des variables précédentes  =========== //

  // =========== Lecture des points  =========== //

    double pointsTab[nbpts][dimension];
    importPoints(pointsTab,inputString,nbpts);

    // =========== Fin lecture  =========== //

    // =========== Début Remplissage du vecteur après lecture =========== //

    std::vector<VecX<DIM>> points = (*extractSP(v));

    // =========== Fin Remplissage du vecteur après lecture  =========== //

    // =========== Calcul de l'erreur d'intégration sur les 1024 gaussiennes pour un pointset  =========== //
    double tabPtsValGauss[gaussianSubSetSize+1];

    double initialSE = mseOfAPointsetOnAllGaussian(&points,&sigma,&shift,&anal,gaussianSubSetSize,tabPtsValGauss,integrandType);
    newPointHolder<dimension> mseTab[nbThrow];

    // =========== Initialisation matrice pour parcours aléatoire des points et des sous-sets d gaussiennes  =========== //

    std::vector<int> rAM = randomAccessMatriceGenerator(nbpts,limit);

    std::vector<int> rAMGaussiennes;


    double prevMSE = 0.;
    for (int  iter_over_pointset = 0;  iter_over_pointset < niters;  iter_over_pointset++) {
      rAM = randomAccessMatriceGenerator(nbpts,limit);
      if (iter_over_pointset % ( total_N_integrands / gaussianSubSetSize ) == 0) {
        rAMGaussiennes = randomAccessMatriceGenerator(( total_N_integrands / gaussianSubSetSize ),0);
      }

      initializeGaussianVectors<dimension>(&sigma,&shift,&anal,rAMGaussiennes.at((iter_over_pointset % (( total_N_integrands / gaussianSubSetSize )))),integrandType,gaussianSubSetSize);
      mseOfAPointsetOnAllGaussian(&points,&sigma,&shift,&anal,gaussianSubSetSize,tabPtsValGauss,integrandType);
      for (int i_pts = 0; i_pts < nbpts-limit; i_pts++) {
        std::cout <<  rAM.at(i_pts) << std::endl;
          #pragma omp parallel for
          for (int i_pt_in_tile = 0; i_pt_in_tile < nbThrow; i_pt_in_tile++) {
            generator.seed((i_pt_in_tile*1234+5678)+std::chrono::system_clock::now().time_since_epoch().count());
            mseTab[i_pt_in_tile].index = rAM.at(i_pts);
            for (int currentDim = 1; currentDim <= dimension; currentDim++) {
                  mseTab[i_pt_in_tile].point[currentDim-1] = points[rAM.at(i_pts)][currentDim-1];
            }
                  double rand = distribution(generator);
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
          mseTab[i_pt_in_tile].apportOfNewPoint = recalculateGaussianValueAllGauss(points[rAM.at(i_pts)],mseTab[i_pt_in_tile].point, &sigma,&shift,tabPtsValGauss,gaussianSubSetSize,nbpts,&anal,integrandType);
            //std::cout << "X ==> " <<  (mseTab[i_pt_in_tile].point)[0] << " et Y ==> " << (mseTab[i_pt_in_tile].point)[1] << std::endl;
          }

              //std::cout <<  mseTab[0].apportOfNewPoint << std::endl;
          newPointHolder<dimension> theChosenOne = *std::min_element(mseTab+0,mseTab+nbThrow,compareTwoNewPointHolder<dimension>);
          if (theChosenOne.apportOfNewPoint < tabPtsValGauss[gaussianSubSetSize]) {

            
            double delta = fabs(theChosenOne.apportOfNewPoint - prevMSE);
            if (rAM.at(i_pts) == 0) {
              std::cout << outputString << " Iteration " << iter_over_pointset << " : " << " MSE : " << initialSE << " -> " << theChosenOne.apportOfNewPoint << " \t delta : "<< delta << std::endl;
            }
            prevMSE = theChosenOne.apportOfNewPoint;
            tabPtsValGauss[gaussianSubSetSize] = theChosenOne.apportOfNewPoint;
            changeAllValueGaussTab(points[theChosenOne.index],theChosenOne.point, &sigma,&shift, tabPtsValGauss,gaussianSubSetSize,nbpts,integrandType);
            for (int i = 0; i < dimension; i++) {
            points[theChosenOne.index][i] = theChosenOne.point[i];
            }
            injectSP(v,&points);
          }
      }
    }

    exportTiles(v,outputString);
    if (outputStringNextStep.compare("OptimizedPts_next_step.dat") != 0) {
      exportTilesNextStep(v,outputStringNextStep,restOfTheDocument);
    }
    //exportPoints(&points,outputString);
    return tabPtsValGauss[gaussianSubSetSize];
  // return 0.0;
}

int main(int argc, char const *argv[]) {

  /* ----------- Début Initialisation des variables ----------- */
  int nbpts = 2;
  int nbThreads = omp_get_max_threads() >= 64 ? 64 : omp_get_max_threads();
  std::string inputString ="pts.dat";
  std::string outputString ="OptimizedPts.dat";
  std::string outputStringNextStep ="OptimizedPts_next_step.dat";
  size_t niters = 1024*1024;
  int gaussianSubSetSize = 4096;
  int nbThrow = 64;
  int integrandType = 1;
  std::string outputStringMSE ="MSE.dat";
  int limit = 1;
  std::vector<std::string>* restOfTheDocument = new std::vector<std::string>;

  //[> ----------- Fin Initialisation des variables ----------- <]

  //// =========== Début CLI11 Configuration =========== //

      CLI::App app { "OptimME" };

      app.add_option("--nbPoints",nbpts,"Number of Points, default: "+std::to_string(nbpts))->required();
      app.add_option("-t,--nbThreads",nbThreads,"Number of threads used , default: "+std::to_string(nbThreads));	//->check(CLI::Range(1,omp_get_max_threads()));
      app.add_option("-n,--iterationNumber",niters,"Number of iterations over the pointset, default: "+std::to_string(niters))->check(CLI::PositiveNumber);
      app.add_option("-i,--input",inputString,"Path to input file, default: "+inputString)->check(CLI::ExistingFile)->required();
      app.add_option("-o,--output",outputString,"Path to output file, default: "+outputString);
      app.add_option("-g",gaussianSubSetSize,"Number of (gaussian) integrands to iterate over in one iteration, default: "+std::to_string(gaussianSubSetSize));
      app.add_option("-m,--writeMSE",outputStringMSE,"If precised, will write the MSE of the pointset by appending it to this file, default: "+outputStringMSE);
      app.add_option("--integrandType",integrandType,"Type of the integrand to compute MSE : 1|-> HeaviSide, 2|-> SoftEllipses, 3|-> HardEllipses, 4|-> SoftRectangles, 5|-> HardRectangles. Default: "+std::to_string(integrandType));
      app.add_option("-l,--limit",limit,"The limit number behind which we don't modify the points. If left by default, will take every point. Default: "+std::to_string(limit));
      app.add_option("--outputNextStep",outputStringNextStep,"Path to output file fot next step, i.e. with all the 729 points. If left by default, it will not write it. default: "+outputStringNextStep);

      CLI11_PARSE(app, argc, argv)

// Initialisation du nombre d'integrandes en fonction du type d'integrande choisi

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

      //// =========== Fin CLI11 Configuration =========== //
                                          //[>***<]
      //// =========== Début OpenMP Configuration =========== //

    omp_set_dynamic(0);
    omp_set_num_threads(nbThreads);

      //// =========== Fin OpenMP Configuration =========== //

    std::vector<Tiles<DIM>>* v = importTiles<DIM>(inputString,nbpts,restOfTheDocument);

        //for (auto i : (*v)) {
          //std::cout << i << std::endl;
        //}
    double val = optimPointME<DIM>(v,nbpts,inputString,niters,nbThrow,outputString,gaussianSubSetSize,integrandType,limit-1,outputStringNextStep,restOfTheDocument);
        // =========== Ecriture de l'erreur associée à un pointset  =========== //
        // std::cout << val << '\n';
    if (outputStringMSE.compare("MSE.dat") != 0) {
      exportValue(nbpts,val,outputStringMSE);
    }
        
        std::cout << DIM << std::endl;
  return 0;
}
