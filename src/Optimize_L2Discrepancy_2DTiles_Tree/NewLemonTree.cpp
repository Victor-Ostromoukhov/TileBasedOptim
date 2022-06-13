#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <omp.h>
#include <algorithm>
#include "Lib/Node/Node.hpp"
#include "Lib/Tiles/Tiles.hpp"
#include "CLI11.hpp"

struct newPointHolder {
  std::string placementX;
  std::string placementY;
  double apportOfNewPoint;
} ;

bool compareTwoNewPointHolder(const newPointHolder &a,const newPointHolder &b){
  if (a.apportOfNewPoint < b.apportOfNewPoint) {
    return true;
  }
  return false;
}

struct point {
  int index;
  double xco;
  double yco;
};

/* ----------------- Compute Discrepancy ----------------- */

double evaluateFirstPartPointContribution(std::vector<point>* v,int oldpointindice){
  double somme2 = 0.0;
  double produit = 1.0;
  int nb_pts = v->size();
  for (int i = 1; i <= nb_pts; i++) {
    produit = 1.0;
    produit *= (1-(std::max(v->at(i-1).xco,v->at(oldpointindice).xco)));
    produit *= (1-(std::max(v->at(i-1).yco,v->at(oldpointindice).yco)));
    somme2 += produit;
  }

    for (int i = 1; i <= nb_pts; i++) {
      produit = 1.0;
      if (i-1 != oldpointindice) {
        produit *= (1-(std::max(v->at(i-1).xco,v->at(oldpointindice).xco)));
        produit *= (1-(std::max(v->at(i-1).yco,v->at(oldpointindice).yco)));
        somme2 += produit;
      }
      }

  return ((1.0/double ((pow(nb_pts,2)))) * somme2);
}

double evaluateSecondPartPointContribution(std::vector<point>* v,int oldpointindice){
  double produit = 1.0;
  int nb_pts = v->size();
  produit *= (1 -  pow((v->at(oldpointindice).xco),2));
  produit *= (1 -  pow((v->at(oldpointindice).yco),2));
  return  (   (  (  pow(2, -2+1 )  )   / double(nb_pts) ) * produit  );
}

double computel2StarDiscPointContribution(std::vector<point>* v,int index){
      double firstPart = evaluateFirstPartPointContribution(v,index);
      double secondPart = evaluateSecondPartPointContribution(v, index) ;
      double thirdPart = 1.0/9.0;
      return sqrt(firstPart - secondPart + thirdPart);
}

/* ==================== Evaluate the discrepancy of a PointSet ==================== */

double evaluateFirstPart(std::vector<point>* v, int nb_pts,int dimension){
  double somme1 = 0.0;
  double somme2 = 0.0;
  double produit = 0.0;
  for (int i = 1; i <= nb_pts; i++) {
    somme2 = 0.0;
    for (int j = 1; j <= nb_pts; j++) {
      produit = 1.0;
      produit *= (1 - (std::max(v->at(i-1).xco,v->at(j-1).xco) ) );
      produit *= (1 - (std::max(v->at(i-1).yco,v->at(j-1).yco) ) );
      somme2 += produit;
    }
    somme1 += somme2;
    }
  return ((1.0/double ((pow(nb_pts,2)))) * somme1);
}

double evaluateSecondPart(std::vector<point>* v, int nb_pts,int dimension){
  double somme = 0.0;
  double produit = 0.0;
  for (int i = 1; i <= nb_pts; i++) {
      produit = 1.0;
      produit *= (1 -  pow((v->at(i-1).xco),2));
      produit *= (1 -  pow((v->at(i-1).yco),2));
      somme += produit;
    }
  return  (   (  (  pow(2, -dimension+1 )  )   / nb_pts ) * somme  );
}

double computel2StarDisc(std::vector<point>* v){
    int dimension = 2;
    double firstPart = evaluateFirstPart(v, v->size(),2);
    double secondPart = evaluateSecondPart(v,  v->size(),2) ;
    double thirdPart = + pow(  3  , -dimension  );
    return sqrt(firstPart - secondPart + thirdPart);
}

/* ==================== Transform a pointset or a point ==================== */

std::vector<point>* FromTernPointToPointNormalized(std::vector<PointTernaire>* TernPointVector){
  std::vector<point>* vectorOfPoints = new std::vector<point>;
  point p;
  for ( auto point : (*TernPointVector)){
    p.index = point.index;
    p.xco = point.getXcoAsDouble();
    p.yco = point.getYcoAsDouble();
    vectorOfPoints->push_back(p);
  }
  return vectorOfPoints;
}

std::vector<PointTernaire>* TreeVersaPelis(std::vector<PointTernaire>* vectorOfPoints, Node* rootOfPermTree){
  std::string toRet = "";
  std::vector<PointTernaire>* VersaPelisVector = new std::vector<PointTernaire>;

  Node* currentNode = rootOfPermTree;
  // std::cout << "Please" << '\n';
  for (int i = 0; i < (int)vectorOfPoints->size(); i++) {
    toRet = "";
    currentNode = rootOfPermTree;
    // std::cout << (int)vectorOfPoints->at(i).xco.length() << '\n';
    for (int j = 0; j < (int)vectorOfPoints->at(i).xco.length(); j++) {
      // std::cout << "Explain" << '\n';
      // std::cout << std::stoi(vectorOfPoints->at(i).xco.substr(j,1)) << std::endl;
      // std::cout << "Intervalle" << '\n';
      toRet+= std::to_string(currentNode->children[std::stoi(vectorOfPoints->at(i).getXcoAsConsommable().substr(j,1))]->permutation);
      // currentNode->children[];
      // std::cout << std::stoi(vectorOfPoints->at(i).xco.substr(j,1)) << '\n';
      // std::cout << "Me" << '\n';
      currentNode = currentNode->children[std::stoi(vectorOfPoints->at(i).getXcoAsConsommable().substr(j,1))];
      // std::cout << "PlAse" << '\n';
    }
    std::cout << "Point "  << vectorOfPoints->at(i).index << " Nouveau Point x " << toRet << '\n';
    VersaPelisVector->push_back(PointTernaire(vectorOfPoints->at(i).index,toRet,vectorOfPoints->at(i).yco));
  }
  return VersaPelisVector;
}

std::vector<PointTernaire>* TreesVersaPelis(std::vector<PointTernaire>* vectorOfPoints, Node* treeHolder[2]){
  std::string toRetX = "";
  std::string toRetY = "";

  Node* currentNodeX = treeHolder[0];
  Node* currentNodeY = treeHolder[1];

  std::vector<PointTernaire>* VersaPelisVector = new std::vector<PointTernaire>;

  for (int i = 0; i < (int)vectorOfPoints->size(); i++) {

    toRetX = "";
    toRetY = "";

    currentNodeX = treeHolder[0];
    currentNodeY = treeHolder[1];

    for (int j = 0; j < (int)vectorOfPoints->at(i).xco.length(); j++) {
      toRetX += std::to_string(currentNodeX->children[std::stoi(vectorOfPoints->at(i).getXcoAsConsommable().substr(j,1))]->permutation);
      currentNodeX = currentNodeX->children[std::stoi(vectorOfPoints->at(i).getXcoAsConsommable().substr(j,1))];
    }

    for (int k = 0; k < (int)vectorOfPoints->at(i).yco.length(); k++) {
      toRetY += std::to_string(currentNodeY->children[std::stoi(vectorOfPoints->at(i).getYcoAsConsommable().substr(k,1))]->permutation);
      currentNodeY = currentNodeY->children[std::stoi(vectorOfPoints->at(i).getYcoAsConsommable().substr(k,1))];
    }

    VersaPelisVector->push_back(PointTernaire(vectorOfPoints->at(i).index,toRetX,toRetY));
  }
  return VersaPelisVector;
}

double TernaryPointToNormalizedDouble(std::string ternaryPoint){
  double doublePointToRet = 0.0;
  int i = 1;
  for( auto bit : ternaryPoint){
    doublePointToRet += (double)(bit-48) * pow(3,-i);
    i++;
  }
  return doublePointToRet;
}

/* ==================== Initializing the trees  ==================== */

  /* =============== Initializing the children  =============== */
void addChildrenRecursive(Node* currentNode,int octave){
  if (octave == 0) return;
  for (int currentOffspring = 0; currentOffspring < 3; currentOffspring++) {
    currentNode->addChild(currentOffspring,currentOffspring);
    addChildrenRecursive(currentNode->children[currentOffspring],octave-1);
  }
}

void addChildren(Node* treeHolder[2],int octave){
  Node* currentNode ;
  for (int numberOfTrees = 0; numberOfTrees < 2; numberOfTrees++) {
    currentNode = treeHolder[numberOfTrees];
    addChildrenRecursive(currentNode,octave);
  }
}

/* =============== Initializing the branches  =============== */

void addBranchesRecursive(Node* currentNode,int octave){
  if (octave == 0) return;
  for (int currentOffspring = 0; currentOffspring < 3; currentOffspring++) {
    currentNode->addChild(currentOffspring,currentOffspring);
    addChildrenRecursive(currentNode->children[currentOffspring],octave-1);
  }
}

void addBranches(Node* treeHolder[2],int octave,std::vector<PointTernaire>* vectorOfPoints){
  for ( auto point : (*vectorOfPoints) ){
      treeHolder[0]->addBranch(point.getXcoAsConsommable(),generateAPlacement(point.getXcoAsConsommable().substr(0,octave),treeHolder[0],octave),treeHolder[0]);
      treeHolder[1]->addBranch(point.getYcoAsConsommable(),generateAPlacement(point.getYcoAsConsommable().substr(0,octave),treeHolder[1],octave),treeHolder[1]);

  }
}

/* ==================== Main function, where all of the magic operates ==================== */

void NewLemonTree(Node* treeHolder[2],std::vector<PointTernaire>* vectorOfPoints, int octave,int nbOfIteration,int nbThrow,std::string outputString){
    // J'ai :
      // Les deux arbres
      // Le vecteur de points ternaires
    // Je récupère le nombre de points
    int nb_pts = vectorOfPoints->size();
    // Je transforme les points ternaires normaux en points ternaires versaPelis
    std::vector<PointTernaire>* vectorOfPointsVersaPelis =  TreesVersaPelis(vectorOfPoints,treeHolder);

    // Je récupère le vecteur de points VP normalisés
    std::vector<point>* vectorOfPointsNormalized = FromTernPointToPointNormalized(vectorOfPointsVersaPelis);

    // Je déclare le tableau qui stocke les essais pour le parallélisme
    newPointHolder discTab[nbThrow];

    // Je calcule la discrépance de ce pointset

    double discrepancy = computel2StarDisc(vectorOfPointsNormalized);
    double initialDiscrepancy = discrepancy;

    for (int currentIteration = 0; currentIteration < nbOfIteration; currentIteration++) {
      // Pour chacun des points
      for (int currentPoint = 0; currentPoint < nb_pts; currentPoint++) {
        // Je calcule l'ancienne discrépance en currentPoint

        double discrepancyEnCurrPoint = computel2StarDiscPointContribution(vectorOfPointsNormalized,currentPoint);

        // #pragma omp for
          for (int currentTryInTile = 0; currentTryInTile < nbThrow; currentTryInTile++) {
            //Je re-remplis le vecteur de points normalisés initialisé vide par le private de omp
            std::vector<point>* vectorOfPointsNormalizedLocal(vectorOfPointsNormalized);
            // Je génère des nouvelles coordonnées x et y
            std::string placementX =  generateAPlacement(vectorOfPoints->at(currentPoint).getXcoAsConsommable().substr(0,octave),treeHolder[0],octave);
            std::string placementY =  generateAPlacement(vectorOfPoints->at(currentPoint).getYcoAsConsommable().substr(0,octave),treeHolder[1],octave);

            double placementXD = TernaryPointToNormalizedDouble(placementX);
            double placementYD = TernaryPointToNormalizedDouble(placementY);

            // Je les remplace dans le vecteur de points VPN

            vectorOfPointsNormalizedLocal->at(currentPoint).xco = placementXD;
            vectorOfPointsNormalizedLocal->at(currentPoint).yco = placementYD;

            double discrepancyEnCurrPointNew = computel2StarDiscPointContribution(vectorOfPointsNormalizedLocal,currentPoint);

            // Je stocke les valeurs d'importance
            discTab[currentTryInTile].placementX = placementX;
            discTab[currentTryInTile].placementY = placementY;
            discTab[currentTryInTile].apportOfNewPoint = discrepancyEnCurrPointNew;

          }
        newPointHolder theChosenOne = *std::min_element(discTab+0,discTab+nbThrow,compareTwoNewPointHolder);

        // Je calcule la diff de discpréance

        double newDisc = sqrt(pow(discrepancy,2) - pow(discrepancyEnCurrPoint,2) + pow(theChosenOne.apportOfNewPoint,2));
        if (newDisc < discrepancy) { // Si c'est mieux
            // Je coupe la branche du point et je rajoute la permutaiton
            std::cout << "Amelioration Iteration " << currentIteration <<": Initial discrepancy : " << initialDiscrepancy << ", Point " << currentPoint << " " << discrepancy << " --> " << newDisc << " Delta => " << discrepancy-newDisc <<'\n';
            discrepancy = newDisc;
            // std::cout << "/* message */" <<vectorOfPoints->at(currentPoint).getXcoAsConsommable().substr(0,octave)<< '\n';
            treeHolder[0]->cutBranch(vectorOfPoints->at(currentPoint).getXcoAsConsommable().substr(0,octave),treeHolder[0]);
            treeHolder[1]->cutBranch(vectorOfPoints->at(currentPoint).getYcoAsConsommable().substr(0,octave),treeHolder[1]);
            // std::cout << "m" << vectorOfPoints->at(currentPoint).getXcoAsConsommable() << '\n';
            treeHolder[0]->addBranch(vectorOfPoints->at(currentPoint).getXcoAsConsommable(),theChosenOne.placementX,treeHolder[0]);
            treeHolder[1]->addBranch(vectorOfPoints->at(currentPoint).getYcoAsConsommable(),theChosenOne.placementY,treeHolder[1]);

            vectorOfPointsNormalized->at(currentPoint).xco = TernaryPointToNormalizedDouble(theChosenOne.placementX);
            vectorOfPointsNormalized->at(currentPoint).yco = TernaryPointToNormalizedDouble(theChosenOne.placementY);

        }else{ // Sinon
            // Je remet le bon point dans le 2VPN
            vectorOfPointsNormalized = FromTernPointToPointNormalized(TreesVersaPelis(vectorOfPoints,treeHolder));
        }
      }

      vectorOfPointsVersaPelis = TreesVersaPelis(vectorOfPoints,treeHolder);
      vectorOfPointsNormalized = FromTernPointToPointNormalized(vectorOfPointsVersaPelis);
    }

    // EXPORTY
    exportPoints(vectorOfPointsVersaPelis,outputString);
    // for ( auto pt : (*vectorOfPointsNormalized)){
    //     std::cout << pt.index << " " << std::setprecision (17) << pt.xco << " " << pt.yco << '\n';
    //   }
}

int main(int argc, char const *argv[]) {

  srand (time(NULL));
  int nbIterationsPerTile = 1;
  int nbThreads = omp_get_max_threads() == 64 ? 64 : omp_get_max_threads();
  size_t niters = 1024*1024;
  std::string inputString ="pts.dat";
  std::string outputString ="OptimizedPts.dat";
  int octave = 2;

  CLI::App app { "Optimizes a set of 2D tiles in ordre to reduce its L2 discrepancy using ternary trees" };

  app.add_option("-t,--nbThreads",nbThreads,"Number of threads used , default: "+std::to_string(nbThreads) );
  app.add_option("-n,--iterationNumber",niters,"Number of iterations over the pointset, default: "+std::to_string(niters))->check(CLI::PositiveNumber);
  app.add_option("-i,--input",inputString,"Path to input file, default: "+inputString)->check(CLI::ExistingFile)->required();
  app.add_option("-o,--output",outputString,"Path to output file, default: "+outputString);
  app.add_option("--octave",octave,"Current octave level, default:" +std::to_string(octave));

    CLI11_PARSE(app, argc, argv)
    // =========== OpenMP Configuration =========== //
    omp_set_dynamic(0);
    omp_set_num_threads(nbThreads);

    std::vector<PointTernaire>* vectorOfPoints = extractSamplingPoints(importTiles(inputString,octave));

    Node* treeHolder[2];

    treeHolder[0] = new Node(0,1);
    treeHolder[1] = new Node(0,1);

    addChildren(treeHolder,octave);

    addBranches(treeHolder,octave,vectorOfPoints);

    // treeHolder[0]->printTT("",treeHolder[0],false);

    NewLemonTree(treeHolder,vectorOfPoints,octave,niters,nbIterationsPerTile,outputString);

  return 0;
}
