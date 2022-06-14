#include <iostream>
#include "Node.hpp"
#include <vector>
#include <random>
#include <algorithm>

Node::Node(){
  permutation = 0;
  direction = 0;
}

Node::Node(unsigned char directionParam, unsigned char permutationParam ){
  direction = directionParam;
  permutation = permutationParam;
  for (int i = 0; i < 3; i++) {
    children[i] = nullptr;
  }
}

Node::~Node(){
  for (int i = 0; i < 3; i++) {
    if (children[i] != nullptr) {
      delete(children[i]);
    }
  }
}

int Node::countChildren(){
  int numberOfChildrenToReturn = 0;
  for (int i = 0; i < 3; i++) {
    if (children[i] != nullptr) {
      numberOfChildrenToReturn++;
    }
  }
  return numberOfChildrenToReturn;
}

void Node::addChild(unsigned char directionParam, unsigned char permutationParam){
  if (children[directionParam] == nullptr) {
    children[directionParam] = new Node(directionParam,permutationParam);
  }
}

void Node::removeChildByDirection(unsigned char directionParam){
  if (children[directionParam] != nullptr) {
    delete(children[directionParam]);
    children[directionParam] = nullptr;
  }
}

void Node::printTT(const std::string& prefix, Node* node, bool isLeft){

  std::cout << prefix;
  if (node == nullptr) {
    std::cout << (isLeft ? "├───*" : "└───*" ) << std::endl;
    return;
  }
  std::cout << (isLeft ? "├─"+std::to_string(node->direction)+"──" : "└─"+std::to_string(node->direction)+"──" );

  std::cout << std::to_string(node->permutation) << std::endl;

  for (int i = 0; i < 3; i++) {
          if (i != 2) {
            printTT( prefix + (isLeft ? "│     " : "      "), node->children[i], true);
          }else{
            printTT( prefix + (isLeft ? "│     " : "      "), node->children[i], false);
          }
        }
}

std::string Node::printForSerialization(){
  return std::to_string(direction) + "-" + std::to_string(permutation) + "(";
}

void Node::cutBranch(std::string branchToCut, Node* root){

  // std::cout << "Pour " << branchToCut << '\n';
  Node* currentNode = root;
  // std::cout << "curr"  << currentNode->children[0]<< '\n';
  // std::cout << "curr"  << currentNode->children[1]<< '\n';
  // std::cout << "curr"  << currentNode->children[2]<< '\n';


  for (int currentTrit = 0; currentTrit < (int) branchToCut.length() ; currentTrit++) {
    // std::cout << "1" << '\n';
    // std::cout << "fr" << std::stoi(branchToCut.substr(currentTrit,1)) << '\n';
    currentNode = currentNode->children[std::stoi(branchToCut.substr(currentTrit,1))];
    // std::cout << "curr"  << currentNode->children[0]<< '\n';
    // std::cout << "curr"  << currentNode->children[1]<< '\n';
    // std::cout << "curr"  << currentNode->children[2]<< '\n';
    // std::cout << "curr march pas" << currentNode->children[std::stoi(branchToCut.substr(branchToCut.length() -1 ,1))] << '\n';
    if (currentNode == nullptr) {
      // std::cout << "2" << '\n';
      return;
    }
  }
  for (int currentChild = 0; currentChild < 3; currentChild++) {
    if (currentNode->children[currentChild] != nullptr) {
      delete(currentNode->children[currentChild]);
      currentNode->children[currentChild] = nullptr;
      // std::cout << "de "<< branchToCut.substr(branchToCut.length() - 1,1) << '\n';
      // std::cout << "3" << '\n';
      // return;
    }
  }
  // std::cout << "curr march pas 2 " << currentNode->children[std::stoi(branchToCut.substr(branchToCut.length() -1 ,1))] << '\n';




}

void Node::addBranch(std::string branchToAddDir, std::string branchToAddPerm, Node* root){
  Node* currentNode = root;
  for (int currentTritDir = 0; currentTritDir < (int) branchToAddDir.length() ; currentTritDir++) {
    if (currentNode->children[std::stoi(branchToAddDir.substr(currentTritDir,1))] == nullptr) {
      currentNode->children[std::stoi(branchToAddDir.substr(currentTritDir,1))] = new Node((unsigned char)branchToAddDir[currentTritDir] - 48,(unsigned char)branchToAddPerm[currentTritDir] - 48);
    }
    currentNode = currentNode->children[std::stoi(branchToAddDir.substr(currentTritDir,1))];
  }
}

unsigned char Node::availablePermutation(){
  bool takenPermutation[3] = {false,false,false};
  for (int childIndex = 0; childIndex < 3; childIndex++){
    if (children[childIndex] != nullptr) {
      takenPermutation[((int) (children[childIndex]->permutation))] = true;
    }
  }
  std::vector<unsigned char> availablePermutationVector;
  for (int permIndex = 0; permIndex < 3; permIndex++){
    if (takenPermutation[permIndex] == false) {
      availablePermutationVector.push_back(permIndex);
    }
  }
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(availablePermutationVector.begin(), availablePermutationVector.end(),g);
  return availablePermutationVector.at(0);
}



std::string serialize(Node* root){
  std::string toExport = "";
  toExport += root->printForSerialization();
  for (int i = 0; i < 3; i++) {
    if (i < (int) root->countChildren()) {
      toExport += serialize(root->children[i]);
    }else{
      toExport += "*";
    }
    if(i != 2)
      toExport += ",";
  }
  toExport += ")";
  return toExport;
}

std::string parse(std::string delimiter,std::string &stringToParse){
  std::string toRet = "";
  size_t pos = 0;
  pos = stringToParse.find(delimiter);
  toRet = stringToParse.substr(0,pos);
    stringToParse = stringToParse.substr(pos + 1, std::string::npos);
  return toRet;
}

void unserialize(std::string &str,Node* root){
    std::string currentNode = "";
    if(str.empty())
        return;
    int dir = std::stoi(str.substr(2,1));
    root->addChild(dir,std::stoi(str.substr(0,1)));
    str = str.substr(4, std::string::npos);
    for (int i = 0; i < 3; i++) {
      if (str.substr(0,1).compare("*") == 0) {
        parse(")",str);
        parse(",",str);
        return;
      }
      unserialize(str,root->children[dir]);
    }
    return;
}

Node* beginUnserialize(std::string str){
  std::string currentNode = "";
  if(str.empty())
      return nullptr;
  currentNode = parse("(",str);
  Node* root = new Node(std::stoi(currentNode.substr(0,1)),std::stoi(currentNode.substr(2,1)));
  for (int NumberOfChildren = 0; NumberOfChildren < 3; NumberOfChildren++) {
    unserialize(str,root);
  }
  std::cout << "Deserialization is a success ... I Hope" << '\n';
  return root;

}

std::string generateAPlacement(std::string numberTernary, Node* root,int octave){
  Node* currentNode = root;
    std::string toRet = "";
// std::cout << "FTR" << numberTernary << '\n';
    for (int currentTrit = 0; currentTrit < octave; currentTrit++) {
      toRet+= std::to_string(currentNode->children[std::stoi(numberTernary.substr(currentTrit,1))]->permutation);
      currentNode = currentNode->children[std::stoi(numberTernary.substr(currentTrit,1))];
    }
// std::cout << "G" << '\n';

    // std::cout << numberTernary.substr(currentTrit,1) << " Gty " << '\n';
    // std::cout << currentNode << '\n';
    // std::cout << "GTY" << '\n';
    // std::cout << currentNode->children[std::stoi(numberTernary.substr(currentTrit,1))] << '\n';
    // std::cout << "currentNode->children[std::stoi(numberTernary.substr(currentTrit,1))]" << '\n';
  // std::cout << currentNode << '\n';
  toRet += std::to_string(currentNode->availablePermutation());
  // std::cout << toRet << '\n';
  while (toRet.length() < 20) {
    toRet += std::to_string(rand() % 3);
    // std::cout << "GG" << '\n';
  }
  return toRet;
}
