#include <iostream>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <fstream>
#include "PointTernaire.hpp"



std::string paddingNumber(std::string stringTernary,int totalPaddingNeeded,bool returnBool,bool returnBool2){
  if ((int)stringTernary.length() == totalPaddingNeeded) return stringTernary;
  if (returnBool) {
    reverse(stringTernary.begin(),stringTernary.end());
  }
  while ((int)stringTernary.length() < totalPaddingNeeded) {
    stringTernary += "0";
  }
  if (returnBool2) {
    reverse(stringTernary.begin(),stringTernary.end());
  }
  return stringTernary;
}

PointTernaire::PointTernaire(){
  index = 0;
  xco = "";
  yco = "";
}

PointTernaire::PointTernaire(int indexParam ,std::string xcoParam ,std::string ycoParam){
  index = indexParam;
  xco = xcoParam;
  yco = ycoParam;
}

PointTernaire::PointTernaire(int indexParam ,int xcoParamAsInt ,int ycoParamAsInt,int octave){
  index = indexParam;
  bool testx = true;
  bool testy = true;
  for (int i = 1; i <= octave; i++) {
    if (xcoParamAsInt < pow(3,i) && testx) {
      xco = paddingNumber(decimalIntToTernaryString(xcoParamAsInt),octave,true,true);
      xco = paddingNumber(xco,20,false,true);
      testx = false;
    }
    if (ycoParamAsInt < pow(3,i) && testy) {
      yco = paddingNumber(decimalIntToTernaryString(ycoParamAsInt),octave,true,true);
      yco = paddingNumber(yco,20,false,true);
      testy = false;
    }
  }



}

PointTernaire::~PointTernaire(){

}

double PointTernaire::getXcoAsDouble(){
  double doubleXcoToRet = 0.0;
  int i = 1;
  for( auto bit : xco){
    doubleXcoToRet += (double)(bit-48) * pow(3,-i);
    i++;
  }
  return doubleXcoToRet;
}

double PointTernaire::getYcoAsDouble(){
  double doubleYcoToRet = 0.0;
  int i = 1;
  for( auto bit : yco){
    doubleYcoToRet += (double)(bit-48) * pow(3,-i);
    i++;
  }
  return doubleYcoToRet;
}

std::string PointTernaire::getXcoAsConsommable(){
  return std::string(xco.rbegin(),xco.rend());
  // return xco;
}

std::string PointTernaire::getYcoAsConsommable(){
  return std::string(yco.rbegin(),yco.rend());
  // return yco;
}



int PointTernaire::ternaryStringToDecimalInt(std::string str){
     int base = 3;
     int len = str.size();
     int power = 1;
     int num = 0;
     int i;
     for (i = len - 1; i >= 0; i--) {
        if (std::stoi(str.substr(i,1)) >= base) {
           printf("Invalid Number");
           return -1;
        }
        num += std::stoi(str.substr(i,1)) * power;
        power = power * base;
     }
     return num;
}

std::string PointTernaire::printAsDouble(bool printIndex){
  std::string toRet = "";
  if (printIndex) {
    toRet += std::to_string(index);
    toRet += "\t";
  }
  toRet += std::to_string(getXcoAsDouble());
  toRet += "\t";
  toRet += std::to_string(getYcoAsDouble());
  return toRet;
}



std::string decimalIntToTernaryString(int decimal){
  int base = 3;
  if(decimal == 0) return "0";
  char NUMS[] = "0123456789ABCDEF"; // Characters that may be used
  std::string result = ""; // Create empty string ready for data to be appended
  do{
      result.push_back(NUMS[decimal%base]);
      // Append the appropriate character in NUMS based on the equation decimal%base

      decimal /= base; // Calculate new value of decimal
  }while(decimal != 0); // Do while used for slight optimisation

  return std::string(result.rbegin(), result.rend());
}


void exportPoints(std::vector<PointTernaire>* vectorOfPoints,std::string outputString){
  std::ofstream o;
  o.open(outputString);
  for (std::vector<PointTernaire>::iterator it = vectorOfPoints->begin();it != vectorOfPoints->end();  it++) {
    o << std::setprecision (17) <<(*it).getXcoAsDouble() << "\t" << (*it).getYcoAsDouble() <<  '\n';
  }
  o.close();
}
