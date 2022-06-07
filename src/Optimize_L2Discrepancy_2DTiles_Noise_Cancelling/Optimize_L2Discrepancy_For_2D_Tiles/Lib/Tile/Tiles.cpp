#include "Tiles.h"
#include <iostream>

Tiles::Tiles(void)
{

}

Tiles::Tiles(std::string bSF,Points_2D sP,Points_2D pRF,Points_2D pv1, Points_2D pv2, std::string aPv2)
{
	beforeSamplingPoint = bSF;
	samplingPoint = sP;
	previousRefPoint = pRF;
	previousv1 = pv1;
	previousv2 = pv2;
	afterPreviousv2 = aPv2;
}

Tiles::~Tiles(void)
{

}
Points_2D Tiles::getSamplingPoint(){
	return samplingPoint;
};

void Tiles::setSamplingPoint(Points_2D newSP){
	samplingPoint = newSP;
}
Points_2D Tiles::getPreviousRefPoint(){
	return previousRefPoint;
};
Points_2D Tiles::getPreviousv1(){
	return previousv1;
};
Points_2D Tiles::getPreviousv2(){
	return previousv2;
};
