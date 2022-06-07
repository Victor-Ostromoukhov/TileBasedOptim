#include "Points_2D.h"
#include <iostream>

Points_2D::Points_2D(void)
{

}

Points_2D::Points_2D(double x1, double y1)
{
	x=x1;
	y=y1;
}

Points_2D::~Points_2D(void)
{
}

double Points_2D::get_pos_x ()
{
	return x;
}
double Points_2D::get_pos_y ()
{
	return y;
}
double Points_2D::get_pos_dim(int dimension){
  if (dimension == 1) {
    return x;
  }
  if (dimension == 2) {
    return y;
  }
	return x;
}


void Points_2D::set_pos_x (double x1)
{
	x=x1;
}
 void Points_2D::set_pos_y (double y1)
{
	y=y1;
}
