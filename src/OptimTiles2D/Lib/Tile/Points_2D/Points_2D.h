#include <iostream>
class Points_2D
{

private:
	double x;
	double y;

public:
	Points_2D();
	Points_2D(double,double);
	~Points_2D(void);

	double get_pos_x ();
	double get_pos_y ();
  double get_pos_dim(int);

	void set_pos_x (double);
	void set_pos_y (double);
  friend std::ostream& operator<<(std::ostream& os, Points_2D& p) {
		os  << p.get_pos_x() << "\t"<< p.get_pos_y();
	 	return os;
  }
};
