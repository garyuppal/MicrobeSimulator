#ifndef SPHERE_H
#define SPHERE_H

#include <deal.II/base/point.h>
using dealii::Point;

namespace MicrobeSimulator{

template<int dim>
class Sphere{
public:
	Sphere();
	Sphere(Point<dim> c, double r);

	// accessors:
	Point<dim> getCenter() const;
	double getRadius() const;

	// mutators:
	void setCenter(Point<dim> c);
	void setRadius(double r);

	double distance_from_border(const Point<dim>& p) const;

	void reflectPoint(const Point<dim>& old_point,
	                  Point<dim>& new_point,
	                  const double buffer = 0.) const;

private:
	Point<dim> center;
	double radius;

}; // class Sphere{}

// IMPLEMENTATION
// -------------------------------------------------------------------
template<int dim>
Sphere<dim>::Sphere()
	: center(), radius(0.0)
{}

template<int dim>
Sphere<dim>::Sphere(Point<dim> c, double r)
	: center(c), radius(r)
{}

// accessors:
template<int dim>
Point<dim> 
Sphere<dim>::getCenter() const
{
	return center;
}

template<int dim>
double 
Sphere<dim>::getRadius() const
{
	return radius;
}

// mutators:
template<int dim>
void 
Sphere<dim>::setCenter(Point<dim> c)
{
	center = c;
}

template<int dim>
void 
Sphere<dim>::setRadius(double r)
{
	radius = r;
}

template<int dim>
double 
Sphere<dim>::distance_from_border(const Point<dim>& p) const
{
	return 0;
}

template<int dim>
void 
Sphere<dim>::reflectPoint(const Point<dim>& old_point,
                Point<dim>& new_point,
                const double buffer) const
{

}



} // close namespace
#endif 