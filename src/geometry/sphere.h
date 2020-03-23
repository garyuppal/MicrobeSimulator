#ifndef SPHERE_H
#define SPHERE_H

#include <deal.II/base/point.h>
using dealii::Point;

#include <deal.II/base/tensor.h>
using dealii::Tensor;

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

	bool isInSphere(const Point<dim>& location,
            const double buffer = 0.) const;

	Tensor<1, dim> getSphereNormalVector(
		const Point<dim>& intersection_point) const;
private:
	Point<dim> center;
	double radius;

	// reflecting off spheres: 
	Point<dim> getLineSphereIntersection(const Point<dim>& oldPoint, 
	                                const Point<dim>& newPoint, 
	                                const double buffer = 0.) const;
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
	// distance from center, minus radius
	return p.distance(center) - radius;
}

template<int dim>
void 
Sphere<dim>::reflectPoint(const Point<dim>& oldPoint,
                Point<dim>& newPoint,
                const double buffer) const
{
	const Point<dim> intersection = 
		getLineSphereIntersection(oldPoint,newPoint,buffer);

	const Tensor<1,dim> normal = this->getSphereNormalVector(intersection); 

	const Tensor<1,dim> incident = newPoint - oldPoint;
	const Tensor<1,dim> transmitted = newPoint - intersection;

	Tensor<1,dim> reflected_point;
	reflected_point = incident - 2.0*( incident*normal )*normal;

	// rescale:
	reflected_point *= (transmitted.norm())/(reflected_point.norm());

	// recenter (shift vector origin)
	newPoint = intersection + reflected_point;
}

template<int dim>
bool 
Sphere<dim>::isInSphere(const Point<dim>& location,
            const double buffer) const
{
	if(this->distance_from_border(location) + buffer < 1e-8)
		return true;
	return false;
}

template<int dim>
Tensor<1, dim> 
Sphere<dim>::getSphereNormalVector(const Point<dim>& intersection_point) const
{
	Tensor<1,dim> normal = intersection_point - center; 
	// rescale (normalize):
	normal /= normal.norm(); 
	return normal;
}

template<int dim>
Point<dim> 
Sphere<dim>::getLineSphereIntersection(const Point<dim>& oldPoint, 
                                const Point<dim>& newPoint, 
                                const double buffer) const
{
	const double effective_radius = radius + buffer;
	// direction of line:  ( line = oldPoint + d*direction)
	Tensor<1,dim> direction = newPoint - oldPoint;
	direction /= direction.norm(); // unit vector

	// Joachimsthal's Equation:
	// d = -b +/- sqrt[ b^2 - c] ... a == 1
	const double b = direction*( oldPoint - center );
	const double c = (oldPoint - center)*(oldPoint - center) 
		- effective_radius*effective_radius;
	const double discriminant = b*b - c;

	if(discriminant < 0)
		throw std::runtime_error("Error: Line does not intersect sphere");

	if(discriminant == 0)
		return oldPoint + (-b)*direction;

	const Point<dim> first_intersection = oldPoint + (-b + std::sqrt(discriminant))*direction;
	const Point<dim> second_intersection = oldPoint + (-b - std::sqrt(discriminant))*direction;

	// pick point closest to old point:
	if( oldPoint.distance(first_intersection) < oldPoint.distance(second_intersection) )
		return first_intersection;
	else
		return second_intersection;
}

} // close namespace
#endif 