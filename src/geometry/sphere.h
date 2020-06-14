#ifndef SPHERE_H
#define SPHERE_H

#include <deal.II/base/point.h>
using dealii::Point;

#include <deal.II/base/tensor.h>
using dealii::Tensor;

#include <sstream>

/* Check to see if all modifications can be made solely in this class
* for exterior type boundaries
*
* geometry methods to be checked are:
* checkBoundaries() -- need to have correct reflect point and
	// need to have distance from border give proper distance

need to check reflecting point

*/

namespace MicrobeSimulator{

/** \brief Enum type for obstacle. Whether it is an interior or
* bounding exterior obstacle.
*/
enum class ObstacleType : int
{
	INTERIOR = 0, EXTERIOR = 1
};

/** \brief Get string for obstacle type enum */
std::string getObstacleTypeString(ObstacleType ot)
{
	if(ot == ObstacleType::INTERIOR)
		return "Interior";
	else //if(ot == ObstacleType::EXTERIOR)
		return "Exterior";
}

/** \brief Sphere boundary class */
/** Can act as an interior obstacle or exterior boundary encompasing
* microbes in simulation such as for a vortex
*/
template<int dim>
class Sphere{
public:
	Sphere();
	Sphere(Point<dim> c, double r);
	Sphere(Point<dim> c, double r, ObstacleType ot);

	// accessors:
	Point<dim> getCenter() const;
	double getRadius() const;
	ObstacleType getObstacleType() const;

	// mutators:
	void setCenter(Point<dim> c);
	void setRadius(double r);
	void setObstacleType(ObstacleType ot);

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
	ObstacleType obstacle_type;

	// reflecting off spheres: 
	Point<dim> getLineSphereIntersection(const Point<dim>& oldPoint, 
	                                const Point<dim>& newPoint, 
	                                const double buffer = 0.) const;
}; // class Sphere{}


// IMPLEMENTATION
// -------------------------------------------------------------------

/** \brief Default constructor */
template<int dim>
Sphere<dim>::Sphere()
	: center(), radius(0.0), obstacle_type(ObstacleType::INTERIOR)
{}

/** \brief Constructor with dimensions, default to INTERIOR type */
template<int dim>
Sphere<dim>::Sphere(Point<dim> c, double r)
	: center(c), radius(r), obstacle_type(ObstacleType::INTERIOR)
{}

/** \brief Constructor specifying dimensions and type */
template<int dim>
Sphere<dim>::Sphere(Point<dim> c, double r, ObstacleType ot)
	: center(c), radius(r), obstacle_type(ot)
{}

// accessors:
/** \brief Returns center point of sphere */
template<int dim>
Point<dim> 
Sphere<dim>::getCenter() const
{
	return center;
}

/** \brief Return sphere radius */
template<int dim>
double 
Sphere<dim>::getRadius() const
{
	return radius;
}

/** \brief Returns obstacle type, INTERIOR or EXTERIOR */
template<int dim>
ObstacleType
Sphere<dim>::getObstacleType() const
{
	return obstacle_type;
}

// mutators:

/** \brief Set center point of sphere */
template<int dim>
void 
Sphere<dim>::setCenter(Point<dim> c)
{
	center = c;
}

/** \brief Set radius of sphere */
template<int dim>
void 
Sphere<dim>::setRadius(double r)
{
	radius = r;
}

/** \brief Set obstacle type of sphere */
template<int dim>
void 
Sphere<dim>::setObstacleType(ObstacleType ot)
{
	obstacle_type = ot;
}

/** \brief  Find distance from border of sphere */
/** Gives distance from point outside of sphere to boundary for interior type object
* and distance from point inside to boundary for exterior type object. Therefore,
* the distance for a point inside for an interior type sphere will be negative
* and the distance for a point ouside an exterior type will be negative. 
*/
template<int dim>
double 
Sphere<dim>::distance_from_border(const Point<dim>& p) const
{

	if(obstacle_type == ObstacleType::EXTERIOR)
		return radius - p.distance(center);

	// distance from center, minus radius
	return p.distance(center) - radius; 
}

/** \brief Reflects point off of sphere */
/** OldPoint is an original location and newPoint is an attempt at a new location.
* If the new location crosses the obstacle boundary, the point will be reflected back
* in a direction depending on if the obstacle is of INTERIOR or EXTERIOR type.
*/
template<int dim>
void 
Sphere<dim>::reflectPoint(const Point<dim>& oldPoint,
                Point<dim>& newPoint,
                const double buffer) const
{
	/** @todo think about if interior, also note normal is swapped if interior */

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

/** \brief Returns true if point is inside an interior sphere or if point is
outside an exterior type sphere */
template<int dim>
bool 
Sphere<dim>::isInSphere(const Point<dim>& location,
            const double buffer) const
{
	if(this->distance_from_border(location) - buffer < 1e-8)
		return true;
	return false;
}

/** \brief Returns vector normal to sphere surface pointing outward for
* an interior sphere and inward for exterior sphere
*/
template<int dim>
Tensor<1, dim> 
Sphere<dim>::getSphereNormalVector(const Point<dim>& intersection_point) const
{
	Tensor<1,dim> normal = intersection_point - center; 
	// rescale (normalize):
	normal /= normal.norm(); 

	if(obstacle_type == ObstacleType::EXTERIOR)
		normal *= -1; // flip if exterior obstacle
	return normal;
}

/** \brief Helper method to calculate reflection of point off of sphere. Calculates
* the intersection of a line with the sphere and returns the point closest to the 
* end point of the line, that is the oldPoint. 
*/
template<int dim>
Point<dim> 
Sphere<dim>::getLineSphereIntersection(const Point<dim>& oldPoint, 
                                const Point<dim>& newPoint, 
                                const double buffer) const
{
	const double effective_radius = (obstacle_type == ObstacleType::EXTERIOR)?
									radius - buffer :
									radius + buffer;

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
	{
		std::ostringstream msg;

		msg << "Error: Line does not intersect sphere: \n"
			<< "old point: " << oldPoint
			<< " new point: " << newPoint
			<< " buffer: " << buffer 
			<< "\n sphere: center = " << center
			<< " radius = " << radius << "\n";
		throw std::runtime_error(msg.str());
	}

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