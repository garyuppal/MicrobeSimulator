#pragma once

#include <deal.II/base/point.h>
using dealii::Point;

#include <deal.II/base/tensor.h>
using dealii::Tensor;

#include <iostream>

namespace MicrobeSimulator{

/** @todo have obstacle base class */

/** Line class for boundaries in 2D domain */
template<int dim>
class Line{
public:
	/** \brief Direction from line in which points are valid */
	enum Orientation{
		ABOVE, BELOW
	};

	Line(const Point<dim>& lft, const Point<dim>& rgt);
	Line(const Point<dim>& lft, const Point<dim>& rgt, Orientation ori);

	// accessors:
	Point<dim> getLeftPoint() const;
	Point<dim> getRightPoint() const;
	Orientation getOrientation() const;

	// methods: (should maybe override virtual methods)
	double distance_from_line(const Point<dim>& p) const; 
	// signed distance, return infinity if not orthogonal to line segment

	bool is_in_bounds(const Point<dim>& p, const double buffer=0) const;
	bool isInBoundingBox(const Point<dim>& p) const;

	
	Tensor<1, dim> getNormalVector(const Point<dim>& p) const;
	
	void reflectPoint(const Point<dim>& old_point, 
				Point<dim>& new_point, const double buffer=0.) const;

	void print(std::ostream& out) const;
	void printInfo(std::ostream& out) const;

private:
	Point<dim> left;
	Point<dim> right;

	Orientation	orientation;

	// store locally mx+b formulation of line
	double slope;
	double intercept;

	// slope of any line perpendicular to this
	double perpendicular_slope;

	Tensor<1, dim> normal; // pointed towards interior

	// for infinity, use:
	// std::numeric_limits<double>::infinity()

	// set slopes:
	void init();

	// orthogonal projection of point onto line:
	Point<dim> getProjectedPoint(const Point<dim>& p) const;
	double getRelativePosition(const Point<dim>& p) const; // return inf if off line, otherwise +-1
	bool isOffLine(const Point<dim>& p) const;

};

// IMPL
// ----------------------------------------------------------------

/** \brief Constructor from points, default orientation is above */
template<int dim>
Line<dim>::Line(const Point<dim>& lft, const Point<dim>& rgt)
	:
	left(lft),
	right(rgt),
	orientation(ABOVE)
{
	init();
}

/** \brief Constructor from points and given orientation */
template<int dim>
Line<dim>::Line(const Point<dim>& lft, const Point<dim>& rgt, Orientation ori)
	:
	left(lft),
	right(rgt),
	orientation(ori)
{
	init();
}

/** \brief Set local values of line slope and slope of perpendicular line */
template<int dim>
void
Line<dim>::init()
{
	// 1) check left < right
	if( left[0] > (right[0] + 1e-14) )
		throw std::runtime_error("lines should be initialized with left x-coordinate "
			"less than or equal to right x-coordinate");

	// 2) set slope and perpendicular_slope
	if( std::fabs(right[0] - left[0]) < 1e-8 )
	{
		if( left[1] > right[1] )
			throw std::runtime_error("vertical lines should be initialized with left y-coordinate "
				"less than or equal to right y-coordinate");

		slope = std::numeric_limits<double>::infinity();

		intercept = std::numeric_limits<double>::infinity(); // or nan, or store x intercept?

	} // if x same, infinite slope
	else
	{
		slope = (right[1] - left[1])/(right[0] - left[0]); // Dy/Dx

		// can use either point, using left here
		intercept = left[1] - slope*left[0];
	} 

	// 3) set perpendicular slope:
	if( std::fabs(slope) < 1e-8)
	{
		perpendicular_slope = std::numeric_limits<double>::infinity();

		// normal points up unless below:
		normal[0] = 0;
		if(orientation == BELOW)
			normal[1] = -1;
		else
			normal[1] = 1;
	}
	else if(std::isinf(slope))
	{
		perpendicular_slope = 0;

		// normal points left (-1) unless below orientation:
		normal[1] = 0;
		if(orientation == BELOW)
			normal[0] = 1;
		else
			normal[0] = -1;
	}
	else
	{
		// std::cout << "setting perpendicular_slope..." << std::endl;
		perpendicular_slope = -1.0/slope;
		// std::cout << "slope is: " << slope << " perp slope is " << perpendicular_slope << std::endl;

		normal[0] = -1;
		normal[1] = -perpendicular_slope; // should have right sign?

		if(perpendicular_slope < 0)
			normal = -normal;

		normal /= std::fabs(normal.norm());

		if(orientation == ABOVE)
			normal = -normal;
	}
}

// accessors:
template<int dim>
Point<dim> 
Line<dim>::getLeftPoint() const
{
	return left;
}

template<int dim>
Point<dim> 
Line<dim>::getRightPoint() const
{
	return right;
}

template<int dim>
typename Line<dim>::Orientation 
Line<dim>::getOrientation() const
{
	return orientation;
}

// helper methods:

/** \brief Return orthogonal projection of point on line */
template<int dim>
Point<dim> 
Line<dim>::getProjectedPoint(const Point<dim>& p) const
{
	// if zero slope, projected point is x component:
	if( std::fabs(slope) < 1e-8)
		return Point<dim>(p[0], left[1]); // left[1] == right[1]

	// if slope infinite, projected point is y component:
	if( std::isinf(slope) )
		return Point<dim>(left[0], p[1]); // left[0] == right[0]

	// given perpendicular line, find intercept for eqn:
	const double b = p[1] - perpendicular_slope*p[0]; // unless slope is inf...

	// get projected point:
	const double x = (intercept - b)/(perpendicular_slope - slope); 
	const double y = slope*x + intercept;

	return Point<dim>(x,y);
}

/** \brief Return +1 if above line segment, -1 if below, and infinity if off line */
/** could also return 0 if on line */
template<int dim>
double 
Line<dim>::getRelativePosition(const Point<dim>& p) const
{
	const Point<dim> projected = getProjectedPoint(p);

	if(isOffLine(projected))
		return std::numeric_limits<double>::infinity();

	// if point is on line:
	if(p.distance(projected) < 1e-14)
		return 0; 

	// get vectors and return cross product sign
	const Tensor<1, dim> line = right - left;
	const Tensor<1, dim> point = p - projected;
	// z component of cross product (line x point):
	const double sign = ( ((line[0]*point[1] - line[1]*point[0]) > 0) ? 1.0 : -1.0 );

	if(orientation == BELOW)
		return -sign;

	return sign;
}

/** \brief Return true if point is not in orthogonal sweep of line segment */
template<int dim>
bool
Line<dim>::isOffLine(const Point<dim>& p) const
{
	// check if given point lies between left and right endpoints
	// given point is already projected onto the line, 
	// it suffices to check if x component is in between left and right (*** need points in right order...)
	// for finite slope, or check y component if infinite

	if( std::isinf(slope) )
		if( p[1] < left[1] || p[1] > right[1] )
			return true;
		else
			return false;

	if( p[0] < left[0] || p[0] > right[0] )
		return true;

	return false;
}

// methods:

/** \brief Normal distance from line to given point */
template<int dim>
double 
Line<dim>::distance_from_line(const Point<dim>& p) const 
{
	const double rel_pos = getRelativePosition(p);
	if( std::isinf(rel_pos) )
		return std::numeric_limits<double>::infinity(); // in case get 0 * inf

	return rel_pos*p.distance(getProjectedPoint(p)); 
}

/** \brief Return whether or not given point is on valid side of line */
template<int dim>
bool 
Line<dim>::is_in_bounds(const Point<dim>& p, const double buffer) const
{
	// first simply check line segment's bounding box 
		// (also to help with possible clashing boundary lines)
	if( !isInBoundingBox(p) )
		return true; // is not in bounding box, then is in bounds relative to this segment

	const double dist = distance_from_line(p);

	// is in bounds relative to this line if dist is infinite, 
	// (no way to handle point if infinite)
	if( std::isinf(dist) )
		return true;

	// if finite, is in bounds only if positive and greater than buffer
	if( dist < buffer)
		return false;

	return true;
}

template<int dim>
bool
Line<dim>::isInBoundingBox(const Point<dim>& p) const
{
	// what if slope is zero or infinite... // pass for now, so is in box
	if( ( std::fabs(slope) < 1e-8) || ( std::fabs(perpendicular_slope) < 1e-8) )
		return true;

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
	{
		const double upper = std::max( left[dim_itr], right[dim_itr] );
		const double lower = std::min( left[dim_itr], right[dim_itr] );

		// check if between left and right for all dimensions
		if( (p[dim_itr] < lower) || (p[dim_itr] > upper) )
			return false;
	}

	return true;
}

// /** \brief Return unit vector normal to line in direction of given point,
// *  with sign given from orientation */
template<int dim>
Tensor<1, dim> 
Line<dim>::getNormalVector(const Point<dim>& /*p*/) const
{
	return normal;
}
// 	Tensor<1, dim> normal;
// 	const double rel_pos = getRelativePosition(p);

// 	// use perpendicular slope:
	
// 	// if inf:
// 	if( std::isinf(slope) )
// 	{
// 		if( std::fabs(rel_pos) < 1e-8 )
// 		{
// 			if(orientation == BELOW)
// 			{
// 				normal[0] = -1;
// 				normal[1] = 0;
// 				return normal;
// 			}
// 			normal[0] = 1;
// 			normal[1] = 0;
// 			return normal;
// 		}
// 		normal[0] = 1;
// 		normal[1] = 0;
// 		return rel_pos*normal; //Tensor<1, dim>(1, 0);
// 	}

// 	// if zero:
// 	if( std::fabs(slope) < 1e-8 )
// 	{
// 		if(std::fabs(rel_pos) < 1e-8)
// 		{
// 			if(orientation == BELOW)
// 			{
// 				normal[0] = 0;
// 				normal[1] = -1;
// 				return normal; //Tensor<1, dim>(0, -1);
// 			}
// 			normal[0] = 0;
// 			normal[1] = 1;
// 			return normal; //Tensor<1, dim>(0, 1);
// 		}
// 		normal[0] = 0;
// 		normal[1] = 1;
// 		return rel_pos*normal; //Tensor<1, dim>(0, 1);
// 	}

// 	// else finite:
// 	// y = mx + b, b doesnt matter here
// 	// use x = 1, so y = perpendicular_slope
// 	// Tensor<1, dim> normal; //(1, perpendicular_slope);
// 	normal[0] = 1;
// 	normal[1] = perpendicular_slope;
// 	normal /= normal.norm(); // normalize

// 	if(std::fabs(rel_pos) < 1e-8)
// 	{
// 		if(orientation == BELOW)
// 			return -normal;
// 		return normal;
// 	}

// 	return rel_pos*normal;
// } // getNormalVector()



// 	// parallel vector is from left to right (right-left)
// 	for(unsigned int dim = 0; dim < 2; ++dim)
// 		normal[dim] = right[dim] - left[dim]; // parallel 

// 	// // rotate:
// 	// // 0 -1					 0 1
// 	// // 1  0 for above and 	-1 0  for below
// 	const double temp = normal[0];
// 	normal[0] = -normal[1];
// 	normal[1] = temp;

// 	// normalize:
// 	normal /= normal.norm(); // can also store this...

// 	if(getRelativePosition(p) < 1e-8)
// 		return normal; // still may want normal vector if on line

// 	return getRelativePosition(p)*normal; 
// } // getNormalVector()

/** \brief Reflect point from line */
template<int dim>
void 
Line<dim>::reflectPoint(const Point<dim>& /* old_point */, 
			Point<dim>& new_point, const double buffer) const
{
	// assuming already checked that old point and new point 
	// are on different sides of the line

	// static int count = 0;

	// if(slope != 0)
	// 	std::cout << "old_point: " << old_point << " try: " << new_point << std::endl;

	const double delta = distance_from_line(new_point) - buffer; // should be negative
	// Tensor<1, dim> normal = getNormalVector(new_point); // should be in direction of old_point
	new_point = new_point - 2.0*delta*normal;

	// if(slope != 0)
	// {
	// 	std::cout << "new_point: " << new_point << " delta: " << delta << " normal: " << normal << std::endl;
		
	// 	++count; 

	// 	if(count > 100)
	// 		throw std::runtime_error("check");
	// }
}

/** \brief Print line end points and orientation */
template<int dim>
void
Line<dim>::print(std::ostream& out) const
{
	const std::string ori_str = 
		( (orientation==ABOVE) ? "ABOVE" : "BELOW" );

	out << left << " " << right << " " << ori_str << std::endl; 
}

/** \brief Print detailed line info */
template<int dim>
void
Line<dim>::printInfo(std::ostream& out) const
{
	const std::string ori_str = 
		( (orientation==ABOVE) ? "ABOVE" : "BELOW" );

	out << "---------------------------------------------------" << std::endl
		<< "\t LINE:" << std::endl
		<< "---------------------------------------------------" << std::endl
		<< "Left end: " << left << std::endl
		<< "Right end: " << right << std::endl
		<< "Orientation: " << ori_str << std::endl
		<< "Slope: " << slope << std::endl
		<< "Intercept: " << intercept << std::endl
		<< "Perpendicular slope: " << perpendicular_slope << std::endl
		<< "Normal: " << normal << std::endl
		<< "---------------------------------------------------"
		<< std::endl << std::endl;
}

} // CLOSE NAMESPACE
/* line.h */