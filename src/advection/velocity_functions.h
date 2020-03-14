#ifndef MICROBESIMULATOR_VELOCITY_FUNCTIONS_H
#define MICROBESIMULATOR_VELOCITY_FUNCTIONS_H

/** A collection of analytic velocity functions for ``simple geometries''
* usually handled with AdvectionHandler class
* ... for cylindrical and vortex velocities, still need to implement 
* cylindrical type geometry and mesh for bacteria and chemicals...
*/

#include "./velocity_interface.h"

namespace MicrobeSimulator{ namespace Velocity{
	using namespace dealii;

// ********************************************************************
/** CONSTANT FLOW
*/

template<int dim>
class Constant : public VelocityInterface<dim>{
public:
	Constant();
	Constant(double r);

	Tensor<1, dim> value(const Point<dim>& location) const override;
	double get_maximum_velocity(double max_coordinate) const override;
	
	double getFlowRate() const;
	void setFlowRate(double r);

private:
	double 	flow_rate;
};

// IMPLEMENTATION:
// --------------------------------------------------------------------
template<int dim>
Constant<dim>::Constant()
	:
	flow_rate(0)
{}

template<int dim>
Constant<dim>::Constant(double r)
	:
	flow_rate(r)
{}

template<int dim>
Tensor<1, dim>
Constant<dim>::value(const Point<dim>& /* location */) const
{
	Tensor<1, dim> v;
	v[0] = flow_rate;
	return v;
}

template<int dim>
double 
Constant<dim>::get_maximum_velocity(double /* max_coordinate */) const
{
	return flow_rate;
}

template<int dim>
double 
Constant<dim>::getFlowRate() const
{
	return flow_rate;
}

template<int dim>
void 
Constant<dim>::setFlowRate(double r)
{
	flow_rate = r;
}

// ********************************************************************
/** COUETTE FLOW
*/

template<int dim>
class Couette : public VelocityInterface<dim>{
public:
	Couette();
	Couette(double s);

	Tensor<1, dim> value(const Point<dim>& location) const override;
	double get_maximum_velocity(double max_coordinate) const override;

	double getShearRate() const;
	void setShearRate(double s) const;

private:
	double shear; // value at p[0] == 1 (dv/dr)
};

// IMPLEMENTATION:
// --------------------------------------------------------------------
template<int dim>
Couette<dim>::Couette()
	:
	shear(0)
{}

template<int dim>
Couette<dim>::Couette(double r)
	:
	shear(r)
{}

template<int dim>
Tensor<1, dim>
Couette<dim>::value(const Point<dim>& location ) const
{
	Tensor<1, dim> v;
	v[0] = shear*location[1]; // need at least 2 dimensions
	return v;
}

template<int dim>
double 
Couette<dim>::get_maximum_velocity(double max_coordinate) const
{
	return shear*std::fabs(max_coordinate);
}

template<int dim>
double 
Couette<dim>::getShearRate() const
{
	return shear;
}

template<int dim>
void 
Couette<dim>::setShearRate(double s) const
{
	shear = s;
}

// ********************************************************************
/** SQUAREPIPE: HAGEN-POISEUILLE FLOW FOR 2D AND 3D SQUARE PIPES
*/

template<int dim>
class SquarePipe : public VelocityInterface<dim>{
public:
	SquarePipe();
	SquarePipe(double h, double v);

	Tensor<1, dim> value(const Point<dim>& location) const override;
	double get_maximum_velocity(double max_coordinate) const override;

	double getHeight() const;
	double getMaxVelocity() const;

	void setHeight(double h);
	void setMaxVelocity(double v);
private:
	double height;
	double vmax;
};

// IMPLEMENTATION:
// --------------------------------------------------------------------
template<int dim>
SquarePipe<dim>::SquarePipe()
	:
	height(1),
	vmax(0)
{}

template<int dim>
SquarePipe<dim>::SquarePipe(double h, double v)
	:
	height(h),
	vmax(v)
{}

// value for 2d:
template<>
Tensor<1, 2>
SquarePipe<2>::value(const Point<2>& location ) const
{
	Point<2> v;
	v[0] = vmax*(1. - (location[1]*location[1])/(height*height) ); 
	return v;
}

template<int dim>
double 
SquarePipe<dim>::get_maximum_velocity(double /* max_coordinate */) const
{
	return vmax;
}

// value for 3d:
template<>
Tensor<1, 3>
SquarePipe<3>::value(const Point<3>& location ) const
{
	Tensor<1, 3> v;
	v[0] = vmax*(1. - (location[1]*location[1])/(height*height) )
				*(1. - (location[2]*location[2])/(height*height) ); 
	return v;
}

template<int dim>
double
SquarePipe<dim>::getHeight() const
{
	return height;
}

template<int dim>
double
SquarePipe<dim>::getMaxVelocity() const
{
	return vmax;
}

template<int dim>
void 
SquarePipe<dim>::setHeight(double h)
{
	height = h;
}

template<int dim>
void 
SquarePipe<dim>::setMaxVelocity(double v)
{
	vmax = v;
}

// ********************************************************************
/** RANKINE VORTEX: Extruded for 3d
*/

template<int dim>
class RankineVortex : public VelocityInterface<dim>{
public:
	RankineVortex();
	RankineVortex(double gamma, double r);

	Tensor<1, dim> value(const Point<dim>& location) const override;
	double get_maximum_velocity(double max_coordinate) const override;

	double getCirculation() const;
	double getRadius() const;

	void setCirculation(double gamma);
	void setRadius(double r);
private:
	double circulation;
	double radius;
};

// IMPLEMENTATION:
template<int dim>
RankineVortex<dim>::RankineVortex()
	:
	circulation(0),
	radius(0)
{}

template<int dim>
RankineVortex<dim>::RankineVortex(double gamma, double r)
	:
	circulation(gamma),
	radius(r)
{}

template<int dim>
Tensor<1, dim> 
RankineVortex<dim>::value(const Point<dim>& location) const
{
	const double dist = sqrt( location*location );
    const double theta = atan2(location[1], location[0]); 
    const double v = ( (dist <= radius) ? circulation*dist/( 2.*numbers::PI
    															*radius*radius )
    							: circulation/( 2.*numbers::PI*dist ) );
    Tensor<1, dim> return_value;
    return_value[0] = v*cos(theta);
    return_value[1] = -v*sin(theta);

	return return_value;
}

template<int dim>
double 
RankineVortex<dim>::get_maximum_velocity(double /* max_coordinate */) const
{
	return circulation/(2.*numbers::PI*radius);
}

template<int dim>
double 
RankineVortex<dim>::getCirculation() const
{
	return circulation;
}

template<int dim>
double 
RankineVortex<dim>::getRadius() const
{
	return radius;
}

template<int dim>
void 
RankineVortex<dim>::setCirculation(double gamma)
{
	circulation = gamma;
}

template<int dim>
void 
RankineVortex<dim>::setRadius(double r)
{
	radius = r;
}

// ********************************************************************
// 3D RECTANGULAR AND CYLINDRICAL PIPES:

// ********************************************************************
/** RECTANGULAR PIPE: HAGEN-POISEUILLE FLOW FOR 3D RECTANGULAR PIPES
*/

template<int dim>
class RectangularPipe : public VelocityInterface<dim>{
public:
	RectangularPipe();
	RectangularPipe(double h, double d, double v);

	Tensor<1, dim> value(const Point<dim>& location) const override;
	double get_maximum_velocity(double max_coordinate) const override;

	double getHeight() const;
	double getDepth() const;
	double getMaxVelocity() const;

	void setHeight(double h);
	void setDepth(double d);
	void setMaxVelocity(double v);
private:
	double height;
	double depth;
	double vmax;
};

// IMPLEMENTATION:
// --------------------------------------------------------------------
template<int dim>
RectangularPipe<dim>::RectangularPipe()
	:
	height(0),
	depth(0),
	vmax(0)
{
	if(dim !=3)
		throw std::runtime_error("Rectangular pipe for 3 dimensions only, use square or cylindrical pipe for 2d");
}

template<int dim>
RectangularPipe<dim>::RectangularPipe(double h, double d, double v)
	:
	height(h),
	depth(d),
	vmax(v)
{
	if(dim !=3)
		throw std::runtime_error("Rectangular pipe for 3 dimensions only, use square or cylindrical pipe for 2d");
}

template<int dim>
Tensor<1, dim> 
RectangularPipe<dim>::value(const Point<dim>& location) const 
{
	Tensor<1, dim> v;
	v[0] = vmax*(1. - (location[1]*location[1])/(height*height) )
				*(1. - (location[2]*location[2])/(depth*depth) ); 
	return v;
}

template<int dim>
double 
RectangularPipe<dim>::get_maximum_velocity(double /* max_coordinate */) const
{
	return vmax;
}

template<int dim>
double 
RectangularPipe<dim>::getHeight() const
{
	return height;
}

template<int dim>
double 
RectangularPipe<dim>::getDepth() const
{
	return depth;
}

template<int dim>
double 
RectangularPipe<dim>::getMaxVelocity() const
{
	return vmax;
}

template<int dim>
void 
RectangularPipe<dim>::setHeight(double h)
{
	height = h;
}

template<int dim>
void 
RectangularPipe<dim>::setDepth(double d)
{
	depth = d;
}

template<int dim>
void 
RectangularPipe<dim>::setMaxVelocity(double v)
{
	vmax = v;
}

// ********************************************************************
/** CYLINDRICAL PIPE: HAGEN-POISEUILLE FLOW FOR 2D AND 3D CYLINDRICAL PIPES
*/
// this function is the same as square pipe for 2d applications
// only differs in 3d where it assumes a cylindrical geometry
template<int dim>
class CylindricalPipe : public VelocityInterface<dim>{
public:
	CylindricalPipe();
	CylindricalPipe(double r, double v);

	Tensor<1, dim> value(const Point<dim>& location) const override;
	double get_maximum_velocity(double max_coordinate) const override;

	double getRadius() const;
	double getMaxVelocity() const;

	void setRadius(double r);
	void setMaxVelocity(double v);
private:
	double radius;
	double vmax;
};

// IMPLEMENTATION:
// -------------------------------------------------------------------
template<int dim>
CylindricalPipe<dim>::CylindricalPipe()
	:
	radius(0),
	vmax(0)
{}

template<int dim>
CylindricalPipe<dim>::CylindricalPipe(double r, double v)
	:
	radius(r),
	vmax(v)
{}

template<>
Tensor<1, 2>
CylindricalPipe<2>::value(const Point<2>& location) const
{
	Tensor<1, 2> v;
	v[0] = vmax*(1. - ((location[1]*location[1])/(radius*radius)) ); 
	return v;
}

template<>
Tensor<1, 3>
CylindricalPipe<3>::value(const Point<3>& location) const
{
	const double rsqr = ( location[1]*location[1] + 
								location[2]*location[2] );
	Tensor<1, 3> v;
	v[0] = vmax*(1. - ((rsqr)/(radius*radius)) ); 
	return v;
}

template<int dim>
double 
CylindricalPipe<dim>::get_maximum_velocity(double /* max_coordinate */) const
{
	return vmax;
}

template<int dim>
double 
CylindricalPipe<dim>::getRadius() const
{
	return radius;
}

template<int dim>
double 
CylindricalPipe<dim>::getMaxVelocity() const
{
	return vmax;
}

template<int dim>
void 
CylindricalPipe<dim>::setRadius(double r)
{
	radius = r;
}

template<int dim>
void 
CylindricalPipe<dim>::setMaxVelocity(double v)
{
	vmax = v;
}


}} // close namespaces
#endif