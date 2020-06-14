#ifndef MICROBESIMULATOR_VELOCITY_FUNCTIONS_H
#define MICROBESIMULATOR_VELOCITY_FUNCTIONS_H

/** @file A collection of analytic velocity functions for ``simple geometries''
* usually handled with AdvectionHandler class
* @todo for cylindrical and vortex velocities, still need to implement 
*/

#include "./velocity_interface.h"

namespace MicrobeSimulator{ namespace Velocity{
	using namespace dealii;

// -------------------------------------------------------------------------------
// CONSTANT
// -------------------------------------------------------------------------------

/** \brief Constant flow in x direction
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

	void printInfo(std::ostream& out) const override;

private:
	double 	flow_rate;
};

// IMPLEMENTATION:
// --------------------------------------------------------------------
/** \brief Constant flow constructor */
template<int dim>
Constant<dim>::Constant()
	:
	flow_rate(0)
{}

/** \brief Construct constant flow with given flow rate */
template<int dim>
Constant<dim>::Constant(double r)
	:
	flow_rate(r)
{}

/** Return vector value of flow rate at given point */
template<int dim>
Tensor<1, dim>
Constant<dim>::value(const Point<dim>& /* location */) const
{
	Tensor<1, dim> v;
	v[0] = flow_rate;
	return v;
}

/** Return maximum possible flow for given function */
template<int dim>
double 
Constant<dim>::get_maximum_velocity(double /* max_coordinate */) const
{
	return flow_rate;
}

/** Get flow rate for constant flow */
template<int dim>
double 
Constant<dim>::getFlowRate() const
{
	return flow_rate;
}

/** Set flow rate of constant flow */
template<int dim>
void 
Constant<dim>::setFlowRate(double r)
{
	flow_rate = r;
}

/** \brief Print info for constant flow */
template<int dim>
void 
Constant<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\tConstant flow:" << std::endl
		<< Utility::short_line << std::endl
		<< "Flow rate: " << flow_rate << std::endl
		<< Utility::short_line << std::endl << std::endl;
}

// -------------------------------------------------------------------------------
// COUETTE FLOW
// -------------------------------------------------------------------------------

/** \brief Couette flow in x direction for linear shear type flow */
template<int dim>
class Couette : public VelocityInterface<dim>{
public:
	Couette();
	Couette(double s);

	Tensor<1, dim> value(const Point<dim>& location) const override;
	double get_maximum_velocity(double max_coordinate) const override;

	double getShearRate() const;
	void setShearRate(double s) const;

	void printInfo(std::ostream& out) const override;

private:
	double shear; // value at p[0] == 1 (dv/dr)
};

// IMPLEMENTATION:
// --------------------------------------------------------------------

/** \brief Default constructor */
template<int dim>
Couette<dim>::Couette()
	:
	shear(0)
{}

/** \brief Construct Couette flow with given shear rate */
template<int dim>
Couette<dim>::Couette(double r)
	:
	shear(r)
{}

/** \brief Return value of Couette flow function at given location */
template<int dim>
Tensor<1, dim>
Couette<dim>::value(const Point<dim>& location ) const
{
	Tensor<1, dim> v;
	v[0] = shear*location[1]; // need at least 2 dimensions
	return v;
}

/** \brief Get maximum possible velocity of Couette flow */
template<int dim>
double 
Couette<dim>::get_maximum_velocity(double max_coordinate) const
{
	return shear*std::fabs(max_coordinate);
}

/** \brief Get shear rate */
template<int dim>
double 
Couette<dim>::getShearRate() const
{
	return shear;
}

/** \brief Set shear rate */
template<int dim>
void 
Couette<dim>::setShearRate(double s) const
{
	shear = s;
}

/** \brief Print info for Couette flow */
template<int dim>
void 
Couette<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\t Couette flow:" << std::endl
		<< Utility::short_line << std::endl
		<< "Shear rate: " << shear << std::endl
		<< Utility::short_line << std::endl << std::endl;
}

// -------------------------------------------------------------------------------
// SQUAREPIPE: HAGEN-POISEUILLE FLOW FOR 2D AND 3D SQUARE PIPES
// -------------------------------------------------------------------------------

/** \brief Square pipe type flow */
/** Hagen-Poiseuille flow for 2D and 3D with square cross section in x direction */
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

	void printInfo(std::ostream& out) const override;

private:
	double height;
	double vmax;
};

// IMPLEMENTATION:
// --------------------------------------------------------------------
/** \brief Default constructor for square pipe flow */
template<int dim>
SquarePipe<dim>::SquarePipe()
	:
	height(1),
	vmax(0)
{}

/** \brief Square pipe flow constructor with given height and max flow rate */
template<int dim>
SquarePipe<dim>::SquarePipe(double h, double v)
	:
	height(h),
	vmax(v)
{}

/** \brief Specialized value for 2D flow */
template<>
Tensor<1, 2>
SquarePipe<2>::value(const Point<2>& location ) const
{
	Point<2> v;
	v[0] = vmax*(1. - (location[1]*location[1])/(height*height) ); 
	return v;
}

/** \brief Specialized value for 3D flow */
template<>
Tensor<1, 3>
SquarePipe<3>::value(const Point<3>& location ) const
{
	Tensor<1, 3> v;
	v[0] = vmax*(1. - (location[1]*location[1])/(height*height) )
				*(1. - (location[2]*location[2])/(height*height) ); 
	return v;
}

/** \brief Get maximum possible flow rate */
template<int dim>
double 
SquarePipe<dim>::get_maximum_velocity(double /* max_coordinate */) const
{
	return vmax;
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

/** \brief Print info for square pipe type flow */
template<int dim>
void 
SquarePipe<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\t Square pipe (Hagen-Poiseuille) flow:" << std::endl
		<< Utility::short_line << std::endl
		<< "Maximum velocity: " << vmax << std::endl
		<< "Height: " << height << std::endl
		<< Utility::short_line << std::endl << std::endl;
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

	void printInfo(std::ostream& out) const override;
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
    return_value[0] = -v*sin(theta);
    return_value[1] = v*cos(theta);

    // std::cout << "location: " << location 
    // 	<< "dist: " << dist 
    // 	<< "theta: " << theta
    // 	<< "v: " << v 
    // 	<< "return_value: " << return_value << std::endl;

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

template<int dim>
void 
RankineVortex<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\t Rankine vortex flow:" << std::endl
		<< Utility::short_line << std::endl
		<< "Circulation: " << circulation << std::endl
		<< "Radius: " << radius << std::endl
		<< Utility::short_line << std::endl << std::endl;
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

	void printInfo(std::ostream& out) const override;
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

template<int dim>
void 
RectangularPipe<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\t Rectangular pipe flow (3D):" << std::endl
		<< Utility::short_line << std::endl
		<< "Maximum velocity: " << vmax << std::endl
		<< "Height: " << height << std::endl
		<< "Depth: " << depth << std::endl
		<< Utility::short_line << std::endl << std::endl;
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

	void printInfo(std::ostream& out) const override;
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

template<int dim>
void 
CylindricalPipe<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\t Cylindrical pipe flow:" << std::endl
		<< Utility::short_line << std::endl
		<< "Maximum velocity: " << vmax << std::endl
		<< "Radius: " << radius << std::endl
		<< Utility::short_line << std::endl << std::endl;
}


}} // close namespaces
#endif