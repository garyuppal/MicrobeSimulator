#ifndef MICROBESIMULATOR_TEST_FUNCTIONS_H
#define MICROBESIMULATOR_TEST_FUNCTIONS_H

#include "../utility/utility.h"
#include "../utility/parameter_handler.h"
#include "../utility/command_line_parser.h"

#include "../geometry/geometry.h"

#include <array>
#include <algorithm>
#include <memory>

#include <deal.II/base/function.h>

namespace MicrobeSimulator{

	/** \brief exact functions for testing ...*/
namespace TestFunctions{

/** \brief A gaussian function */
template<int dim>
class Gaussian : public dealii::Function<dim>{
public:
	Gaussian(const Point<dim>& c, double a, double w);

	double value(const Point<dim>& p,
		const unsigned int component = 0) const override;
private:
	Point<dim> center;
	double amplitude;
	double width;
};

template<int dim>
Gaussian<dim>::Gaussian(const Point<dim>& c, double a, double w)
	:
	center(c), amplitude(a), width(w)
{}

template<int dim>
double
Gaussian<dim>::value(const Point<dim>& p,
		const unsigned int /* component */) const
{
	return amplitude*std::exp( -(p-center)*(p-center)/width );
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

/** \brief Time dependent gaussian */
/** with possible flow rate, diffusion, decay?, time evolution */
template<int dim> 
class GaussianSolution : public dealii::Function<dim>{
public:
	GaussianSolution(const ParameterHandler& prm, 
		const Geometry<dim>& geo, 
		unsigned int id);

	static void declare_parameters(ParameterHandler& prm);

	// accessors:
	double value(const Point<dim>& p,
		const unsigned int component = 0) const override;

	double getTime() const;

	void setTime(double t);
private:
	Point<dim> bottom_left;
	Point<dim> top_right;

	// actual source:
	Point<dim> center;
	double amplitude;
	double width;
	double velocity;
	double time;

	std::array<BoundaryCondition, dim> boundary_conditions; // for adding possible mirror charges

	double gauss_value(const Point<dim>& p) const; 
	Point<dim> image_point(const Point<dim>& p, unsigned int dimension, unsigned int side) const;
};

// IMPL
// ----------------------------------------------------------------------------------------------------
template<int dim>
GaussianSolution<dim>::GaussianSolution(const ParameterHandler& prm, 
	const Geometry<dim>& geo, 
	unsigned int id)
{
	bottom_left = geo.getBottomLeftPoint();
	top_right = geo.getTopRightPoint();
	boundary_conditions = geo.getBoundaryConditions();

	const std::string section = "Gaussian solutions";

	amplitude = prm.get_double_vector(section, "Amplitude")[id];
	width = prm.get_double_vector(section, "Width")[id];
	velocity = prm.get_double_vector(section, "Width")[id];

	time = 1.0;

	center = prm.get_point_list(section, "Center")[id];
}

template<int dim>
void 
GaussianSolution<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Gaussian solutions");
		prm.declare_entry("Amplitude", "{0,0}", Patterns::List(Patterns::Double()));
		prm.declare_entry("Width", "{0,0}", Patterns::List(Patterns::Double()));
		prm.declare_entry("Velocity", "{0,0}", Patterns::List(Patterns::Double()));
		prm.declare_entry("Center", "{{0,0},{0,0}}", 
			Patterns::List(Patterns::List(Patterns::Double())));
	prm.leave_subsection();
}

// accessors:
template<int dim>
double 
GaussianSolution<dim>::value(const Point<dim>& p,
	const unsigned int /* component */) const
{
	double result = 0;
	result += gauss_value(p);

	for(unsigned int i = 0; i < dim; ++i)
		if(boundary_conditions[i] == BoundaryCondition::REFLECT)
			result += gauss_value(image_point(p, i, 0)) + gauss_value(image_point(p, i, 1));

	return result;
}

template<int dim>
double 
GaussianSolution<dim>::gauss_value(const Point<dim>& p) const
{	
	// add velocity in x direction
	Point<dim> vel_shift = ( (dim==2) ? Point<dim>(velocity*time,0) : Point<dim>(velocity*time,0,0) );

	return amplitude*
			std::pow( (1./(4.*Utility::PI*time)), (double)dim/2.0 )*
			std::exp( -(p-center-vel_shift)*(p-center-vel_shift)/(4.*time*width) );
}

template<int dim>
Point<dim> 
GaussianSolution<dim>::image_point(const Point<dim>& p, unsigned int dimension, unsigned int side) const
{
	Point<dim> image(p);
	// acually wnt two images *** add another parameter to this fucntion or take two points as reference...
	if(side==0)
	{
		const double delta = image[dimension] - bottom_left[dimension];
		image[dimension] = image[dimension] - 2.*delta;
	}
	else if(side==1)
	{
		const double delta = top_right[dimension] - image[dimension];
		image[dimension] = image[dimension] + 2.*delta;
	}
	else
	{
		throw std::runtime_error("not a valid option for side");
	}

	return image;
}

template<int dim>
double 
GaussianSolution<dim>::getTime() const
{
	return time;
}

template<int dim>
void 
GaussianSolution<dim>::setTime(double t)
{
	time = t;
}

} // NAMESPACE: TestFunctions

namespace ExactFunctions{

template<int dim>
class Constant : public dealii::Function<dim>{
public:
	Constant();
	Constant(double v);

	double value(const Point<dim>& p,
		const unsigned int component = 0) const override;

private:
	double val;
};

// IMPL
// ------------------------------------------------------------
template<int dim>
Constant<dim>::Constant()
	:
	val(0)
{}

template<int dim>
Constant<dim>::Constant(double v)
	:
	val(v)
{}

template<int dim>
double
Constant<dim>::value(const Point<dim>& /*p*/,
		const unsigned int /* component */) const
{
	return val;
}


}} // CLOSE NAMESPACES
#endif