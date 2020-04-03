#ifndef MICROBESIMULATOR_TEST_FUNCTIONS_H
#define MICROBESIMULATOR_TEST_FUNCTIONS_H

#include "../utility/parameter_handler.h"
#include "../utility/command_line_parser.h"

#include <array>
#include <algorithm>
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


}} // CLOSE NAMESPACES
#endif