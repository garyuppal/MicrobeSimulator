#ifndef MICROBESIMULATOR_EXACT_SOLUTIONS_H
#define MICROBESIMULATOR_EXACT_SOLUTIONS_H

#include <deal.II/base/function.h>

namespace MicrobeSimulator{ namespace ExactFunctions{
	using namespace dealii;

/** SCALAR GAUSSIAN FUNCTION
*/
template<int dim>
class Gaussian : public Function<dim>
{
public:
	Gaussian();
	Gaussian(const Point<dim>& c, double h, double w);

    virtual double value (const Point<dim>   &p,
                      const unsigned int  component = 0) const;

    // accessors:
    Point<dim>	getCenter() const;
    double 		getHeight() const;
    double 		getWidth() const;

    // mutators:
    void 	setCenter(const Point<dim>& c);
    void	setHeight(double h);
    void 	setWidth(double w);

private:
	Point<dim>		center;
	double 			height;
	double 			width;
};

// IMPLEMENTATION:
// ---------------------------------------------------------------------
template<int dim>
Gaussian<dim>::Gaussian()
	:
	center(),
	height(1),
	width(1)
{}

template<int dim>
Gaussian<dim>::Gaussian(const Point<dim>& c, double h, double w)
	:
	center(c),
	height(h),
	width(w)
{}

template<int dim>
double 
Gaussian<dim>::value (const Point<dim>   &p,
                  const unsigned int /* component */) const
{
	return height*std::exp( - (p-center)*(p-center) / width );
}

// accessors:
template<int dim>
Point<dim>	
Gaussian<dim>::getCenter() const
{
	return center;
}

template<int dim>
double 		
Gaussian<dim>::getHeight() const 
{
	return height;
}

template<int dim>
double 		
Gaussian<dim>::getWidth() const 
{
	return width;
}

// mutators:
template<int dim>
void 	
Gaussian<dim>::setCenter(const Point<dim>& c) 
{
	center = c;
}

template<int dim>
void	
Gaussian<dim>::setHeight(double h) 
{
	height = h;
}

template<int dim>
void 	
Gaussian<dim>::setWidth(double w) 
{
	width = w;
}

// ------------------------------------------------------------------------------

/** SINE FUNCTION
*/


}} // close namespaces
#endif