#ifndef MICROBESIMULATOR_DISCRETE_FUNCTION_H
#define MICROBESIMULATOR_DISCRETE_FUNCTION_H

#include <deal.II/base/point.h>
using dealii::Point; // make indep of deal ii -- define your own point class...


namespace MicrobeSimulator{ 

	template<int dim>
	class DiscreteFunction{
	public:
		DiscreteFunction() {}
		virtual ~DiscreteFunction() {}

		virtual double value(const Point<dim>& ) const = 0;

	// private:
		// Point<dim> center;
	};

}

#endif


