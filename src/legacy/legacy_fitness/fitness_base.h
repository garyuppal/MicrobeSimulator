#ifndef MICROBE_SIMULATOR_FITNESS_BASE_H
#define MICROBE_SIMULATOR_FITNESS_BASE_H

#include <deal.II/base/point.h>
using dealii::Point;

#include <array>

// #include "../utility/argparser.h"
// #include "../chemicals/chemicals.h"


namespace MicrobeSimulator{

	template<int dim, int NumberChemicals>
	class FitnessBase{
	public:
		FitnessBase();
		virtual ~FitnessBase() {}

		virtual double value(const Point<dim>& location, 
			const std::array<double, NumberChemicals>& secretion_rates) const; 

	// can include pointers to chemicals as private variables in specific instantiation
	};


	template<int dim, int NumberChemicals>
	FitnessBase<dim, NumberChemicals>::FitnessBase()
	{}


	template<int dim, int NumberChemicals>
	double
	FitnessBase<dim, NumberChemicals>::value(const Point<dim>& /* location */, 
			const std::array<double, NumberChemicals>& /* secretion_rates */) const
	{
		return 0.;
	}


}


#endif