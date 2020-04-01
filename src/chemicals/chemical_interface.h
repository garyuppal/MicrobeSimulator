#ifndef MICROBE_SIMULATOR_CHEMICAL_INTERFACE_H
#define MICROBE_SIMULATOR_CHEMICAL_INTERFACE_H

#include <deal.II/base/point.h>
using dealii::Point;
#include <deal.II/base/function.h>

#include <vector>

namespace MicrobeSimulator{ namespace Chemicals{

	template<int dim>
	class ChemicalInterface
	{
	public:
		ChemicalInterface() {}
		virtual ~ChemicalInterface() {} 

		// interface functions:
		virtual double value(const Point<dim>& p) const = 0;  
		virtual std::vector<double> 
			value_list(const std::vector<Point<dim> >& points) const = 0;

		virtual void project_function(const Function<dim>& initial_condition) = 0;

		// no secretion:
		virtual void update() = 0;
		// secretion, but no control:
		virtual void update(const std::vector<Point<dim> >& locations, 
							const std::vector<double>& amounts) = 0;
		// secretion and function: (could do zero amount to do just function)
		virtual void update(const std::vector<Point<dim> >& locations, 
							const std::vector<double>& amounts,
							const Function<dim>& control_function) = 0;

		virtual double getTimeStep() const = 0;

		virtual double getMass() const = 0;

		virtual void print(std::ostream& out, unsigned int chem_id) const = 0; 
		virtual void printInfo(std::ostream& out) const = 0;
	}; 

}}// close namespace
#endif // chemical_interface.h