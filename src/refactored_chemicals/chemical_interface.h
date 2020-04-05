#ifndef MICROBE_SIMULATOR_REFACTORED_CHEMICAL_INTERFACE_H
#define MICROBE_SIMULATOR_REFACTORED_CHEMICAL_INTERFACE_H

#include <deal.II/base/point.h>
using dealii::Point;

#include <deal.II/base/function.h>
using dealii::Function;

#include <vector>

#include "../utility/parameter_handler.h"
#include "./control_functions.h"

namespace MicrobeSimulator{ 
	/** \brief Namespace for chemical related classes and methods */
	namespace RefactoredChemicals{

/** \brief Base class for chemicals. Provides interface needed for chemical
* fields regardless of how they are implemented */
template<int dim>
class ChemicalInterface
{
public:
	ChemicalInterface(const ParameterHandler& prm, double dt, unsigned int id);
	virtual ~ChemicalInterface() {} 

	// interface functions:
	virtual double value(const Point<dim>& p) const = 0;  
	virtual std::vector<double> 
		value_list(const std::vector<Point<dim> >& points) const = 0;

	// project continuous function onto discrete field
	virtual void project_function(const Function<dim>& initial_condition) = 0;

	// updating methods:

	// no secretion:
	virtual void update() = 0;
	// secretion, but no control:
	virtual void update(const std::vector<Point<dim> >& locations, 
						const std::vector<double>& amounts) = 0;
	// secretion and function: (could do zero amount to do just function)
	virtual void update(const std::vector<Point<dim> >& locations, 
						const std::vector<double>& amounts,
						const Function<dim>& control_function) = 0;

	// accessors:
	double getDiffusionConstant() const;
	double getDecayConstant() const;
	double getTimeStep() const;
	unsigned int getID() const;

	// integrate field over total volume 
	virtual double getMass() const = 0;

	// output
	virtual void print(std::ostream& out) const = 0; 
	virtual void printInfo(std::ostream& out) const = 0;

protected:
	double diffusion_constant;
	double decay_constant;
	double time_step;
	unsigned int chemical_id;
}; 


// IMPL
// ------------------------------------------------------------------------------------

/** \brief ChemicalInterface base constructor */
template<int dim>
ChemicalInterface<dim>::ChemicalInterface(
	const ParameterHandler& prm, double dt, unsigned int id)
	:
	diffusion_constant(0), decay_constant(0), time_step(dt), chemical_id(id)
{
	/** @todo wont work if we want a mix of implementations, 
	* need to modify parameters to have set for each chemical */
	const std::string section = "Chemicals";

	diffusion_constant = prm.get_double_vector(section, "Diffusion")[id];
	decay_constant = prm.get_double_vector(section, "Decay rate")[id];
}

/** \brief Return diffusion constant */
template<int dim>
double 
ChemicalInterface<dim>::getDiffusionConstant() const 
{
	return diffusion_constant;
}

/** \brief Return decay constant */
template<int dim>
double 
ChemicalInterface<dim>::getDecayConstant() const
{
	return decay_constant;
}

/** \brief Return time step */
template<int dim>
double 
ChemicalInterface<dim>::getTimeStep() const 
{
	return time_step;
}

/** \brief Return chemical ID */
template<int dim>
unsigned int 
ChemicalInterface<dim>::getID() const
{
	return chemical_id;	
}

}}// close namespace
#endif // chemical_interface.h