#ifndef MICROBESIMULATOR_REFACTORED_CONTROL_FUNCTIONS_H
#define MICROBESIMULATOR_REFACTORED_CONTROL_FUNCTIONS_H

#include <deal.II/base/function.h>
#include "../utility/utility.h"
#include "../utility/parameter_handler.h"

namespace MicrobeSimulator{ namespace RefactoredChemicals{

// ---------------------------------------------------------------------
// BASE
// ---------------------------------------------------------------------
template<int dim>
class TimedFunction : public Function<dim>{
public:
	TimedFunction() : Function<dim>(), function_time(0) {}
	virtual ~TimedFunction() {}
	virtual void update_time(double dt) {function_time += dt;}
	virtual void printInfo(std::ostream&) const {}
protected:
	double function_time;
};	

// ---------------------------------------------------------------------
// SQUARE PULSE
// ---------------------------------------------------------------------
template<int dim>
class SquarePulse : public TimedFunction<dim>{
public:
	SquarePulse(double a, double on, double off, double del);
	// void update_time(); // if using this then need to be virtual and part of interface...
	double value(const Point<dim>& p, 
		const unsigned int component = 0) const override;

	void printInfo(std::ostream& out) const override;
private:
	double amplitude;
	double on_period;
	double off_period;
	double delay;

	// double time; // internal variable updated to get value, should update before calling value

	bool functionIsOn() const;
};

// IMPL
// ---------------------------------------------------------------------
template<int dim>
SquarePulse<dim>::SquarePulse(double a, double on, double off, double del)
	:
	TimedFunction<dim>(),
	amplitude(a),
	on_period(on),
	off_period(off),
	delay(del)
{}

template<int dim>
double 
SquarePulse<dim>::value(const Point<dim>& /*p*/, 
		const unsigned int /* component */) const 
{
	if( functionIsOn() )
		return amplitude;

	return 0;
}

template<int dim>
bool 
SquarePulse<dim>::functionIsOn() const
{
	const double t = (this->function_time) - delay;

	if(t > 0)
	{
		const double remainder = fmod(t, on_period + off_period);
		if(remainder <= on_period)
			return true;
	}

	return false;
}

template<int dim>
void 
SquarePulse<dim>::printInfo(std::ostream& out) const 
{
	out << "\n\n" << Utility::medium_line << std::endl
		<< "\t\t SQUARE PULSE FUNCTION:" << std::endl
		<<  Utility::medium_line << std::endl
   	  	<< "Height: " << amplitude << std::endl
   		<< "On period: " << on_period << std::endl
   		<< "Off period: " << off_period << std::endl
   		<< "Delay: " << delay << std::endl
		<< std::endl << Utility::medium_line 
		<< std::endl << std::endl << std::endl;
}

// ---------------------------------------------------------------------
// SINUSOIDAL
// ---------------------------------------------------------------------

// template<int dim>
// class Sinusoidal : public Function<dim>{
// public:
// 	Sinusoidal() : Function<dim>() {}
// 	double value(const Point<dim>& p, 
// 		const unsigned int component = 0) const override;
// };

// template<int dim>
// double 
// Sinusoidal<dim>::value(const Point<dim>& /* p */, 
// 	const unsigned int /* component */) const 
// {
// 	return 0;
// }







// --------------------------------------------------------------------------------------------
// CONTROL FUNCTION HANDLER CLASS:
// --------------------------------------------------------------------------------------------

/** \brief Class to handle control functions for chemical fields
*/
/**
* 
* @todo Implement more than just the square pulse maybe
*/
template<int dim>
class Controls{
public:
	Controls();

	static void declare_parameters(ParameterHandler& prm);

	void setup(const ParameterHandler& prm);

	const std::string& getControlType(unsigned int i) const;
	const TimedFunction<dim>& operator[](unsigned int i) const;

	void update_time(double dt);

	bool isActive() const;

	void printInfo(std::ostream& out) const;

private:
	bool active;
	std::vector<std::string> 							control_types;
	std::vector<std::shared_ptr<TimedFunction<dim> > > 	control_functions; 
};

// ---------------------------------------------------------------------------------------------
/** \brief Constructor
*/
template<int dim>
Controls<dim>::Controls()
	:
	active(false)
{}

/** \brief Sets up parameters to be read in for Controls class
*/
template<int dim>
void
Controls<dim>::declare_parameters(ParameterHandler& prm) 
{
	prm.enter_subsection("Controls");
		prm.declare_entry("Type","{None,None}",
							Patterns::List(Patterns::Selection("None|Square pulse")));
		
		prm.enter_subsection("Square pulse"); // currently fixed for 2 chemicals?
			prm.declare_entry("Height","{0,0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("On period","{0,0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Off period", "{0,0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Delay", "{0,0}",Patterns::List(Patterns::Double()));
		prm.leave_subsection();
	prm.leave_subsection();
}

/** \brief Returns string reference to control function type for ith control
*/
template<int dim>
const std::string& 
Controls<dim>::getControlType(unsigned int i) const
{
	return control_types[i];
}

/** \brief Read only access operator for a "TimedFunction<dim>" 
*/
template<int dim>
const TimedFunction<dim>& 
Controls<dim>::operator[](unsigned int i) const
{
	return *(control_functions[i]);
}

/** \brief Sets up control functions from parameters
*/
template<int dim>
void
Controls<dim>::setup(const ParameterHandler& prm)
{
	// get total possible number from number of chemicals	
	// const unsigned int numchem = prm.get_unsigned("Chemicals", "Number chemicals");
	const std::string section = "Controls";

	/** currently up to user to make sure same number of controls as chemicals */
	control_types = prm.get_list(section, "Type");

	for(unsigned int i = 0; i < control_types.size(); ++i)
	{
		// assuming either both are None or both square***
		if( boost::iequals(control_types[i], "Square pulse") )
		{
			active = true; /** @todo Instead of using active, maybe allow for None type
							* or have array of bools ... */
			std::cout << "Using square pulse control function for chemical " << i << std::endl;
			const std::string subsection = section + ".Square pulse";
			const double a = prm.get_double_vector(subsection, "Height")[i];
			const double on = prm.get_double_vector(subsection, "On period")[i];
			const double off = prm.get_double_vector(subsection, "Off period")[i];
			const double delay = prm.get_double_vector(subsection, "Delay")[i];

			control_functions.emplace_back(std::make_shared<SquarePulse<dim> >(a,on,off,delay));
		}
	}
}

/** \brief Increments internal clock for control functions by time step dt
*/
template<int dim>
void 
Controls<dim>::update_time(double dt) 
{
	for(unsigned int i = 0; i < control_functions.size(); ++i)
		control_functions[i]->update_time(dt);	
}

/** \brief Checks if control functions are other than "NONE"
*/
template<int dim>
bool 
Controls<dim>::isActive() const
{
	return active;
}

/** \brief Prints control function information to ostream taken as input
*/
template<int dim>
void
Controls<dim>::printInfo(std::ostream& out) const
{
	/** Takes an ostream reference as input. For example, to print to console
	* just call function with std::cout as the parameter.  
	*/
	for(unsigned int i = 0; i < control_functions.size(); ++i)
		control_functions[i]->printInfo(out);
}


}} // close namespace
#endif