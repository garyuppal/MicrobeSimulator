#pragma once

#include "../utility/parameter_handler.h"
#include "./aging_chemical_interface.h"
#include "./aging_fdm_chemical.h"

#include <memory>


namespace MicrobeSimulator{ namespace Aging{

template<int dim>
class Chemicals{
public:
	Chemicals();

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm, double time_step);

	// read only access:
	const Aging::AgingChemicalInterface<dim>& 
	operator[](unsigned int i) const;

	// project function:
	void project_function(const Function<dim>& initial_condition, unsigned int i);

	void update(const std::vector<Point<dim> >& source_locations, 
						const std::vector<double>& sources,
						const std::vector<Point<dim> >& sink_locations, 
						const std::vector<double>& sinks);

	unsigned int getNumberChemicals() const;

	void output(std::string output_directory, 
		unsigned int save_step_number) const; 
	void output(std::string output_directory, 
		unsigned int save_step_number, unsigned int run_number) const;

	void printInfo(std::ostream& out) const;
private:
	unsigned int number_chemicals;
	std::vector<std::shared_ptr<AgingChemicalInterface<dim> > > 	chemicals; 
};

// IMPL
// -------------------------------------------------------------------
template<int dim>
Chemicals<dim>::Chemicals()
	:
	number_chemicals(0)
{}

template<int dim>
void 
Chemicals<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Chemicals");
		prm.declare_entry("Number chemicals", "1", Patterns::Unsigned());
		prm.declare_entry("Diffusion","{1}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Decay rate","{0}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Save chemicals", "False",Patterns::Bool());
		prm.declare_entry("FDM discretization", "{0.2}", Patterns::List(Patterns::Double()));
		prm.declare_entry("Boundary conditions",
			          "{WRAP,WRAP}",
			          Patterns::List(Patterns::Selection("WRAP|REFLECT|OPEN")),
			          "Boundary conditions.");
		prm.declare_entry("Flow rate", "0", Patterns::Double());
		prm.declare_entry("Bounds", "{{0,0}, {0,0}}",
			Patterns::List(Patterns::List(Patterns::Double())) );
	prm.leave_subsection();
}

template<int dim>
void
Chemicals<dim>::init(const ParameterHandler& prm, double time_step)
{
	const std::string section = "Chemicals";
	number_chemicals = prm.get_unsigned(section, "Number chemicals");

	for(unsigned int i = 0; i < number_chemicals; ++i)
	{
		chemicals.emplace_back(
			std::make_shared<AgingChemical<dim> >(prm, time_step, i));
	}
}

template<int dim>
const Aging::AgingChemicalInterface<dim>& 
Chemicals<dim>::operator[](unsigned int i) const
{
	return *(chemicals[i]);
}

template<int dim>
void 
Chemicals<dim>::project_function(const Function<dim>& initial_condition, unsigned int i)
{
	chemicals[i]->project_function(initial_condition);
}

template<int dim>
void 
Chemicals<dim>::update(const std::vector<Point<dim> >& source_locations, 
					const std::vector<double>& sources,
					const std::vector<Point<dim> >& sink_locations, 
					const std::vector<double>& sinks)
{
	/** same for now
	* @todo generalize to mutliple chemicsls 
	*/
	for(unsigned int i = 0; i < chemicals.size(); ++i)
		chemicals[i]->update(source_locations, sources[i], sink_locations, sinks[i]);
}

template<int dim>
unsigned int 
Chemicals<dim>::getNumberChemicals() const
{
	return chemicals.size();
}

template<int dim>
void 
Chemicals<dim>::output(std::string output_directory, 
	unsigned int save_step_number) const
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)
	{
		std::string outfile = output_directory
							+ "/Chemical_"
							+ dealii::Utilities::int_to_string(i, 2)
							+ "_" 
							+ dealii::Utilities::int_to_string(save_step_number,4)
							+ ".vtk"; // might want to change format dynamically ...
								// can set based on impl given by prm...
		std::ofstream out(outfile);

		chemicals[i]->print(out); 
	}
}

template<int dim>
void 
Chemicals<dim>::output(std::string output_directory, 
	unsigned int save_step_number, unsigned int run_number) const
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)
	{
		std::string outfile = output_directory
							+ "/Chemical_"
							+ dealii::Utilities::int_to_string(i, 2)
							+ "_R" 
							+ dealii::Utilities::int_to_string(run_number,4)
							+ "_" 
							+ dealii::Utilities::int_to_string(save_step_number,4)
							+ ".dat"; // might want to change format dynamically ...
								// can set based on impl given by prm...
		std::ofstream out(outfile);

		chemicals[i]->print(out); 
	}
}

template<int dim>
void
Chemicals<dim>::printInfo(std::ostream& out) const
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)
	{
		out << "Chemical " << i << ":" << std::endl;
		chemicals[i]->printInfo(out);
	}
}

}} // CLOSE NAMESPACE
/* chemicals.h */