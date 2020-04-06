#ifndef MICROBE_SIMULATOR_REFACTORED_CHEMICAL_HANDLER_H
#define MICROBE_SIMULATOR_REFACTORED_CHEMICAL_HANDLER_H

// #include <deal.II/base/timer.h>

// // direct solver:
// #include <deal.II/lac/sparse_direct.h>

// #include "../utility/cell_iterator_map.h" // switch to field class ... or actual map

// #include "./chemical_fe_base.h"

// #include "./control_functions.h"
// #include "../utility/fe_tools.h"

#include "../utility/parameter_handler.h"

#include <memory>
#include "./chemical_interface.h"

// chemical implementations:
// #include "./field_base.h"
#include "./fdm_chemical.h"
#include "./fe_chemical.h"

// finite elements:
// #include "./chemical_fe_base.h" // cg_support and dg_support classes

// TODO, parameter setup and initialization (from simulator)
// store impl type to change saving method...
// initializing time steps and calculating stability....

namespace MicrobeSimulator{ namespace RefactoredChemicals{

template<int dim>
class ChemicalHandler{
public:
	ChemicalHandler();

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm,
			const Geometry<dim>& geo, 
			const Triangulation<dim>& tria,
			const Velocity::AdvectionHandler<dim>& velocity_function,
			double time_step);

	// update: (apply to all), or give index???, 
	// and only supply vector for ith chemical amounts, or ith control function
	void update(); 

	void update(const std::vector<Point<dim> >& locations,
				const std::vector<std::vector<double> >& amounts);

	void update(const std::vector<Point<dim> >& locations,
				const std::vector<std::vector<double> >& amounts,
				const Controls<dim>& control_functions);

	// project function:
	void project_function(const Function<dim>& initial_condition, unsigned int i);

	// read only access:
	const ChemicalInterface<dim>& 
	operator[](unsigned int i) const;

	unsigned int getNumberChemicals() const;

	void output(std::string output_directory, 
		unsigned int save_step_number) const; 

	void output_grid(std::string output_directory, 
		unsigned int save_step_number,
		const Geometry<dim>& geo) const;

	void printInfo(std::ostream& out) const;

//***************************************************
	// for debugging:
	// void test_convergence(const ParameterHandler& prm); 

private:
	// shared pointer to base 
	// std::shared_ptr<Chemical_FE_Base<dim> >					fe_base; // only create if needed
	bool baseInit; // may be good idea to store what chemicals were intialized and what type, 
		// may or may not need separate dg_base (though we'll probably never use this)
	// but then, instead of baseinit, just check if cg or dg type chem has already been declared

	std::shared_ptr<Chemical_FE_Base<dim> >					cg_base; // needs triangulation to initialize
	std::vector<std::shared_ptr<ChemicalInterface<dim> > > 	chemicals; 
};

// IMPL
// -------------------------------------------------------------
/** \brief Constructor for chemical handler */
template<int dim>
ChemicalHandler<dim>::ChemicalHandler()
	:
	baseInit(false)
{}

/** \brief Declare parameters needed for chemicals */
template<int dim>
void 
ChemicalHandler<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Chemicals");
		prm.declare_entry("Number chemicals","2",Patterns::Unsigned());
		prm.declare_entry("Implementation","{FE, FE}",
						Patterns::List(Patterns::Selection("FE|DG|FDM")),
						"Implementation method for chemicals.");
		prm.declare_entry("Diffusion",
							"{5,15}",
							Patterns::List(Patterns::Double()));
		prm.declare_entry("Decay rate",
							"{50,15}",
							Patterns::List(Patterns::Double()));
		prm.declare_entry("Save chemicals", "False",Patterns::Bool());
		prm.declare_entry("Save type","VTK", Patterns::Selection("VTK|Grid|Both"));
		prm.declare_entry("Time step factor", "1", Patterns::Double());
		prm.declare_entry("Viscosity beta","1",Patterns::Double());

		prm.declare_entry("FDM discretization","{0.2, 0.2}",
						Patterns::List(Patterns::Double()) ); // figure out later
		// prm cell size (*** can put all needed parameters in arrays)
		// may be confusing if using mixed types ... but shouldn't really be doing that anyway..
	prm.leave_subsection();
}

/** \brief Initialize chemicals from system parameters */
template<int dim>
void 
ChemicalHandler<dim>::init(const ParameterHandler& prm,
							const Geometry<dim>& geo, 
							const Triangulation<dim>& tria,
							const Velocity::AdvectionHandler<dim>& velocity_function,
							double time_step)
{
	std::vector<std::string> 
		impl_types = prm.get_string_vector("Chemicals", "Implementation");

	for(unsigned int i = 0; i < impl_types.size(); ++i)
	{
		if( boost::iequals(impl_types[i], "FDM") )
		{
			chemicals.emplace_back( 
				std::make_shared<FDMChemical<dim> >(prm, geo, velocity_function, time_step, i) );
		}
		else if( boost::iequals(impl_types[i], "FE") )
		{
			// make sure to only initialize once
			if(!baseInit)
			{
				cg_base = std::make_shared<Chemical_FE_Base<dim> >(tria); // need to include velocity function? // could also init i suppose
				cg_base->setup(velocity_function, geo.getBoundaryConditions()); // move to constructor ...
				baseInit = true; 
			}
			chemicals.emplace_back(
				std::make_shared<FE_Chemical<dim> >(prm, time_step, i, cg_base) );
		}
		else if( boost::iequals(impl_types[i], "DG") )
		{
			// make a dg base
			// fe_base = std::make_shared<CG_Support<dim> >(prm, tria, velocity_function);
			throw std::runtime_error("still need to implement dg chems");
		}
		else
		{
			throw std::runtime_error("invalid chemical implementation");
		}
	}
}

/** \brief Return const reference to ith chemical for read access */
template<int dim>
const ChemicalInterface<dim>&
ChemicalHandler<dim>::operator[](unsigned int i) const
{
	return *(chemicals[i]);
}

/** \brief Update all chemicals by their time step */
template<int dim>
void
ChemicalHandler<dim>::update()
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)	
		chemicals[i]->update();
}

/** \brief Update all chemicals by their time step 
* adding sources of given amplitudes at given locations */
template<int dim>
void 
ChemicalHandler<dim>::update(const std::vector<Point<dim> >& locations,
				const std::vector<std::vector<double> >& amounts)
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)
		chemicals[i]->update(locations, amounts[i]); 
}

/** \brief Update all chemicals by their time step 
* adding sources of given amplitudes at given locations 
* and adding control functions */
template<int dim>
void 
ChemicalHandler<dim>::update(const std::vector<Point<dim> >& locations,
				const std::vector<std::vector<double> >& amounts,
				const Controls<dim>& control_functions)
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)
		chemicals[i]->update(locations, amounts[i], control_functions[i]); 
}

/** \brief Project function to ith chemcical */
template<int dim>
void 
ChemicalHandler<dim>::project_function(const Function<dim>& initial_condition,
	unsigned int i)
{
	chemicals[i]->project_function(initial_condition);
}

/** \brief Get number of active chemicals */
template<int dim>
unsigned int 
ChemicalHandler<dim>::getNumberChemicals() const
{
	return chemicals.size();
}

/** \brief Output chemical to file */
template<int dim>
void 
ChemicalHandler<dim>::output(std::string output_directory, 
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

/** \brief Output grid of chemical values to file */
template<int dim>
void 
ChemicalHandler<dim>::output_grid(std::string output_directory, 
		unsigned int save_step_number,
		const Geometry<dim>& geo) const
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)
	{
		std::string outfile = output_directory
							+ "/Chemical_"
							+ dealii::Utilities::int_to_string(i, 2)
							+ "_" 
							+ dealii::Utilities::int_to_string(save_step_number,4)
							+ ".dat"; // might want to change format dynamically ...
		std::ofstream out(outfile);

		// get points
		std::vector<Point<dim> > querry_points = geo.getQuerryPoints(0.2); 
		// get from geometry. using 0.2 resolution as default
		std::vector<double> values = chemicals[i]->value_list(querry_points);

		// print to file:
		for(unsigned int q = 0; q < querry_points.size(); ++q)
			out << querry_points[q] << " " << values[q] << std::endl;
	}
}

/** \brief Print out information for each chemical */
template<int dim>
void 
ChemicalHandler<dim>::printInfo(std::ostream& out) const
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)
	{
		out << "Chemical " << i << ":" << std::endl;
		chemicals[i]->printInfo(out);
	}
}

}} // CLOSE NAMESPACE
#endif