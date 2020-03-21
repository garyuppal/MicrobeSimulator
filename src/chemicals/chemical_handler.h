#ifndef MICROBE_SIMULATOR_CHEMICAL_HANDLER_H
#define MICROBE_SIMULATOR_CHEMICAL_HANDLER_H

#include <deal.II/base/timer.h>

// direct solver:
#include <deal.II/lac/sparse_direct.h>

#include "../utility/cell_iterator_map.h" // switch to field class ... or actual map

#include "./chemical_fe_base.h"
#include "./chemical_interface.h"
#include "./fe_chemical.h"

#include "./control_functions.h"
#include "../utility/fe_tools.h"

#include "../utility/parameter_handler.h"

#include <memory>

namespace MicrobeSimulator{ namespace Chemicals{

template<int dim>
class ChemicalHandler{
public:
	ChemicalHandler();
	ChemicalHandler(unsigned int numchem);

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm,
			const Geometry<dim>& geometry, 
			const Triangulation<dim>& tria,
			const Velocity::AdvectionHandler<dim>& velocity_function,
			double chemical_time_step);

	// ~ChemicalHandler();

	// std::shared_ptr<const ChemicalInterface<dim> > getTestPointer() const;
	// do we still want to keep fe_base ? , could leave unitialized
	
	// std::vector<std::shared_ptr<const ChemicalInterface<dim> > >
	// getChemicalPointers() const;

	// update:
	void update(); 

	void update(const std::vector<Point<dim> >& locations,
				const std::vector<std::vector<double> >& amounts);

	void update(const std::vector<Point<dim> >& locations,
				const std::vector<std::vector<double> >& amounts,
				const Controls<dim>& control_functions);

	void project_function(unsigned int c_id, const Function<dim>& initial_condition);

	// read only access: (maybe allow write access too...)
	const ChemicalInterface<dim>& 
	operator[](unsigned int i) const;

	// double value(unsigned int chem_id) const; 
	unsigned int getNumberChemicals() const;

	void output(std::string output_directory, 
		unsigned int save_step_number) const; 

	void output_grid(std::string output_directory, 
		unsigned int save_step_number,
		const Geometry<dim>& geo) const; 

	void printInfo(std::ostream& out) const;
private:
	// shared pointer to base 
	std::shared_ptr<Chemical_FE_Base<dim> >					fe_base; // only create if needed

	std::shared_ptr<FETools::PointCellMap<dim> >			point_cell_map;

	std::vector<std::shared_ptr<ChemicalInterface<dim> > > 	chemicals; 
	// later make so can create fdm and fem chemicals

	// along side vector of control types...
	// std::vector<std::shared_ptr<TimedFunction<dim> > > 		control_functions; 
};

// IMPL
// -------------------------------------------------------------
template<int dim>
ChemicalHandler<dim>::ChemicalHandler()
{}

template<int dim>
ChemicalHandler<dim>::ChemicalHandler(unsigned int numchem)
{
	chemicals.reserve(numchem);
	for(unsigned int i = 0; i < numchem; ++i)
		chemicals.emplace_back( std::make_shared<FE_Chemical<dim> >() );
}

template<int dim>
void 
ChemicalHandler<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Chemicals");
		prm.declare_entry("Number chemicals","2",Patterns::Unsigned());
		prm.declare_entry("Implementation","{}",
						Patterns::List(Patterns::Selection("FE|DG|FDM")),
						"Implementation for chemicals. ...");
		prm.declare_entry("Diffusion",
							"{5,15}",
							Patterns::List(Patterns::Double()));
		prm.declare_entry("Decay rate",
							"{50,15}",
							Patterns::List(Patterns::Double()));
		prm.declare_entry("Save chemicals", "False",Patterns::Bool());
		prm.declare_entry("Grid save","False",Patterns::Bool());
		prm.declare_entry("Time step factor", "1", Patterns::Double());
		prm.declare_entry("Viscosity beta","1",Patterns::Double());
	prm.leave_subsection();
}

template<int dim>
void 
ChemicalHandler<dim>::init(const ParameterHandler& prm,
			const Geometry<dim>& geometry, 
			const Triangulation<dim>& tria,
			const Velocity::AdvectionHandler<dim>& velocity_function,
			double chemical_time_step)
{
	const std::string section = "Chemicals";
	bool using_base = false;

	const unsigned int numchem = prm.get_unsigned(section, "Number chemicals");
	// parameter handler should check dimensions ...
	// std::vector<std::string> impl = prm.get_string_vector(section, "Implementation");
	const std::vector<double> diff = prm.get_double_vector(section,"Diffusion");
	const std::vector<double> decay = prm.get_double_vector(section, "Decay rate");
	const double vbeta = prm.get_double(section, "Viscosity beta");
	chemicals.reserve(numchem);
	for(unsigned int i = 0; i < numchem; ++i)
	{
		// for now using FE chemicals, otherwise check impl
		if(using_base == false)
		{
			fe_base = std::make_shared<Chemical_FE_Base<dim> >(tria); // default degree of elements = 1
			
			std::cout << "calling setup..." << std::endl;
			fe_base->setup(velocity_function, geometry.getBoundaryConditions() ); // also add bcs !!!
	
			point_cell_map = std::make_shared<FETools::PointCellMap<dim> >();
			point_cell_map->initialize(geometry, fe_base->get_dof_handler(), 0.1 /* source resolution */);
			// can do initialization in constructor		

			fe_base->attach_point_cell_map(point_cell_map); // probably want to use shared pointer in fe base class
			using_base = true;
		}

		chemicals.emplace_back(std::make_shared<FE_Chemical<dim> >(
				fe_base, diff[i], decay[i], vbeta, chemical_time_step) );
	}

	// SquarePulse(double a, double on, double off, double del)

}

// --------------------------------------------------------------------
template<int dim>
void
ChemicalHandler<dim>::update()
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)	
		chemicals[i]->update();
}

template<int dim>
void 
ChemicalHandler<dim>::update(const std::vector<Point<dim> >& locations,
				const std::vector<std::vector<double> >& amounts)
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)
		chemicals[i]->update(locations, amounts[i]); 
}

template<int dim>
void 
ChemicalHandler<dim>::update(const std::vector<Point<dim> >& locations,
				const std::vector<std::vector<double> >& amounts,
				const Controls<dim>& control_functions)
{
	for(unsigned int i = 0; i < chemicals.size(); ++i)
		chemicals[i]->update(locations, amounts[i], control_functions[i]); 
}

template<int dim>
void 
ChemicalHandler<dim>::project_function(unsigned int c_id, 
	const Function<dim>& initial_condition)
{
	chemicals[c_id]->project_function(initial_condition);
}

///--------------------------------------------------------------------

template<int dim>
unsigned int 
ChemicalHandler<dim>::getNumberChemicals() const
{
	return chemicals.size();
}

template<int dim>
const ChemicalInterface<dim>&
ChemicalHandler<dim>::operator[](unsigned int i) const
{
	return *(chemicals[i]);
}

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
		std::ofstream out(outfile);

		chemicals[i]->print(out, i); 
	}
}

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