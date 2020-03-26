#ifndef MICROBESIMULATOR_FULL_SIMULATOR_H
#define MICROBESIMULATOR_FULL_SIMULATOR_H

#include <deal.II/grid/tria.h>
using dealii::Triangulation;

#include "../bacteria/bacteria_handler.h" 
#include "../bacteria/bacteria_fitness.h"
#include "../chemicals/chemical_handler.h"
#include "../advection/advection_handler.h"
#include "../geometry/geometry.h"
#include "../geometry/geometry_builder.h" 
#include "../utility/parameter_handler.h"
#include "../utility/command_line_parser.h"

#include <array>
#include <algorithm>
#include <deal.II/base/function.h>

/** @file
* \brief Simulator file
* Simulation class defined here
*/

/** \brief Main namespace for MicrobeSimulator */
namespace MicrobeSimulator{

	namespace SimulationTools{
		/** \brief Output vector to file */
		void output_vector(const std::vector<double>& vec, 
			const std::string& data_name,
			const std::string& output_directory)
		{
			std::string outFile = output_directory + "/" + data_name + ".dat";
			std::ofstream out(outFile);

			for(unsigned int i = 0; i < vec.size(); ++i)	
				out << vec[i] << std::endl;
		}
	} // close SimulationTools namespace


/** \brief Simulator class
* This class handles reading in parameters, constructing, and executing the simulation.
*/
template<int dim>
class FullSimulator{
public:
	FullSimulator(const CommandLineParameters& cmd_prm);

	void run();
private:
	ParameterHandler prm;

	Triangulation<dim>						triangulation;
	Geometry<dim> 							geometry;
	Velocity::AdvectionHandler<dim>			velocity_function;
	Chemicals::ChemicalHandler<dim> 		chemicals;

	Chemicals::Controls<dim>				control_functions;

	Bacteria::BacteriaHandler<dim>			bacteria;
	Bacteria::OR_Fitness<dim>				fitness_function;
		/** @todo type of fitness given by parameter as well */

	// SYSTEM CONSTANTS:
	std::string 							output_directory;
	double 									run_time;
	double 									time;
	double 									chemical_time_step;
	double 									bacteria_time_step;
	double 									save_period;
	unsigned int 							save_step_number;
	unsigned int 							time_step_number;
	unsigned int 							bacteria_time_step_multiplier;
	unsigned int 							number_reintroduced; // for adding in new groups

	// output:
	void output_bacteria() const;
	void output_chemicals(bool isGridSave) const;

	// update:
	void update_bacteria();
	void reintro_bacteria(); 

	// setup:
	void declare_parameters();
	void setup_system(); 
	void assign_local_parameters();
	void setup_time_steps();
	double get_chemical_time_step();
	void setup_fitness();

	// Simulators:
	void run_microbes();
};

// IMPL
// -----------------------------------------------------------------

/** \brief Construct simulator */
template<int dim>
FullSimulator<dim>::FullSimulator(const CommandLineParameters& cmd_prm)
	:
	prm(cmd_prm.getParameterFile(), cmd_prm.getJobID()),
	triangulation(Triangulation<dim>::maximum_smoothing),
	output_directory(cmd_prm.getOutputDirectory()),
	run_time(0),
	time(0),
	chemical_time_step(1),
	bacteria_time_step(1),
	save_step_number(0),
	time_step_number(0),
	bacteria_time_step_multiplier(1),
	number_reintroduced(0)
{}

/** \brief Run simulation */
template<int dim>
void
FullSimulator<dim>::run()
{
	setup_system();
	run_microbes();
}

/** \brief Run evolving microbe simulation */
template<int dim>
void
FullSimulator<dim>::run_microbes()
{
	std::cout << std::endl << std::endl
		<< "Starting microbe simulation" << std::endl
		<< Utility::long_line << std::endl
		<< Utility::long_line << std::endl << std::endl;

	// save period:
	const unsigned int modsave
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );
	const bool isSavingChemicals = prm.get_bool("Chemicals","Save chemicals");
	const bool isGridSave = prm.get_bool("Chemicals", "Grid save");

	// spread out initial bacteria:
	const unsigned int intial_spread = 15;
	for(unsigned int i = 0; i < intial_spread; ++i)
		bacteria.move(chemical_time_step, geometry, velocity_function);

	// check after intial spread:
	output_bacteria();
	++save_step_number;

	const bool reintro = prm.get_bool("Bacteria","Reintroducing");
	bool dont_kill = true;

	do{
		if(reintro)
			reintro_bacteria();

		// update time:
		time += chemical_time_step;
		++time_step_number;

		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;
			if(isSavingChemicals)
				output_chemicals(isGridSave);
			output_bacteria();
			++save_step_number;
		}

		/** @todo give option to have some controls be "none", or add constant zero type
		*/
		if(control_functions.isActive())
		{
			control_functions.update_time(chemical_time_step); // increment internal clock
			chemicals.update(bacteria.getAllLocations(), bacteria.getAllRates(), control_functions);
		}
		else
		{
			chemicals.update(bacteria.getAllLocations(), bacteria.getAllRates());
		}

		if(time_step_number % bacteria_time_step_multiplier == 0)
			update_bacteria();

		if(!reintro)
		{
		   	if( !bacteria.isAlive() )
		   	{
		   		dont_kill = false;
		   		std::cout << "\n\nEverybody died!" << std::endl;
		   	}
	    }
	}while( (time < run_time) && dont_kill );

	if(reintro)
		SimulationTools::output_vector(bacteria.get_pg_rates(),
			"pg_rates",output_directory);
} // run_microbes()

/** \brief Method to reintoduce new microbe groups into simulation */
template<int dim>
void
FullSimulator<dim>::reintro_bacteria()
{
	const unsigned int intro_number = 
		std::floor(time/(prm.get_double("Bacteria","Reintroduction period")));

	if(intro_number > number_reintroduced)
	{
		bacteria.reintro(prm, geometry);
		++number_reintroduced;
	}	
}


// UPDATE:
// ---------------------------------------------------

/** \brief Update bacteria */
/** Do random walk, reproduce, and mutate */
template<int dim>
void
FullSimulator<dim>::update_bacteria()
{
	bacteria.move(bacteria_time_step, geometry, velocity_function);
	bacteria.reproduce(bacteria_time_step, fitness_function);
	bacteria.mutate(bacteria_time_step);
}

// OUTPUT:
// ------------------------------------------------------

/** \brief Output bacteria to file in simulation output directory */
template<int dim>
void
FullSimulator<dim>::output_bacteria() const
{
	std::string outfile = output_directory
						+ "/bacteria_"
						+ dealii::Utilities::int_to_string(save_step_number,4)
						+ ".dat";
	std::ofstream out(outfile);
	bacteria.print(out);
}

/** \brief Output chemical fields to file in simulation output directory */
template<int dim>
void
FullSimulator<dim>::output_chemicals(bool isGridSave) const
{
	if(isGridSave)
		chemicals.output_grid(output_directory, save_step_number, geometry);
	else
		chemicals.output(output_directory, save_step_number);
}


// SETUP:
// ------------------------------------------------------------------------

/** \brief Setup simulation system */
/** Read in parameters, construct and intialize required objects,
* and output simulation information to simulation directory 
*/
template<int dim>
void
FullSimulator<dim>::setup_system()
{
	std::cout << "\n\nSETTING UP SYSTEM\n" << std::endl;

	std::cout << "...setting up parameters" << std::endl;

	declare_parameters(); 
	prm.parse_parameter_file();
	prm.print_simple(std::cout);
	std::ofstream out(output_directory + "/parameters.dat");
	prm.print(out);
	prm.printLoopedParameterGrid(std::cout);
	assign_local_parameters();

	std::cout << "...setting up geometry" << std::endl;
		GeometryTools::GeometryBuilder<dim> geo_bldr(prm);
		geo_bldr.printInfo(std::cout);
		geo_bldr.build_geometry(geometry);
		geometry.printInfo(std::cout);
		geometry.outputGeometry(output_directory);

	std::cout << "...setting up grid" << std::endl;
		geo_bldr.build_grid(geometry, triangulation);
		GridGenerationTools::output_grid(output_directory,"before_stokes_grid",triangulation);

	std::cout << "...setting up velocity" << std::endl;
		velocity_function.init(prm,geometry,triangulation,output_directory);
		GridGenerationTools::output_grid(output_directory,"after_stokes_grid",triangulation);

	std::cout << "...setting up time steps" << std::endl;
		setup_time_steps();
		std::cout << "\t...using chemical time step: " << chemical_time_step << std::endl;
		std::cout << "\t...using bacteria time step " << bacteria_time_step << std::endl;

	std::cout << "...setting up chemicals" << std::endl;
		chemicals.init(prm, geometry,
			triangulation, velocity_function, chemical_time_step);
		chemicals.printInfo(std::cout);

	std::cout << "...setting up control functions" << std::endl;
		control_functions.setup(prm);
		control_functions.printInfo(std::cout);

	std::cout << "...setting up bacteria" << std::endl;
		bacteria.init(prm, geometry);
		bacteria.printInfo(std::cout);

	std::cout << "...setting up fitness" << std::endl;
		setup_fitness();
		fitness_function.printInfo(std::cout);
	std::cout << std::endl << std::endl;
} // setup_system()

/** \brief Store local parameters read in from file */
template<int dim>
void
FullSimulator<dim>::assign_local_parameters()
{
	run_time = prm.get_double("Run time");
	save_period = prm.get_double("Save period");
}

/** \brief Setup suitable time steps for bacteria and chemicals */
template<int dim>
void
FullSimulator<dim>::setup_time_steps()
{
	chemical_time_step = get_chemical_time_step();
	const double bacteria_max_time_step = 0.01;
	bacteria_time_step_multiplier = std::max(
			static_cast<int>( std::floor(bacteria_max_time_step / chemical_time_step) ) , 1);
	bacteria_time_step = bacteria_time_step_multiplier*chemical_time_step;
}

/** \brief Get stable chemical time step from CFL condition */
template<int dim>
double
FullSimulator<dim>::get_chemical_time_step()
{
	const double min_time_step = 0.001;
	double maximal_velocity = velocity_function.get_maximum_velocity(0); 
	maximal_velocity = std::max(maximal_velocity, 0.1);

	double cfl_time_step = dealii::GridTools::minimal_cell_diameter(triangulation)
							/ maximal_velocity;
	std::cout << "\tCFL_TIME_STEP: " << cfl_time_step << std::endl;

	cfl_time_step = std::min( min_time_step,
		(prm.get_double("Chemicals", "Time step factor"))*cfl_time_step );

	return cfl_time_step;
}

/** \brief Setup fitness object for bacteria reproduction */
/** @todo Move this to a fitness_handler class ... */
template<int dim>
void
FullSimulator<dim>::setup_fitness()
{
	const std::string section = "Fitness";
	fitness_function.attach_chemicals(chemicals);

	// this needs to be a base function if using polymorphism: (can also use constructor)
	double sc = prm.get_double(section, "Secretion cost");
	std::cout << "secretion cost obtained: " << sc << std::endl;

	fitness_function.setup_fitness_constants(
			prm.get_double_vector(section, "Chemical fitness"),
			prm.get_double(section, "Secretion cost"),
			prm.get_double_vector(section, "Chemical saturation"));
}

/** \brief Declare local parameters needed for simulator 
* to be read in from file
*/
template<int dim>
void 
FullSimulator<dim>::declare_parameters()
{
	prm.declare_entry("Simulator type",
						"Bacteria",
						Patterns::Selection("Bacteria|IC_AS|Aging|Chemotaxis"),
						"Simulator type.");
	prm.declare_entry("Time step", "1", Patterns::Double());
	prm.declare_entry("Run time","0",Patterns::Double());
	prm.declare_entry("Save period","1",Patterns::Double());
	prm.declare_entry("Output directory","./",Patterns::Anything());

	// declare class based parameters:
	Velocity::AdvectionHandler<dim>::declare_parameters(prm);
	GeometryTools::GeometryBuilder<dim>::declare_parameters(prm);
	Chemicals::ChemicalHandler<dim>::declare_parameters(prm);
	Chemicals::Controls<dim>::declare_parameters(prm);
	Bacteria::BacteriaHandler<dim>::declare_parameters(prm);
	Bacteria::Fitness::declare_parameters(prm);
}


} // close namespace
#endif
