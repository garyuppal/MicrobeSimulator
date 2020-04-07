#ifndef MICROBESIMULATOR_SIMULATOR_H
#define MICROBESIMULATOR_SIMULATOR_H

#include <deal.II/grid/tria.h>
using dealii::Triangulation;

#include "../bacteria/bacteria_handler.h" 
#include "../bacteria/bacteria_fitness.h"
#include "../advection/advection_handler.h"
#include "../geometry/geometry.h"
#include "../geometry/geometry_builder.h" 
#include "../utility/parameter_handler.h"
#include "../utility/command_line_parser.h"

#include "./test_functions.h"

// #include "../chemicals/chemical_handler.h"
#include "../refactored_chemicals/chemical_handler.h"

#include <array>
#include <algorithm>

#include <deal.II/base/function.h>
#include <deal.II/base/convergence_table.h>

/** \brief Main namespace for MicrobeSimulator */
namespace MicrobeSimulator{ 

/** \brief Namespace for simulator */
namespace Simulation{

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

/** \brief Simulator class
* This class handles reading in parameters, constructing, and executing the simulation.
*/
template<int dim>
class Simulator{
public:
	Simulator(const CommandLineParameters& cmd_prm);

	// if adding other, multiple types of simulation, can maybe have cmd_parser select which on,
	// or create simulator interface and separate classes...
	void run();
	void run_convergence_check();
private:
	ParameterHandler prm;

	Triangulation<dim>						triangulation;
	Geometry<dim> 							geometry;
	Velocity::AdvectionHandler<dim>			velocity_function;

	RefactoredChemicals::ChemicalHandler<dim> chemicals;
	RefactoredChemicals::Controls<dim>			control_functions;

	Bacteria::BacteriaHandler<dim>			bacteria;
	Bacteria::TestNewFitness::Fitness_Function<dim>	 fitness_function;

	// SYSTEM CONSTANTS:
	std::string 							output_directory;
	double 									run_time;
	double 									time;
	double 									chemical_time_step; // get from handler?
	double 									bacteria_time_step; // get from bacteria?
	double 									save_period;
	unsigned int 							save_step_number;
	unsigned int 							time_step_number;
	unsigned int 							bacteria_time_step_multiplier;
	unsigned int 							number_reintroduced; // for adding in new groups


	// output:
	void output_bacteria() const;
	void output_chemicals() const;
	void output_chemical_mass() const;
	void output_local_chemicals() const;


	// update:
	void update_bacteria();
	void reintro_bacteria(); 

	// setup:
	void declare_parameters();
		void setup_parameters();
		void setup_geometry_grid();
		void setup_velocity();
		void setup_chemicals();
		void setup_bacteria();

	void setup_system(); 
	void assign_local_parameters();
	void setup_time_steps();
	double get_chemical_time_step();
	void setup_fitness();

	// Simulators:
	void run_microbes();

	// FOR CONVERGENCE TEST:
	void reset_system(const unsigned int cycle, 
		std::vector<TestFunctions::GaussianSolution<dim> >& gaussian_functions);
	void solve_chemicals();
	void process_solution(const unsigned int cycle, 
		dealii::ConvergenceTable& convergence_table);
};

// IMPL
// -----------------------------------------------------------------

/** \brief Construct simulator */
template<int dim>
Simulator<dim>::Simulator(const CommandLineParameters& cmd_prm)
	:
	prm(cmd_prm.getParameterFile(), cmd_prm.getJobID()),
	triangulation(Triangulation<dim>::maximum_smoothing),
	// test_ff(chemicals, prm),
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
Simulator<dim>::run()
{
	setup_parameters();
	setup_system();
	run_microbes();
}

/** \brief Run evolving microbe simulation */
template<int dim>
void
Simulator<dim>::run_microbes()
{
	std::cout << std::endl << std::endl
		<< "Starting microbe simulation" << std::endl
		<< Utility::long_line << std::endl
		<< Utility::long_line << std::endl << std::endl;

	// save period:
	const unsigned int modsave
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );
	const bool isSavingChemicals = prm.get_bool("Chemicals","Save chemicals");
	const bool recordMass = prm.get_bool("Debug","Record chemical mass");
	const bool trackMicrobeChem = prm.get_bool("Debug","Track microbe chemicals");

	// const bool isGridSave = prm.get_bool("Chemicals", "Grid save");

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
				output_chemicals();
			if(recordMass)
				output_chemical_mass();
			if(trackMicrobeChem)
				output_local_chemicals();

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
		output_vector(bacteria.get_pg_rates(),
			"pg_rates",output_directory);
} // run_microbes()

/** \brief Method to reintoduce new microbe groups into simulation */
template<int dim>
void
Simulator<dim>::reintro_bacteria()
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
Simulator<dim>::update_bacteria()
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
Simulator<dim>::output_bacteria() const
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
Simulator<dim>::output_chemicals() const
{
	std::string save_type = prm.get_string("Chemicals","Save type");

	if( boost::iequals(save_type, "VTK") )
	{
		chemicals.output(output_directory, save_step_number);
	}
	else if( boost::iequals(save_type, "Grid"))
	{
		chemicals.output_grid(output_directory, save_step_number, geometry);
	}
	else if( boost::iequals(save_type, "Both") )
	{
		chemicals.output(output_directory, save_step_number);
		chemicals.output_grid(output_directory, save_step_number, geometry);
	}
}

/** \brief Output integral of chemical field over space 
* (primarily for debugging) */
template<int dim>
void
Simulator<dim>::output_chemical_mass() const
{
	std::string outfile = output_directory
					+ "/chemical_masses_"
					+ dealii::Utilities::int_to_string(save_step_number,4)
					+ ".dat";
	std::ofstream out(outfile);

	for(unsigned int i = 0; i < chemicals.getNumberChemicals(); ++i)
		out << chemicals[i].getMass() << std::endl;
}

/** \brief Output local chemical concentration for each microbe 
* (primarily for debugging) */
template<int dim>
void
Simulator<dim>::output_local_chemicals() const
{
	const unsigned int nc = chemicals.getNumberChemicals();

	std::vector< Point<dim> > locations = bacteria.getAllLocations();
	for(unsigned int c = 0; c < nc; ++c)
	{
		std::vector<double> local_chems( chemicals[c].value_list(locations) );

		std::string outFile = output_directory 
				+ "/local_chem" 
				+ dealii::Utilities::int_to_string(c, 2)
				+ "_"
				+ dealii::Utilities::int_to_string(save_step_number, 4)
				+ ".dat";
		std::ofstream out(outFile);

		for(unsigned int b = 0; b < local_chems.size(); ++b)
			out << locations[b] << " " << local_chems[b] << std::endl;
	}
}

// SETUP:
// ------------------------------------------------------------------------

/** \brief Setup system parameters and read in from file */
template<int dim>
void
Simulator<dim>::setup_parameters()
{
	std::cout << "...setting up parameters" << std::endl;

	declare_parameters(); 
	prm.parse_parameter_file();
	prm.print_simple(std::cout);
	std::ofstream out(output_directory + "/parameters.dat");
	prm.print(out);
	prm.printLoopedParameterGrid(std::cout);
	assign_local_parameters();
}

/** \brief setup and output geometry and mesh grid */
template<int dim>
void
Simulator<dim>::setup_geometry_grid()
{
	std::cout << "...setting up geometry" << std::endl;
		GeometryTools::GeometryBuilder<dim> geo_bldr(prm);
		geo_bldr.printInfo(std::cout);
		geo_bldr.build_geometry(geometry);
		geometry.printInfo(std::cout);
		geometry.outputGeometry(output_directory);

	std::cout << "...setting up grid" << std::endl;
		geo_bldr.build_grid(geometry, triangulation);
		GridGenerationTools::output_grid(output_directory,"before_stokes_grid",triangulation);
}

/** \brief Setup up velocity */
template<int dim>
void
Simulator<dim>::setup_velocity()
{
	std::cout << "...setting up velocity" << std::endl;
		velocity_function.init(prm,geometry,triangulation,output_directory);
		GridGenerationTools::output_grid(output_directory,"after_stokes_grid",triangulation);
}

/** \brief Setup up chemicals and control functions */
template<int dim>
void
Simulator<dim>::setup_chemicals()
{
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
}

/** \brief Setup up bacteria */
template<int dim>
void
Simulator<dim>::setup_bacteria()
{
	std::cout << "...setting up bacteria" << std::endl;
		bacteria.init(prm, geometry);
		bacteria.printInfo(std::cout);

	std::cout << "...setting up fitness" << std::endl;
		setup_fitness();
		fitness_function.printInfo(std::cout);
}


/** \brief Setup simulation system */
/** Read in parameters, construct and intialize required objects,
* and output simulation information to simulation directory 
*/
template<int dim>
void
Simulator<dim>::setup_system()
{
	std::cout << "\n\nSETTING UP SYSTEM\n" << std::endl;

// setup geometry and grid:
	setup_geometry_grid();

// setup velocity:
	setup_velocity();

// setup chemicals:
	setup_chemicals();

// setup bacteria:
	setup_bacteria();

	std::cout << std::endl << std::endl;
} // setup_system()

/** \brief Store local parameters read in from file */
template<int dim>
void
Simulator<dim>::assign_local_parameters()
{
	run_time = prm.get_double("Run time");
	save_period = prm.get_double("Save period");
}

/** \brief Setup suitable time steps for bacteria and chemicals */
template<int dim>
void
Simulator<dim>::setup_time_steps()
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
Simulator<dim>::get_chemical_time_step()
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
Simulator<dim>::setup_fitness()
{
	fitness_function.init(prm, chemicals);
}

/** \brief Declare local parameters needed for simulator 
* to be read in from file
*/
template<int dim>
void 
Simulator<dim>::declare_parameters()
{
	prm.declare_entry("Simulator type",
						"Bacteria",
						Patterns::Selection("Bacteria|IC_AS|Aging|Chemotaxis"),
						"Simulator type.");
	prm.declare_entry("Time step", "1", Patterns::Double());
	prm.declare_entry("Run time","0",Patterns::Double());
	prm.declare_entry("Save period","1",Patterns::Double());
	prm.declare_entry("Output directory","./",Patterns::Anything());

	// parameters for debugging: (only need to read in if in debug mode)
	prm.enter_subsection("Debug");
		prm.declare_entry("Chemicals only","False",Patterns::Bool());
		prm.declare_entry("Grid only","False",Patterns::Bool());
		prm.declare_entry("Record chemical mass","False",Patterns::Bool());
		prm.declare_entry("Track microbe chemicals","False",Patterns::Bool());
		prm.declare_entry("Convergence check","False",Patterns::Bool());


		prm.enter_subsection("Gaussian");
			prm.declare_entry("Centers",
								"{{}}",
								Patterns::List(Patterns::List(Patterns::Double())));
			prm.declare_entry("Amplitudes","{0,0}",
								Patterns::List(Patterns::Double()));
			prm.declare_entry("Widths","{0,0}",
								Patterns::List(Patterns::Double()));
		prm.leave_subsection();

	prm.leave_subsection();

	// declare class based parameters:
	Velocity::AdvectionHandler<dim>::declare_parameters(prm);
	GeometryTools::GeometryBuilder<dim>::declare_parameters(prm);
	Chemicals::ChemicalHandler<dim>::declare_parameters(prm);
	Chemicals::Controls<dim>::declare_parameters(prm);
	Bacteria::BacteriaHandler<dim>::declare_parameters(prm);
	Bacteria::Fitness::declare_parameters(prm);
}


// METHODS TO CHECK CONVERGENCE:
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------

// TODO: get exact and initial solutions from parameters (depends on bcs too)
// resetting system, redoing grid ...
// outputting final results

/** \brief CONVERGENCE CHECKS WITH MESH REFINEMENT */
/** Solutions are time dependent so we may want to querry accuracy at different 
* time steps.
*/
template<int dim>
void
Simulator<dim>::run_convergence_check()
{	
	TestFunctions::GaussianSolution<dim>::declare_parameters(prm);

	setup_parameters();
	setup_geometry_grid();
	setup_velocity();
	setup_chemicals();

	// convergence tables for each chemical:
	std::vector<dealii::ConvergenceTable> convergence_tables; // one for each chemical
	for(unsigned int i = 0; i < chemicals.getNumberChemicals(); ++i)
		convergence_tables.emplace_back( dealii::ConvergenceTable() );

	std::cout << std::endl << std::endl
		<< "Running Convergence checks..." << std::endl
		<< Utility::long_line << std::endl
		<< Utility::long_line << std::endl << std::endl;

	// solutions (for each chemical):
	std::vector<TestFunctions::GaussianSolution<dim> > gaussian_functions;
	for(unsigned int i = 0; i < chemicals.getNumberChemicals(); ++i)
		gaussian_functions.emplace_back(
			TestFunctions::GaussianSolution<dim>(prm, geometry, i) );

	// try one cycle first
	unsigned int cycle = 0; // for loop over total
	{
		std::cout << "CYCLE " << cycle << ":" << std::endl
			<< Utility::medium_line << std::endl;

		reset_system(cycle, gaussian_functions);	// ***refine mesh and project intial function 
		solve_chemicals(); // evolve chemicals upto given time ... need cycle for saving?



		for(unsigned int i = 0; i < chemicals.getNumberChemicals(); ++i)
		{
			gaussian_functions[i].setTime(run_time + 1.0); // start at time 1.0
			chemicals[i].process_solution(convergence_tables[i], gaussian_functions[i], cycle);
		}

		std::cout << std::endl << std::endl;
	} // for all cycles

	// OUTPUT
    std::cout << std::endl;
    for(unsigned int i = 0; i < chemicals.getNumberChemicals(); ++i)
    {
        convergence_tables[i].write_text(std::cout);

        std::string conv_filename = "convergence_chemical";

        conv_filename += "_" + dealii::Utilities::int_to_string(i, 3);
        conv_filename += ".tex";
        std::ofstream table_file(conv_filename);
        convergence_tables[i].write_tex(table_file);
    }
} // run_convergence_check()

template<int dim>
void 
Simulator<dim>::reset_system(const unsigned int cycle, 
	std::vector<TestFunctions::GaussianSolution<dim> >& gaussian_functions)
{
	// use cycle to refine mesh ...
	time = 0;
	time_step_number = 0;

	// recalc time steps??? or use min over refinements (could start with finest)

	// intialize chemicals:
	std::string section = "Debug.Gaussian";

	for(unsigned int i = 0; i < chemicals.getNumberChemicals(); ++i)
	{
		gaussian_functions[i].setTime(1.0);
		
		chemicals.project_function(gaussian_functions[i], i);
	} // for all chemicals

}

template<int dim>
void 
Simulator<dim>::solve_chemicals()
{
	std::cout << "SOLVING..." << std::endl
		<< Utility::short_line << std::endl;

	// save period:
	const unsigned int modsave
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );
	const bool isSavingChemicals = prm.get_bool("Chemicals","Save chemicals");
	const bool recordMass = prm.get_bool("Debug","Record chemical mass");

	do{
		// update time:
		time += chemical_time_step;
		++time_step_number;

		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;
			if(isSavingChemicals)
				output_chemicals();
			if(recordMass)
				output_chemical_mass();

			++save_step_number;
		}

		chemicals.update(); 

	}while( time < run_time );

}

template <int dim>
void 
Simulator<dim>::process_solution(const unsigned int cycle, 
	dealii::ConvergenceTable& convergence_table)
{

}



}} // CLOSE NAMESPACES
#endif
