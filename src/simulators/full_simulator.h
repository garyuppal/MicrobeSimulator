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

#include "./test_functions.h"

#include <array>
#include <algorithm>
#include <deal.II/base/function.h>

/** @file
* \brief Simulator file
* Simulation class defined here
* @todo add parameter handler pattern checking
* @todo add list option to parameter looping
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
	Bacteria::TestNewFitness::Fitness_Function<dim>	 test_ff;

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

	void run_chemicals_only(); // for debugging purposes
		void intialize_chemicals(); 

	void run_convergence_check();
	// // DEBUGGING METHODS:
	// void setup_debug();
	// void run_debug(); 
};

// IMPL
// -----------------------------------------------------------------

/** \brief Construct simulator */
template<int dim>
FullSimulator<dim>::FullSimulator(const CommandLineParameters& cmd_prm)
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
FullSimulator<dim>::run()
{
	setup_parameters();

	const bool gridOnly = prm.get_bool("Debug", "Grid only");
	const bool chemicalsOnly = prm.get_bool("Debug", "Chemicals only");
	const bool testConvergence = prm.get_bool("Debug", "Convergence check");

	// options here to run other simulations/tests:
	if(gridOnly)
	{
		setup_geometry_grid();
	}
	else if(chemicalsOnly)
	{
		setup_geometry_grid();
		setup_velocity();
		setup_chemicals();
		run_chemicals_only();
	}
	else if(testConvergence)
	{
		setup_geometry_grid();
		setup_velocity();
		setup_chemicals();
		run_convergence_check();
	}
	else
	{
		setup_system();
		run_microbes();
	}
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
		SimulationTools::output_vector(bacteria.get_pg_rates(),
			"pg_rates",output_directory);
} // run_microbes()

template<int dim>
void
FullSimulator<dim>::run_chemicals_only()
{
	std::cout << std::endl << std::endl
		<< "Starting chemicals debugging simulation" << std::endl
		<< Utility::long_line << std::endl
		<< Utility::long_line << std::endl << std::endl;

	// save period:
	const unsigned int modsave
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );
	const bool isSavingChemicals = prm.get_bool("Chemicals","Save chemicals");
	const bool recordMass = prm.get_bool("Debug","Record chemical mass");

	// intialize chem field:
	intialize_chemicals();

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

		// add option of only controls...
		chemicals.update(); 

	}while( time < run_time );

}

template<int dim>
void
FullSimulator<dim>::intialize_chemicals()
{
	std::string section = "Debug.Gaussian";
	const unsigned int one = 0;
	const unsigned int two = 1;

	Point<dim> center_one = prm.get_point_list(section, "Centers")[one];
	Point<dim> center_two = prm.get_point_list(section, "Centers")[two];

	double amp_one = prm.get_double_vector(section, "Amplitudes")[one];
	double amp_two = prm.get_double_vector(section, "Amplitudes")[two];

	double wid_one = prm.get_double_vector(section, "Widths")[one];
	double wid_two = prm.get_double_vector(section, "Widths")[two];

	TestFunctions::Gaussian<dim> gauss_one(center_one, amp_one, wid_one);
	TestFunctions::Gaussian<dim> gauss_two(center_two, amp_two, wid_two);

	chemicals.project_function(one, gauss_one);
	chemicals.project_function(two, gauss_two);
}

/** \brief CONVERGENCE CHECKS WITH MESH REFINEMENT */
/** Solutions are time dependent so we may want to querry accuracy at different 
* time steps.
*/
template<int dim>
void
FullSimulator<dim>::run_convergence_check()
{
	ConvergenceTable convergence_table;

	std::cout << std::endl << std::endl
		<< "Running Convergence checks..." << std::endl
		<< Utility::long_line << std::endl
		<< Utility::long_line << std::endl << std::endl;



}

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
FullSimulator<dim>::output_chemicals() const
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
FullSimulator<dim>::output_chemical_mass() const
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
FullSimulator<dim>::output_local_chemicals() const
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
FullSimulator<dim>::setup_parameters()
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
FullSimulator<dim>::setup_geometry_grid()
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
FullSimulator<dim>::setup_velocity()
{
	std::cout << "...setting up velocity" << std::endl;
		velocity_function.init(prm,geometry,triangulation,output_directory);
		GridGenerationTools::output_grid(output_directory,"after_stokes_grid",triangulation);
}

/** \brief Setup up chemicals and control functions */
template<int dim>
void
FullSimulator<dim>::setup_chemicals()
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
FullSimulator<dim>::setup_bacteria()
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
FullSimulator<dim>::setup_system()
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
template <int dim>
void 
FullSimulator<dim>::process_solution(const unsigned int cycle)
{
	Vector<float> difference_per_cell(triangulation.n_active_cells());
	VectorTools::integrate_difference(dof_handler,
	                              solution,
	                              Solution<dim>(),
	                              difference_per_cell,
	                              QGauss<dim>(fe->degree + 1),
	                              VectorTools::L2_norm);
	const double L2_error =
		VectorTools::compute_global_error(triangulation,
		                                difference_per_cell,
		                                VectorTools::L2_norm);

	VectorTools::integrate_difference(dof_handler,
	                              solution,
	                              Solution<dim>(),
	                              difference_per_cell,
	                              QGauss<dim>(fe->degree + 1),
	                              VectorTools::H1_seminorm);
	const double H1_error =
		VectorTools::compute_global_error(triangulation,
		                                difference_per_cell,
		                                VectorTools::H1_seminorm);

	const QTrapez<1>     q_trapez;
	const QIterated<dim> q_iterated(q_trapez, fe->degree * 2 + 1);
	VectorTools::integrate_difference(dof_handler,
	                              solution,
	                              Solution<dim>(),
	                              difference_per_cell,
	                              q_iterated,
	                              VectorTools::Linfty_norm);
	const double Linfty_error =
		VectorTools::compute_global_error(triangulation,
		                                difference_per_cell,
		                                VectorTools::Linfty_norm);

	const unsigned int n_active_cells = triangulation.n_active_cells();
	const unsigned int n_dofs         = dof_handler.n_dofs();

	std::cout << "Cycle " << cycle << ':' << std::endl
	      << "   Number of active cells:       " << n_active_cells
	      << std::endl
	      << "   Number of degrees of freedom: " << n_dofs << std::endl;
	convergence_table.add_value("cycle", cycle);
	convergence_table.add_value("cells", n_active_cells);
	convergence_table.add_value("dofs", n_dofs);
	convergence_table.add_value("L2", L2_error);
	convergence_table.add_value("H1", H1_error);
	convergence_table.add_value("Linfty", Linfty_error);
}



} // close namespace
#endif
