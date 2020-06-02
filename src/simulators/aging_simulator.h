#pragma once

#include "../utility/parameter_handler.h"
#include "../utility/command_line_parser.h"

#include "../aging/chemicals.h"
#include "../aging/cells.h"
#include "../aging/fitness.h"

/** @file
* @todo extend model with extra cooperative factors and death of source cells
* @todo print out geometry for easier plotting/analysis
* @todo integrate in with existing code, possibly add FEM methods, 
*	may be faster with implicit solver
*/


namespace MicrobeSimulator{ 
	/** \brief Namespace for aging simulations */
	namespace Aging{

template<int dim>
class Simulator{
public:
	Simulator(const CommandLineParameters& cmd_prm);
	void run();

private:
	ParameterHandler prm;

	// Velocity::AdvectionHandler<dim> velocity_function; // need to modify for 1 dim
	// need geometry for 2 and 3 d

	Aging::Chemicals<dim>	chemicals;
	Aging::Cells<dim>		cells;
	Aging::Fitness<dim>		fitness;

	// SYSTEM CONSTANTS:
	std::string 							output_directory;
	double 									run_time;
	double 									time;
	double 									time_step;
	double 									save_period;
	unsigned int 							save_step_number;
	unsigned int 							time_step_number;
	double 									source_strength;

	// update:
	void update_chemicals();
	void update_cells();

	// output:
	void output_chemicals() const;
	void output_cells() const;

	// setup:
	void declare_parameters();
	void setup(); 
	void setup_parameters();
	void assign_local_parameters();
	void setup_chemicals();
	void setup_cells();
	void setup_fitness();

	void printInfo(std::ostream& out) const;
};

// IMPL
// -------------------------------------------------------------------------

template<int dim>
Simulator<dim>::Simulator(const CommandLineParameters& cmd_prm)
	:
	prm(cmd_prm.getParameterFile(), cmd_prm.getJobID()),
	output_directory(cmd_prm.getOutputDirectory()),
	run_time(0),
	time(0),
	time_step(1),
	save_step_number(0),
	time_step_number(0),
	source_strength(0)
{}

template<int dim>
void 
Simulator<dim>::run()
{
	setup();

	std::cout << std::endl << std::endl
		<< "Starting aging simulation" << std::endl
		<< Utility::long_line << std::endl
		<< Utility::long_line << std::endl << std::endl;

	// save period:
	const unsigned int modsave
		= static_cast<unsigned int>( std::ceil(save_period / time_step) );
	const bool isSavingChemicals = prm.get_bool("Chemicals","Save chemicals");

	// loop over time:
	do{
		// update chemicals:
		update_chemicals();
		update_cells(); 

		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;
			if(isSavingChemicals)
				output_chemicals();
			output_cells();
			++save_step_number;
		}

		// update time:
		time += time_step;
		++time_step_number;

		if(!cells.isAlive())
			std::cout << "Everyone died" << std::endl;

	}while( (time < run_time) && cells.isAlive() );
}

// ---------------------------------------------------------------------------------
// UPDATE
// ---------------------------------------------------------------------------------
template<int dim>
void 
Simulator<dim>::update_chemicals()
{
	// sinks given by cell locations,
	// soruces are constant and depend on mixing parameter
	std::vector<Point<dim> > source_locations = {Point<dim>()};
	std::vector<double> sources(source_locations.size(), source_strength);

	std::vector<Point<dim> > sink_locations = cells.getLocations();
	std::vector<double> sinks(sink_locations.size(), cells.getConsumptionRate());

	chemicals.update(source_locations, sources, sink_locations, sinks);
}

template<int dim>
void 
Simulator<dim>::update_cells()
{
	cells.update(time_step, fitness);
}

// ---------------------------------------------------------------------------------
// OUTPUT
// ---------------------------------------------------------------------------------
template<int dim>
void
Simulator<dim>::output_chemicals() const
{
	chemicals.output(output_directory, save_step_number);
}

template<int dim>
void
Simulator<dim>::output_cells() const
{
	std::string outfile = output_directory
						+ "/cells_"
						+ dealii::Utilities::int_to_string(save_step_number,4)
						+ ".dat";
	std::ofstream out(outfile);
	cells.print(out);	
}

// ---------------------------------------------------------------------------------
// SETUP
// ---------------------------------------------------------------------------------
template<int dim>
void 
Simulator<dim>::setup()
{
	std::cout << "\n\nSETTING UP AGING SIMULATION\n" << std::endl;

	setup_parameters();
	setup_chemicals();
	setup_cells();
	setup_fitness();

	std::cout << std::endl << std::endl;
}

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
	std::ofstream out_grid(output_directory + "/parameter_grid.dat");
	prm.printLoopedParameterGrid(std::cout);
	prm.printLoopedParameterGrid(out_grid);
	assign_local_parameters();
}

template<int dim>
void 
Simulator<dim>::assign_local_parameters()
{
	run_time = prm.get_double("Run time");
	time_step = prm.get_double("Time step");
	save_period = prm.get_double("Save period");
	source_strength = prm.get_double("Source strength");

	printInfo(std::cout);
}

template<int dim>
void 
Simulator<dim>::setup_chemicals()
{
	std::cout << "...setting up chemicals" << std::endl;
	chemicals.init(prm, time_step);
	chemicals.printInfo(std::cout);
}

template<int dim>
void 
Simulator<dim>::setup_cells()
{
	std::cout << "...setting up cells" << std::endl;
	cells.init(prm);
	cells.printInfo(std::cout);
}

template<int dim>
void
Simulator<dim>::setup_fitness()
{
	std::cout << "...setting up fitness" << std::endl;
	fitness.init(prm, chemicals);
	fitness.printInfo(std::cout);	
}

template<int dim>
void 
Simulator<dim>::declare_parameters()
{
	prm.declare_entry("Time step", "1", Patterns::Double());
	prm.declare_entry("Run time","0",Patterns::Double());
	prm.declare_entry("Save period","1",Patterns::Double());
	prm.declare_entry("Source strength","0",Patterns::Double());

	Aging::Chemicals<dim>::declare_parameters(prm);
	Aging::Cells<dim>::declare_parameters(prm);
	Aging::Fitness<dim>::declare_parameters(prm);
}

template<int dim>
void 
Simulator<dim>::printInfo(std::ostream& out) const
{
	out << "\n\n" << Utility::medium_line << std::endl 
		<< "\t\t SIMULATION INFO:" << std::endl
		<< Utility::medium_line << std::endl
		<< "Time step: " << time_step << std::endl
		<< "Run time: " << run_time << std::endl
		<< "Save period: " << save_period << std::endl
		<< "Source strength: " << source_strength << std::endl
		<< std::endl << Utility::medium_line << std::endl
		<< std::endl << std::endl;
}

}} // CLOSE NAMESPACES
/* aging_simulator.h */