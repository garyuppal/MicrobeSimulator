#pragma once

#include "../utility/parameter_handler.h"
#include "../utility/command_line_parser.h"

#include "../aging/chemicals.h"
#include "../aging/cells.h"
#include "../aging/fitness.h"

#include <utility> // pair
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
	unsigned int 							run_number;
	unsigned int 							number_runs;
	double 									reset_period;
	std::vector<double> 					reset_amps;
	

	// update:
	void update_chemicals();
	void update_cells();

	// output:
	void output_chemicals() const;
	void output_cells() const;
	void output_local_chemicals() const;
	void output_vector(const std::vector<double>& vect, const std::string& name) const;

	void output_cell_counts(const std::vector<std::pair<double,  unsigned int> >& ccnts, 
							unsigned int rn);

	// setup:
	void declare_parameters();
	void setup(); 
	void setup_parameters();
	void assign_local_parameters();
	void setup_chemicals();
	void setup_cells();
	void setup_fitness();

	// reset:
	void reset_system();
	void reset_chemicals();

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
	source_strength(0),
	run_number(0),
	number_runs(1),
	reset_period(0)
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
	const unsigned int mod_reset
		= static_cast<unsigned int>( std::ceil(reset_period / time_step) );
	const bool trackMicrobeChem = prm.get_bool("Track microbe chemicals");

	const bool isSavingCells = prm.get_bool("Cells", "Save cells");
	const bool isSavingCellCounts = prm.get_bool("Cells", "Save counts");

	std::vector<double> death_times;
	death_times.reserve(number_runs);

	do{
		std::cout << std::endl << std::endl
			<< Utility::short_line << std::endl
			<< " Starting run " << run_number+1 << " of " << number_runs << std::endl
			<< Utility::short_line << std::endl << std::endl;

		reset_system();

		std::vector< std::pair<double, unsigned int> > cell_counts;
		cell_counts.reserve( std::ceil( run_time / save_period ) );

		// loop over time:
		do{
			// update time:
			time += time_step;
			++time_step_number;

			// update chemicals:
			update_chemicals();
			update_cells(); 

			// reset chemicals as needed:
			if(time_step_number % mod_reset == 0)
				reset_chemicals();

			// output:
			if(time_step_number % modsave == 0)
			{
				std::cout << "saving at time: " << time << std::endl;
				if(isSavingChemicals)
					output_chemicals();
				if(trackMicrobeChem)
					output_local_chemicals();
				if( isSavingCells )
					output_cells();

				if( isSavingCellCounts )
					cell_counts.emplace_back( std::make_pair(time, cells.getSize() ) );
				++save_step_number;
			}



			if(!cells.isAlive())
			{
				std::cout << "Everyone died" << std::endl;
				death_times.emplace_back(time);
			}

		}while( (time < run_time) && cells.isAlive() );

		if( isSavingCellCounts )
			output_cell_counts(cell_counts, run_number);

		++run_number;
	}while(run_number < number_runs); // loop over runs

	output_vector(death_times, "death_times");

}

template<int dim>
void
Simulator<dim>::output_cell_counts(const std::vector<std::pair<double,  unsigned int> >& ccnts, 
							unsigned int rn)
{
	std::string cell_count_str = output_directory
								+ "/cell_counts_R" 
								+ dealii::Utilities::int_to_string(rn,4)
								+ ".dat";
	std::ofstream out(cell_count_str);

	const unsigned int n = ccnts.size();
	for(unsigned int i = 0; i < n; ++i)
		out << ccnts[i].first << " " << ccnts[i].second << std::endl;
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
	// std::vector<Point<dim> > source_locations = {Point<dim>()};
	// std::vector<double> sources(source_locations.size(), source_strength);

	// std::vector<Point<dim> > sink_locations = cells.getLocations();
	// std::vector<double> sinks(sink_locations.size(), cells.getConsumptionRate());

	chemicals.update(cells.getLocations(), cells.getSecretionRates(),
		 cells.getLocations(), cells.getConsumptionRates());
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
	chemicals.output(output_directory, save_step_number, run_number);
}

template<int dim>
void
Simulator<dim>::output_cells() const
{
	std::string outfile = output_directory
						+ "/cells_R"
						+ dealii::Utilities::int_to_string(run_number,4)
						+ "_"
						+ dealii::Utilities::int_to_string(save_step_number,4)
						+ ".dat";
	std::ofstream out(outfile);
	cells.print(out);	
}

template<int dim>
void
Simulator<dim>::output_local_chemicals() const
{
	const unsigned int nc = chemicals.getNumberChemicals();

	std::vector< Point<dim> > locations = cells.getLocations();
	for(unsigned int c = 0; c < nc; ++c)
	{
		std::vector<double> local_chems( chemicals[c].value_list(locations) );

		std::string outFile = output_directory 
				+ "/local_chem" 
				+ dealii::Utilities::int_to_string(c, 2)
				+ "_R"
				+ dealii::Utilities::int_to_string(run_number,4)
				+ "_"
				+ dealii::Utilities::int_to_string(save_step_number, 4)
				+ ".dat";
		std::ofstream out(outFile);

		for(unsigned int b = 0; b < local_chems.size(); ++b)
			out << locations[b] << " " << local_chems[b] << std::endl;
	}
}

template<int dim>
void 
Simulator<dim>::output_vector(
	const std::vector<double>& vect, const std::string& name) const
{
	std::string outFile = output_directory 
		+ "/"
		+ name
		+ ".dat";
	std::ofstream out(outFile);

	for(unsigned int i = 0; i < vect.size(); ++i)
		out << vect[i] << std::endl;
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
Simulator<dim>::reset_system()
{
	time = 0;
	save_step_number = 0;
	time_step_number = 0;
	reset_chemicals();
	// chemicals.init(prm, time_step); // project reset function instead***
	cells.init(prm);
}

template<int dim>
void
Simulator<dim>::reset_chemicals()
{
	const unsigned int n_chem = chemicals.getNumberChemicals();
	// const std::vector<double> amps = prm.get_double_vector("Reset amplitudes");

	for(unsigned int i = 0; i < n_chem; ++i)
		chemicals.project_function(
			ExactFunctions::Constant<dim>(reset_amps[i]) ,i);

	std::cout << "reset chems" << std::endl;
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
	number_runs = prm.get_unsigned("Number runs");
	reset_period = prm.get_double("Chemical reset period");
	reset_amps = prm.get_double_vector("Reset amplitudes");
	
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
	prm.declare_entry("Number runs","1",Patterns::Unsigned());
	
	prm.declare_entry("Chemical reset period","0",Patterns::Double());
	prm.declare_entry("Reset amplitudes","{0}",Patterns::List(Patterns::Double()));

	prm.declare_entry("Track microbe chemicals","False",Patterns::Bool());

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
		<< "Number runs: " << number_runs << std::endl
		<< std::endl << Utility::medium_line << std::endl
		<< std::endl << std::endl;
}

}} // CLOSE NAMESPACES
/* aging_simulator.h */