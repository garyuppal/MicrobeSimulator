#ifndef MICROBESIMULATOR_BACTERIA_SIMULATOR_H
#define MICROBESIMULATOR_BACTERIA_SIMULATOR_H

#include <deal.II/grid/tria.h>
using dealii::Triangulation;

#include "../bacteria/bacteria.h"
#include "../chemicals/chemical_fe_base.h"
#include "../chemicals/fe_chemical.h"
#include "../fitness/fitness_functions.h"
// #include "../advection/stokes_solver.h"
// #include "../advection/velocity_functions.h"
#include "../advection/advection_handler.h"
#include "../geometry/geometry.h"

#include "../utility/argparser.h"
#include "../utility/cell_iterator_map.h"
#include "../utility/grid_generation_tools.h"
#include "../utility/simulation_tools.h"

#include <array>
#include <algorithm>

namespace MicrobeSimulator{ 

template<int dim>
class BacteriaSimulator{
public:
	BacteriaSimulator(const ArgParser& para); // done
	// ~BacteriaSimulator();

	void run_cycles();
	void run();

private:
	static const unsigned int numchem = 2;

	const ArgParser* const								parameters;

	Triangulation<dim>									triangulation;
	Geometry<dim> 										geometry; 
	FETools::PointCellMap<dim>							point_cell_map; 

	Velocity::AdvectionHandler<dim>						velocity_function;

	Chemicals::Chemical_FE_Base<dim>					chemical_fe_base; 
	std::array<Chemicals::FE_Chemical<dim>, numchem> 	chemicals; 
	Bacteria<dim, numchem> 								bacteria; 
	FitnessFunctions::TwoChemFitness<dim> 				fitness_function;

	std::vector<double> 								pg_rates; // for restarting bacteria

	// SYSTEM CONSTANTS:
	std::string 										output_directory;
	double 												run_time;
	double 												time;
	double 												chemical_time_step;
	double 												bacteria_time_step;
	double 												save_period; 
	unsigned int 										save_step_number;
	unsigned int 										time_step_number;
	unsigned int 										bacteria_time_step_multiplier;

	// for adding in new groups:
	unsigned int 										number_reintroduced;

	// METHODS:
	// setup:
	void setup_system(); // done
	void setup_system_constants(); // done
	void setup_time_steps(); // done
	void setup_chemicals(); // done
	void setup_cell_map(); // done
	void setup_fitness(); // done
	double get_chemical_time_step(); // .....

	// run:
	void run_microbes(unsigned int cycle); // done
	void run_microbes();
	void update_chemicals(); // done
	void update_bacteria(); //  -- add mutation...
	void reintro_bacteria(); 
	void reset_time();
	void reset_bacteria();
	void reset_chemicals();
};



// IMPLEMENTATION:
// ----------------------------------------------------------------------------------------
template<int dim>
BacteriaSimulator<dim>::BacteriaSimulator(const ArgParser& para)
	:
	parameters(&para),
	triangulation(Triangulation<dim>::maximum_smoothing),
	// velocity_function(triangulation),
	chemical_fe_base(triangulation),
	output_directory("."),
	time(0),
	save_step_number(0),
	time_step_number(0),
	number_reintroduced(0)
{}

	
template<int dim>
void 
BacteriaSimulator<dim>::setup_system()
{
	std::cout << "\n\nSETTING UP SYSTEM\n" << std::endl;

	setup_system_constants();
	SimulationTools::setup_geometry(*parameters, geometry);
	SimulationTools::output_geometry(output_directory, geometry);
	SimulationTools::setup_grid(*parameters, geometry, triangulation); 
	SimulationTools::output_grid(output_directory,"before_stokes_grid",triangulation);

	SimulationTools::setup_velocity(*parameters, 
									geometry, 
									velocity_function, 
									triangulation);

	SimulationTools::output_grid(output_directory,"after_stokes_grid",triangulation);

	setup_time_steps();
	setup_cell_map();
	// velocity_function.attach_point_cell_map(point_cell_map); // won't work ...
	
	setup_chemicals();
	SimulationTools::setup_bacteria(*parameters, geometry, bacteria);
	setup_fitness();
}


template<int dim>
void
BacteriaSimulator<dim>::setup_time_steps()
{
	std::cout << "...Setting up time steps" << std::endl;

	chemical_time_step = get_chemical_time_step(); // can move all to setup_time_steps();
	std::cout << "...using chemical time step " << chemical_time_step << std::endl;

	const double bacteria_max_time_step = 0.01;
	bacteria_time_step_multiplier = std::max( 
			static_cast<int>( std::floor(bacteria_max_time_step / chemical_time_step) ) , 1);
	bacteria_time_step = bacteria_time_step_multiplier*chemical_time_step;

	std::cout << "...using bacteria time step " << bacteria_time_step << std::endl;
}


template<int dim>
void 
BacteriaSimulator<dim>::setup_system_constants()
{
	std::cout << "...Setting up system constants" << std::endl;
	run_time = parameters->getRunTime();
	save_period = parameters->getSavePeriod();

	output_directory = parameters->getOutputDirectory();
}


template<int dim>
void 
BacteriaSimulator<dim>::setup_chemicals()
{
	std::cout << "...Setting up chemicals" << std::endl;

	// chemical base:
	chemical_fe_base.setup(velocity_function);
	chemical_fe_base.attach_point_cell_map(point_cell_map);

	// public good:
	chemicals[0].reinit(chemical_fe_base,
						 parameters->getGoodDiffusionConstant(),
						 parameters->getGoodDecayConstant(),
						 parameters->getViscosityBeta(),
						 chemical_time_step);
	chemicals[0].printInfo(std::cout);

	// waste chemical:
	chemicals[1].reinit(chemical_fe_base,
						 parameters->getWasteDiffusionConstant(),
						 parameters->getWasteDecayConstant(),
						 parameters->getViscosityBeta(),
						 chemical_time_step);
	chemicals[1].printInfo(std::cout);
}


template<int dim>
void 
BacteriaSimulator<dim>::setup_cell_map()
{
	std::cout << "...Setting up point-cell field" << std::endl;
	point_cell_map.initialize(geometry, 
							chemical_fe_base.get_dof_handler(), 
							parameters->getSourceResolution());

	point_cell_map.printInfo(std::cout); 
}


template<int dim>
void
BacteriaSimulator<dim>::setup_fitness()
{
	std::cout << "...Setting up fitness" << std::endl;

	fitness_function.attach_chemicals(chemicals[0], chemicals[1]);
	
	fitness_function.set_fitness_constants(parameters->getAlphaGood(),
											parameters->getAlphaWaste(),
											parameters->getGoodSaturation(),
											parameters->getWasteSaturation(),
											parameters->getSecretionCost());

	fitness_function.printInfo(std::cout);
}


/** NOTE: probably want to look at diffusion stabilty as well...
* ...also, can possibly use different time steps for different chemicals...
*/
template<int dim>
double 
BacteriaSimulator<dim>::get_chemical_time_step()
{
	const double min_time_step = 0.01;
	double maximal_velocity = velocity_function.get_maximum_velocity(0); // SimulationTools::getMaxVelocity(geometry, velocity_function);
	maximal_velocity = std::max(maximal_velocity, 0.1);

	// const double maximal_velocity = std::max(velocity_function.get_maximal_velocity(), 0.1);
	double cfl_time_step = dealii::GridTools::minimal_cell_diameter(triangulation) 
							/ maximal_velocity;
	std::cout << "CFL_TIME_STEP: " << cfl_time_step << std::endl;

	cfl_time_step = std::min( min_time_step, 
		(parameters->getTimeStepFactor())*cfl_time_step );
	std::cout << "using chem time step: " << cfl_time_step << std::endl << std::endl;

	return cfl_time_step;
}

// UPDATES:
// -------------------------------------------------------------------------------------------
template<int dim>
void 
BacteriaSimulator<dim>::update_chemicals()
{
	chemicals[0].update(bacteria.getLocations(), 
		bacteria.getSecretionRates(0));
	
	chemicals[1].update(bacteria.getLocations(), 
		bacteria.getSecretionRates(1));
}

template<int dim>
void 
BacteriaSimulator<dim>::reintro_bacteria()
{
	const unsigned int intro_number = 
		std::floor(time/(parameters->getReintroductionPeriod()));
	if(intro_number > number_reintroduced)
	{
		const unsigned int number_groups = parameters->getNumberGroups();
		const double left_subdomain_length = parameters->getLeftSubdomainLength();

		std::vector<Point<dim> > group_locations = 
		SimulationTools::get_bacteria_locations(geometry, 
			number_groups, left_subdomain_length); 

		bacteria.reintro(parameters->getNumberBacteria(),
						parameters->getGoodSecretionRate(),
						parameters->getWasteSecretionRate(),
						group_locations,
						geometry,
						bacteria_time_step); // for initial spread
		++number_reintroduced;
	}
}

template<int dim>
void 
BacteriaSimulator<dim>::update_bacteria()
{
	bacteria.randomWalk(bacteria_time_step, geometry, velocity_function, pg_rates); 

	// bacteria.reproduce_simple(bacteria_time_step, fitness_function);
	// bacteria.reproduce(bacteria_time_step, fitness_function);

	const double mutation_rate = parameters->getMutationRate();

	// change mutation to gaussian and get this from parameters separately
	// const double mutation_strength = 0.8*parameters->getGoodSecretionRate();

	const double original_secretion_rate = parameters->getGoodSecretionRate();
	// mutate chemical 0 
	// bacteria.mutate(0, bacteria_time_step, mutation_rate, mutation_strength);

	bacteria.mutate_binary(0,bacteria_time_step, mutation_rate, original_secretion_rate);
}

template<int dim>
void 
BacteriaSimulator<dim>::reset_time()
{
	time = 0;
	time_step_number = 0;
	save_step_number = 0;
}

template<int dim>
void 
BacteriaSimulator<dim>::reset_bacteria()
{	
	// shuffle secretion rates before reinit:
	std::random_shuffle(pg_rates.begin(), pg_rates.end());

	// keep only n of restarting bacteria:
	const unsigned int number_bacteria = parameters->getNumberBacteria();
	const unsigned int n_kept = pg_rates.size();
	const unsigned int n_bact_keep = std::min(n_kept, number_bacteria);	
	pg_rates.resize(n_bact_keep);

	std::cout << "... resetting " << n_bact_keep 
		<< " bacteria from " << n_kept 
		<< " collected from previous cycle"
		<< std::endl;

	// get group locations:
	const unsigned int number_groups = parameters->getNumberGroups();
	const double left_subdomain_length = parameters->getLeftSubdomainLength();

	std::vector<Point<dim> > group_locations = 
		SimulationTools::get_bacteria_locations(geometry, 
			number_groups, left_subdomain_length); 

	// reinitialize bacteria:
	bacteria.reinit(parameters->getBacteriaDiffusionConstant(), 
						pg_rates, // of length number_bacteria
						parameters->getWasteSecretionRate(),
						group_locations);

	// clear rates for next cycle:
	pg_rates.clear();
}


template<int dim>
void 
BacteriaSimulator<dim>::reset_chemicals()
{
	chemicals[0].project_initial_condition(dealii::ZeroFunction<dim>());
	chemicals[1].project_initial_condition(dealii::ZeroFunction<dim>());
}

// RUN:
// -------------------------------------------------------------------------------------------
template<int dim>	
void 
BacteriaSimulator<dim>::run_microbes(unsigned int cycle)
{
	std::cout << std::endl << std::endl
		<< "Starting microbe simulation" << std::endl
		<<"--------------------------------------------------------------------" << std::endl
		<< "\tRUNNING CYCLE: " << cycle << "...\n" << std::endl;

	// save period:
	const unsigned int modsave 
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );

	// spread out intial bacteria...
	const unsigned int intial_spread = 15;
	for(unsigned int i = 0; i < intial_spread; ++i)
		bacteria.randomWalk(chemical_time_step, geometry, velocity_function);

	// check after intial spread:
	SimulationTools::output_bacteria(bacteria,
				output_directory, save_step_number, cycle);
	++save_step_number;

	do{
		if(parameters->isAddingNewGroups())
			reintro_bacteria();

		// update time:
		time += chemical_time_step; 
		++time_step_number;

		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;

			if(parameters->isSavingChemicals())
				SimulationTools::output_chemicals<dim,2>(chemicals,
					output_directory, save_step_number, cycle);

			SimulationTools::output_bacteria(bacteria,
				output_directory, save_step_number, cycle);
			++save_step_number;
		}

		update_chemicals();

		if(time_step_number % bacteria_time_step_multiplier == 0)
			update_bacteria();
		
	   	if( !bacteria.isAlive() )
	   		std::cout << "\n\nEverybody died!" << std::endl;	
	}while( (time < run_time) && bacteria.isAlive() );
} // run_microbes()


template<int dim>
void
BacteriaSimulator<dim>::run_cycles()
{
	std::cout <<"\n\n...running bacteria simulator" << std::endl;

	setup_system();
	run_microbes(0); // initial cycle

	const unsigned int n_cycles = parameters->getNumberRunCycles();
	for(unsigned int cycle = 1; cycle < n_cycles; ++cycle)
	{
		reset_time();
		reset_bacteria(); 
		reset_chemicals();  
		run_microbes(cycle); 
	}
}

template<int dim>	
void 
BacteriaSimulator<dim>::run_microbes()
{
	std::cout << std::endl << std::endl
		<< "Starting microbe simulation" << std::endl
		<<"--------------------------------------------------------------------" 
		<< std::endl << std::endl;

	// save period:
	const unsigned int modsave 
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );

	// spread out intial bacteria...
	const unsigned int intial_spread = 15;
	for(unsigned int i = 0; i < intial_spread; ++i)
		bacteria.randomWalk(chemical_time_step, geometry, velocity_function);

	// check after intial spread:
	SimulationTools::output_bacteria(bacteria,
				output_directory, save_step_number);
	++save_step_number;

	do{
		if(parameters->isAddingNewGroups())
			reintro_bacteria();
		
		// update time:
		time += chemical_time_step; 
		++time_step_number;

		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;

			if(parameters->isSavingChemicals())
				SimulationTools::output_chemicals<dim,2>(chemicals,
					output_directory, save_step_number);

			SimulationTools::output_bacteria(bacteria,
				output_directory, save_step_number);
			++save_step_number;
		}

		update_chemicals();

		if(time_step_number % bacteria_time_step_multiplier == 0)
			update_bacteria();

	   	if( !bacteria.isAlive() )
	   		std::cout << "\n\nEverybody died!" << std::endl;	
	}while( (time < run_time) && bacteria.isAlive() );

	SimulationTools::output_vector(pg_rates,"pg_rates",output_directory);
} // run_microbes()


template<int dim>
void
BacteriaSimulator<dim>::run()
{
	std::cout <<"\n\n...running bacteria simulator" << std::endl;

	setup_system();
	run_microbes(); // initial cycle	
}

} // close namespace
#endif