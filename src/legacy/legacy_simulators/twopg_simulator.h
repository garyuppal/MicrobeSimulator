#ifndef MICROBESIMULATOR_TWOPG_SIMULATOR_H
#define MICROBESIMULATOR_TWOPG_SIMULATOR_H

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
class TwoPGSimulator{
public:
	TwoPGSimulator(const ArgParser& para);

	void run();
private:
	static const unsigned int numchem = 3;

	const ArgParser* const								parameters;

	Triangulation<dim>									triangulation;
	Geometry<dim> 										geometry; 
	FETools::PointCellMap<dim>							point_cell_map; 

	Velocity::AdvectionHandler<dim>						velocity_function;

	Chemicals::Chemical_FE_Base<dim>					chemical_fe_base; 
	std::array<Chemicals::FE_Chemical<dim>, numchem> 	chemicals; 
	Bacteria<dim, numchem> 								bacteria; 
	// FitnessFunctions::TwoChemFitness<dim> 				fitness_function;
	FitnessFunctions::AND_Fitness<dim, numchem>			fitness_function;
	FitnessFunctions::OR_Fitness<dim, numchem>			or_fitness_function;

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

	// METHODS:
	// setup:
	void setup_system(); // done
	void setup_system_constants(); // done **
	void setup_time_steps(); // done **
	void setup_chemicals(); // done ***
	void setup_bacteria(); // done***
	void setup_cell_map(); // done **
	void setup_fitness(); // done -- could have virtual set functions, or better yet, constuctor,
		// so we can use runtime polymorphism... (set for AND currently)
	double get_chemical_time_step(); // done **

	// ** == copied for bacteria_simulator.h (can we put these elsewhere?)
	// and/or create a more general simulator class?

	// *** -> could probably be easily generalized

	// run:
	void run_microbes(); // done
	void update_chemicals(); // 
	void update_bacteria(); //  
};

// IMPLEMENTATION:
// ---------------------------------------------------------------------------------------
template<int dim>
TwoPGSimulator<dim>::TwoPGSimulator(const ArgParser& para)
	:
	parameters(&para),
	triangulation(Triangulation<dim>::maximum_smoothing),
	chemical_fe_base(triangulation),
	output_directory("."),
	time(0),
	save_step_number(0),
	time_step_number(0)
{
	dealii::MultithreadInfo::set_thread_limit(1);
}


template<int dim>
void 
TwoPGSimulator<dim>::setup_system()
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
	
	setup_chemicals();
	setup_bacteria();
	// SimulationTools::setup_bacteria(*parameters, geometry, bacteria);
	setup_fitness();
}

template<int dim>
void
TwoPGSimulator<dim>::setup_time_steps()
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
double 
TwoPGSimulator<dim>::get_chemical_time_step()
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

template<int dim>
void 
TwoPGSimulator<dim>::setup_system_constants()
{
	std::cout << "...Setting up system constants" << std::endl;
	run_time = parameters->getRunTime();
	save_period = parameters->getSavePeriod();

	output_directory = parameters->getOutputDirectory();
}

template<int dim>
void 
TwoPGSimulator<dim>::setup_cell_map()
{
	std::cout << "...Setting up point-cell field" << std::endl;
	point_cell_map.initialize(geometry, 
							chemical_fe_base.get_dof_handler(), 
							parameters->getSourceResolution());

	point_cell_map.printInfo(std::cout); 
}

template<int dim>
void 
TwoPGSimulator<dim>::setup_chemicals()
{
	std::cout << "...Setting up chemicals" << std::endl;

	// chemical base:
	chemical_fe_base.setup(velocity_function, geometry.getBoundaryConditions());
	chemical_fe_base.attach_point_cell_map(point_cell_map);

	// public good 1:
	chemicals[0].reinit(chemical_fe_base,
						 parameters->getGoodDiffusionConstant(),
						 parameters->getGoodDecayConstant(),
						 parameters->getViscosityBeta(),
						 chemical_time_step);
	chemicals[0].printInfo(std::cout);

	// public good 2:
	chemicals[1].reinit(chemical_fe_base,
						 parameters->getGoodDiffusionConstant(),
						 parameters->getGoodDecayConstant(),
						 parameters->getViscosityBeta(),
						 chemical_time_step);
	chemicals[1].printInfo(std::cout);

	// waste chemical:
	chemicals[2].reinit(chemical_fe_base,
						 parameters->getWasteDiffusionConstant(),
						 parameters->getWasteDecayConstant(),
						 parameters->getViscosityBeta(),
						 chemical_time_step);
	chemicals[2].printInfo(std::cout);
}

template<int dim>
void 
TwoPGSimulator<dim>::setup_bacteria()
{
	std::cout << "... setting up bacteria" << std::endl;

	std::vector<Point<2> > initial_locations = parameters->getInitialLocations();

	if(initial_locations.empty())
	{
		unsigned int number_groups = parameters->getNumberGroups();
		if(number_groups == 0)
			number_groups = parameters->getNumberBacteria();
		const double left_length = parameters->getLeftSubdomainLength();
		initial_locations = SimulationTools::get_bacteria_locations(
										geometry, number_groups, left_length);
	}

	bacteria.init(parameters->getBacteriaDiffusionConstant(),
		parameters->getNumberBacteria(),
		std::array<double, 3>(
		{parameters->getGoodSecretionRate(),
			parameters->getGoodSecretionRate(),
		 	parameters->getWasteSecretionRate()}),
		initial_locations); 

	bacteria.printInfo(std::cout);
}

template<int dim>
void
TwoPGSimulator<dim>::setup_fitness()
{
	std::cout << "...Setting up fitness" << std::endl;

	if(parameters->isUsingANDFitness())
	{
		fitness_function.attach_chemicals(
			std::array<Chemicals::ChemicalInterface<dim>*, numchem>(
				{&chemicals[0], &chemicals[1], &chemicals[2]})
			);
		
		fitness_function.set_fitness_constants(parameters->getAlphaGood(),
												parameters->getAlphaWaste(),
												parameters->getGoodSaturation(),
												parameters->getWasteSaturation(),
												parameters->getSecretionCost());

		fitness_function.printInfo(std::cout);
	}
	else
	{
		or_fitness_function.attach_chemicals(
			std::array<Chemicals::ChemicalInterface<dim>*, numchem>(
				{&chemicals[0], &chemicals[1], &chemicals[2]})
			);
		
		or_fitness_function.set_fitness_constants(parameters->getAlphaGood(),
												parameters->getAlphaWaste(),
												parameters->getGoodSaturation(),
												parameters->getWasteSaturation(),
												parameters->getSecretionCost());

		or_fitness_function.printInfo(std::cout);
	}
}


// UPDATE:
// --------------------------------------------------------------------------------
template<int dim>
void 
TwoPGSimulator<dim>::update_chemicals()
{
	for(unsigned int c = 0; c < numchem; ++c)
		chemicals[c].update(bacteria.getLocations(), 
			bacteria.getSecretionRates(c));
}

template<int dim>
void 
TwoPGSimulator<dim>::update_bacteria()
{
	// MOVE:
	bacteria.randomWalk(bacteria_time_step, geometry, velocity_function, pg_rates); 

	// REPRODUCE:
	if(parameters->isUsingANDFitness())
		bacteria.reproduce_simple2(bacteria_time_step, fitness_function); 
	else
		bacteria.reproduce_simple2(bacteria_time_step, or_fitness_function); 
		// *** need to check reprduce function and use as the more general case
	// *** possibly restructure bacteria back to old method .. can we use type punning to get 
		// access without copying??
	// bacteria.reproduce(bacteria_time_step, fitness_function);

	// MUTATE:
	const double mutation_rate = parameters->getMutationRate();
	const double original_secretion_rate = parameters->getGoodSecretionRate();
	
	// switch between two with probability (1/2)
	// *** probably not ``correct'' way to do so, but should be same as done in MATLAB code
	double prob = ((double) rand() / (RAND_MAX));
	if(prob < 0.5)
		bacteria.mutate_binary(0,bacteria_time_step, mutation_rate, original_secretion_rate);
	else
		bacteria.mutate_binary(1,bacteria_time_step, mutation_rate, original_secretion_rate);
}


// RUN:
// ---------------------------------------------------------------------------------

template<int dim>	
void 
TwoPGSimulator<dim>::run_microbes()
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
		// update time:
		time += chemical_time_step; 
		++time_step_number;

		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;

			if(parameters->isSavingChemicals())
				SimulationTools::output_chemicals<dim, 3>(chemicals,
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

} // run_microbes()



template<int dim>
void 
TwoPGSimulator<dim>::run()
{
	std::cout <<"\n\n...running TWOPG bacteria simulator" << std::endl;

	setup_system();
	run_microbes(); 
}


} // close namespace
#endif