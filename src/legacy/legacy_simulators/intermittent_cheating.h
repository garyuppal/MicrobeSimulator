#ifndef MICROBESIMULATOR_INTERMITTENT_CHEATING_H
#define MICROBESIMULATOR_INTERMITTENT_CHEATING_H

#include <deal.II/grid/tria.h>
using dealii::Triangulation;

// #include "../bacteria/bacteria.h"
#include "../chemicals/chemical_fe_base.h"
#include "../chemicals/fe_chemical.h"
#include "../fitness/fitness_functions.h"
#include "../advection/stokes_solver.h"
#include "../geometry/geometry.h"
#include "../bacteria/bacteria_handler.h"

#include "../utility/argparser.h"
#include "../utility/cell_iterator_map.h"
#include "../utility/grid_generation_tools.h"
#include "../utility/simulation_tools.h"

#include <array>

namespace MicrobeSimulator{

template<int dim>
class IntermittentCheating{
public:
	IntermittentCheating(const ArgParser& para); // done

	void run();
private:
	static const unsigned int numchem = 2;

	const ArgParser* const 								parameters;

	Triangulation<dim>									triangulation;
	Geometry<dim> 										geometry; 
	FETools::PointCellMap<dim>							point_cell_map; 

	Velocity::AdvectionHandler<dim>						velocity_function;
	Chemicals::Chemical_FE_Base<dim>					chemical_fe_base; 
	std::array<Chemicals::FE_Chemical<dim>, numchem> 	chemicals; 

	MultiBacteria::BacteriaHandler<dim> 				bacteria;
	FitnessFunctions::TwoChemFitness<dim> 				fitness_function;


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
	void setup_system_constants(); // done 
	void setup_time_steps(); // done
	void setup_chemicals();  // done
	void setup_cell_map(); // done
	void setup_bacteria(); // done
	void setup_fitness(); // done
	double get_chemical_time_step(); // done 

	// run:
	void run_microbes(); // done
	void update_chemicals(); // done
	void update_bacteria();  // done
	void output_bacteria() const; // done
};


// IMPLEMENTATION:
// ----------------------------------------------------------------------------------------
template<int dim>
IntermittentCheating<dim>::IntermittentCheating(const ArgParser& para)
	:
	parameters(&para),
	triangulation(Triangulation<dim>::maximum_smoothing),
	chemical_fe_base(triangulation),
	output_directory("."),
	time(0),
	save_step_number(0),
	time_step_number(0)
{}

template<int dim>
void 
IntermittentCheating<dim>::setup_system()
{
	std::cout << "\n\nSETTING UP SYSTEM\n" << std::endl;

	setup_system_constants();
	SimulationTools::setup_geometry(*parameters, geometry);

	// *** update this with argparser ...
	std::array<BoundaryCondition, dim> bcs = {BoundaryCondition::REFLECT,
		BoundaryCondition::REFLECT}; 
	geometry.setBoundaryConditions(bcs);
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
	setup_bacteria(); // need to implement, local function
	setup_fitness();
}

template<int dim>
void 
IntermittentCheating<dim>::setup_system_constants()
{
	std::cout << "...Setting up system constants" << std::endl;
	run_time = parameters->getRunTime();
	save_period = parameters->getSavePeriod();

	output_directory = parameters->getOutputDirectory();
}

template<int dim>
void 
IntermittentCheating<dim>::setup_time_steps()
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
IntermittentCheating<dim>::get_chemical_time_step()
{
	const double min_time_step = 0.01; // will depend on diffusion parameters...
	double maximal_velocity = velocity_function.get_maximum_velocity(0); // SimulationTools::getMaxVelocity(geometry, velocity_function);
	maximal_velocity = std::max(maximal_velocity, 0.1);

	// const double maximal_velocity = std::max(velocity_function.get_maximal_velocity(), 0.1);
	double cfl_time_step = dealii::GridTools::minimal_cell_diameter(triangulation) 
							/ maximal_velocity;

	cfl_time_step = std::min(min_time_step, cfl_time_step);
	std::cout << "CFL_TIME_STEP: " << cfl_time_step << std::endl;

	return parameters->getTimeStepFactor()*cfl_time_step;
}


template<int dim>
void 
IntermittentCheating<dim>::setup_chemicals()
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
IntermittentCheating<dim>::setup_cell_map()
{
	std::cout << "...Setting up point-cell field" << std::endl;
	point_cell_map.initialize(geometry, 
							chemical_fe_base.get_dof_handler(), 
							parameters->getSourceResolution());

	point_cell_map.printInfo(std::cout); 
}

template<int dim>
void 
IntermittentCheating<dim>::setup_bacteria()
{
	std::cout << "... setting up bacteria" << std::endl;

	std::vector<Point<2> > initial_locations = parameters->getInitialLocations();

	if(initial_locations.empty())
	{
		unsigned int number_groups = parameters->getNumberGroups();
		if(number_groups == 0)
			number_groups = parameters->getNumberBacteria();
		const double left_length = parameters->getLeftSubdomainLength();
		initial_locations = SimulationTools::get_bacteria_locations(geometry, 
			number_groups, left_length);
	}

	bacteria.init(parameters->getBacteriaDiffusionConstant(),
		parameters->getNumberRegularBacteria(),
		parameters->getNumberICBacteria(),
		parameters->getNumberASBacteria(),
		// parameters->getICSwitchingRate(),
		parameters->getICOnPeriod(),
		parameters->getICOffPeriod(),
		parameters->getGoodSecretionRate(),
		parameters->getWasteSecretionRate(),
		initial_locations);

	bacteria.printInfo(std::cout);
}

template<int dim>
void 
IntermittentCheating<dim>::setup_fitness()
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

template<int dim>
void 
IntermittentCheating<dim>::update_chemicals()
{
	chemicals[0].update(bacteria.getAllLocations(), 
		bacteria.getAllGoodRates());
	
	chemicals[1].update(bacteria.getAllLocations(), 
		bacteria.getAllWasteRates());
}

template<int dim>
void 
IntermittentCheating<dim>::update_bacteria()
{
	// move:
	bacteria.randomWalk(bacteria_time_step, geometry, velocity_function); 
	bacteria.updateTimeAndSecretion(bacteria_time_step); // may include altruistic suicide later...

	// reproduce:
	bacteria.reproduce_simple(bacteria_time_step, fitness_function); 

	// mutate:
	const double species_mutation_rate = parameters->getSpeciesMutationRate();
	const double reg_mutation_rate = parameters->getMutationRate();
	const double original_secretion_rate = parameters->getGoodSecretionRate(); 

	bacteria.mutate_binary(bacteria_time_step, species_mutation_rate,
		reg_mutation_rate, original_secretion_rate); 
}

template<int dim>
void 
IntermittentCheating<dim>::output_bacteria() const
{
	std::string outfile = output_directory
							+ "/bacteria_" 
							+ dealii::Utilities::int_to_string(save_step_number,4)
							+ ".dat";
	std::ofstream out(outfile);
	bacteria.print(out);
}


// run:
// ------------------------------------------------------------------------------------
template<int dim>
void 
IntermittentCheating<dim>::run_microbes()
{
	std::cout << std::endl << std::endl
		<< "Starting Intermittent Cheating Microbe Simulation" << std::endl
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
	output_bacteria(); // need local function
	++save_step_number;

	do{
		// update time:
		time += chemical_time_step; 
		++time_step_number;

		// std::cout << "total bacteria: " << bacteria.getTotalNumber() << std::endl;
		
		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;

			if(parameters->isSavingChemicals())
				SimulationTools::output_chemicals<dim,2>(chemicals,
					output_directory, save_step_number);

			output_bacteria();
			++save_step_number;
		}

		update_chemicals();

		if(time_step_number % bacteria_time_step_multiplier == 0)
			update_bacteria();

	   	if( !bacteria.isAlive() )
	   		std::cout << "\n\nEverybody died!" << std::endl;	
	}while( (time < run_time) && bacteria.isAlive() );

}

template<int dim>
void
IntermittentCheating<dim>::run()
{
	std::cout <<"\n\n...running bacteria simulator" << std::endl;

	setup_system();
	run_microbes(); 
}

} // close namespaces
#endif