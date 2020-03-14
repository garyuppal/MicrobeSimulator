#ifndef MICROBESIMULATOR_AGING_SIMULATION_H
#define MICROBESIMULATOR_AGING_SIMULATION_H

#include <deal.II/grid/tria.h>
using dealii::Triangulation;

// #include "../bacteria/bacteria.h"
#include "../bacteria/aging_cells.h"
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

namespace MicrobeSimulator{ namespace Simulators{

// should be able to run simulations for 2d and 3d ...
template<int dim>
class AgingSimulator{
public:
	AgingSimulator(const ArgParser& para);

	void run();

private:
	static const unsigned int numchem = 1;

	const ArgParser* const 					parameters;

	Triangulation<dim> 						triangulation;
	Geometry<dim>							geometry;
	FETools::PointCellMap<dim>				point_cell_map;

	Velocity::AdvectionHandler<dim> 		velocity_function; // use constant flow or square_pipe

	Chemicals::Chemical_FE_Base<dim> 		chemical_fe_base;
	Chemicals::FE_Chemical<dim> 			cooperative_factors;

	MultiBacteria::AgingCells<dim> 			cells;
	FitnessFunctions::AgingDeathProb<dim>	pdeath_function;

	// SYSTEM CONSTANTS:
	std::string 							output_directory;
	double 									run_time;
	double 									time;
	double 									chemical_time_step;
	double 									cell_time_step;
	double 									save_period;
	unsigned int 							save_step_number;
	unsigned int 							time_step_number;
	unsigned int 							cell_time_step_multiplier;

	// METHODS:
	void setup_system();
	void setup_system_constants(); 
	void setup_time_steps();
	void setup_chemicals(); 
	void setup_cell_map(); 
	void setup_aging_cells();
	void setup_pdeath(); 
	double get_chemical_time_step(); 

	// run:
	void run_cells();
	void update_chemicals();
	void output_cells();
	void output_cooperative_factors();
};

// IMPLEMENTATION:
// -----------------------------------------------------------------
template<int dim>
AgingSimulator<dim>::AgingSimulator(const ArgParser& para)
	:
	parameters(&para),
	chemical_fe_base(triangulation),
	output_directory("."),
	time(0),
	save_step_number(0),
	time_step_number(0)
{}

// METHODS:
template<int dim>
void 
AgingSimulator<dim>::setup_system()
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
	// SimulationTools::setup_bacteria(*parameters, geometry, bacteria);
	setup_aging_cells();
	setup_pdeath();
}

template<int dim>
void 
AgingSimulator<dim>::setup_system_constants()
{
	std::cout << "...Setting up system constants" << std::endl;
	run_time = parameters->getRunTime();
	save_period = parameters->getSavePeriod();

	output_directory = parameters->getOutputDirectory();
}

template<int dim>
void 
AgingSimulator<dim>::setup_time_steps()
{
	std::cout << "...Setting up time steps" << std::endl;

	chemical_time_step = get_chemical_time_step(); // can move all to setup_time_steps();
	std::cout << "...using chemical time step " << chemical_time_step << std::endl;

	const double cell_max_time_step = 0.01;
	cell_time_step_multiplier = std::max( 
			static_cast<int>( std::floor(cell_max_time_step / chemical_time_step) ) , 1);
	cell_time_step = cell_time_step_multiplier*chemical_time_step;

	std::cout << "...using bacteria time step " << cell_time_step << std::endl;
}

template<int dim>
void 
AgingSimulator<dim>::setup_chemicals()
{
	std::cout << "...Setting up chemicals" << std::endl;

	// chemical base:
	chemical_fe_base.setup(velocity_function);
	chemical_fe_base.attach_point_cell_map(point_cell_map);

	// public good:
	cooperative_factors.reinit(chemical_fe_base,
						 parameters->getGoodDiffusionConstant(),
						 parameters->getGoodDecayConstant(),
						 parameters->getViscosityBeta(),
						 chemical_time_step);
	cooperative_factors.printInfo(std::cout);
}

template<int dim>
void 
AgingSimulator<dim>::setup_cell_map() // this should be moved to simulation tools...
{
	std::cout << "...Setting up point-cell field" << std::endl;
	point_cell_map.initialize(geometry, 
							chemical_fe_base.get_dof_handler(), 
							parameters->getSourceResolution());

	point_cell_map.printInfo(std::cout); 
}

template<int dim>
void 
AgingSimulator<dim>::setup_aging_cells()
{
	std::vector<Point<dim> > locations = parameters->getInitialLocations();

	if(locations.empty())
		cells.init(parameters->getGoodSecretionRate(),
					geometry.getBottomLeftPoint(),
					geometry.getTopRightPoint(),
					parameters->getCellDensity(),
					parameters->getCoarseGrainResolution());
	else
		cells.init(parameters->getGoodSecretionRate(),
					locations);

	cells.printInfo(std::cout);
}

template<int dim>
void 
AgingSimulator<dim>::setup_pdeath()
{
	std::cout << "...Setting up death rate function" << std::endl;

	pdeath_function.attach_chemicals(cooperative_factors);
	pdeath_function.set_fitness_constants(parameters->getAlphaGood(),
										parameters->getSecretionCost(), // using beta for k ***
										parameters->getGoodSaturation());
}

template<int dim>
double 
AgingSimulator<dim>::get_chemical_time_step() // could probably move this to simulation tools too...
{
	const double min_time_step = 0.01;
	double maximal_velocity = velocity_function.get_maximum_velocity(0); // SimulationTools::getMaxVelocity(geometry, velocity_function);
	maximal_velocity = std::max(maximal_velocity, 0.1);

	// const double maximal_velocity = std::max(velocity_function.get_maximal_velocity(), 0.1);
	double cfl_time_step = dealii::GridTools::minimal_cell_diameter(triangulation) 
							/ maximal_velocity;

	cfl_time_step = std::min(min_time_step,  
		parameters->getTimeStepFactor()*cfl_time_step);
	std::cout << "CFL_TIME_STEP: " << cfl_time_step << std::endl << std::endl;

	return cfl_time_step;
}

// UPDATE AND OUTPUT:
// -----------------------------------------------------------
template<int dim>
void 
AgingSimulator<dim>::update_chemicals()
{
	cooperative_factors.update(cells.getLocations(), cells.getSecretionRateVector()); // secretion rates is scalar...
}

template<int dim>
void				
AgingSimulator<dim>::output_cooperative_factors()
{
	cooperative_factors.output_solution(output_directory, 0, save_step_number);
}

template<int dim>
void 
AgingSimulator<dim>::output_cells()
{
	std::string outfile = output_directory
							+ "/cells_" 
							+ dealii::Utilities::int_to_string(save_step_number,4)
							+ ".dat";
	std::ofstream out(outfile);
	cells.print(out);
}


// RUN:
// ------------------------------------------------------------
template<int dim>
void 
AgingSimulator<dim>::run_cells()
{
	std::cout << std::endl << std::endl
		<< "Starting cell aging simulation" << std::endl
		<<"--------------------------------------------------------------------" 
		<< std::endl << std::endl;

	// save period:
	const unsigned int modsave 
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );

	// checking mass:
	std::vector<double> mass_check;

	do{
		// update time:
		time += chemical_time_step; 
		++time_step_number;

		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;

			if(parameters->isSavingChemicals())
				output_cooperative_factors();

			output_cells();
			++save_step_number;
		}

		mass_check.emplace_back(cooperative_factors.getMass());
		update_chemicals();

		if(time_step_number % cell_time_step_multiplier == 0)
			cells.random_death(cell_time_step, pdeath_function);

	   	if( !cells.isAlive() )
	   		std::cout << "\n\nEverybody died!" << std::endl;	

	} while( (time < run_time) && cells.isAlive() );

	SimulationTools::output_vector(mass_check,"mass_check", output_directory);
	SimulationTools::output_vector( cells.getSecretionRateVector(),
									"secretion rates", output_directory);
}

template<int dim>
void
AgingSimulator<dim>::run()
{
	std::cout <<"\n\n...running aging simulation" << std::endl;

	setup_system();
	run_cells();
}

}} // close namespace
#endif