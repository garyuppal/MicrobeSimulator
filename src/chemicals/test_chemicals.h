#ifndef MICROBESIMULATOR_TEST_CHEMICALS_H
#define MICROBESIMULATOR_TEST_CHEMICALS_H

#include <deal.II/grid/tria.h>

#include "./fe_chemical.h"
#include "../advection/stokes_solver.h"
#include "../geometry/geometry.h"

#include "../utility/argparser.h"
#include "../utility/cell_iterator_map.h"
#include "../utility/grid_generation_tools.h"
#include "../utility/simulation_tools.h"
#include "../utility/exact_functions.h"

namespace MicrobeSimulator{ namespace Chemicals{
	using namespace dealii;

template<int dim>
class Test_Chemicals{
public:
	Test_Chemicals(const ArgParser& para);

	void run_test();
private:
	static const unsigned int numchem = 2;

	const ArgParser* const 					parameters;

	Triangulation<dim>						triangulation;

	Geometry<dim> 							geometry; 
	FETools::PointCellMap<dim>				point_cell_map; 
	Velocity::StokesSolver<dim>				stokes_solution;
	Chemical_FE_Base<dim>					chemical_fe_base; 
	std::array<FE_Chemical<dim>, numchem> 	chemicals; 

	// SYSTEM CONSTANTS:
	std::string 							output_directory;
	double 									run_time;
	double 									time;
	double 									chemical_time_step;
	double 									save_period; 
	unsigned int 							save_step_number;
	unsigned int 							time_step_number;

	// METHODS:
	// setup:
	void setup_system(); // done
	void setup_system_constants(); // done
	void setup_chemicals(); // done
	void setup_cell_map(); // done
	double get_chemical_time_step(); // done
};


// IMPLEMENTATION:
// -----------------------------------------------------------------------------
template<int dim>
Test_Chemicals<dim>::Test_Chemicals(const ArgParser& para)
	:
	parameters(&para),
	triangulation(Triangulation<dim>::maximum_smoothing),
	stokes_solution(triangulation),
	chemical_fe_base(triangulation),
	output_directory("."),
	time(0),
	save_step_number(0),
	time_step_number(0)
{}

template<int dim>
void 
Test_Chemicals<dim>::setup_system()
{
	std::cout << "\n\nSETTING UP CHEMICAL TEST SYSTEM\n" << std::endl;

	setup_system_constants();
	SimulationTools::setup_geometry(*parameters, geometry);
	SimulationTools::output_geometry(output_directory, geometry);
	SimulationTools::setup_grid(*parameters, geometry, triangulation); 
	SimulationTools::output_grid(output_directory,"before_stokes_grid",triangulation);
	SimulationTools::setup_stokes_velocity(parameters->getStokesRefinement(),
										parameters->getMaximumVelocity(),
										geometry.getNumberSpheres(),
										geometry.getNumberRectangles(),
										stokes_solution,
										output_directory);
	SimulationTools::output_grid(output_directory,"after_stokes_grid",triangulation);
	chemical_time_step = get_chemical_time_step();
	std::cout << "... using time step: " << chemical_time_step << std::endl;
	setup_cell_map();
	setup_chemicals();
}

template<int dim>
void 
Test_Chemicals<dim>::setup_system_constants()
{
	std::cout << "...Setting up system constants" << std::endl;
	run_time = parameters->getRunTime();
	save_period = parameters->getSavePeriod();

	output_directory = parameters->getOutputDirectory();
}

template<int dim>
void 
Test_Chemicals<dim>::setup_chemicals()
{
	std::cout << "...Setting up chemicals" << std::endl;

	// chemical base:
	chemical_fe_base.setup(stokes_solution);
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
Test_Chemicals<dim>::setup_cell_map()
{
	std::cout << "...Setting up point-cell field" << std::endl;
	point_cell_map.initialize(geometry, 
							chemical_fe_base.get_dof_handler(), 
							parameters->getSourceResolution());

	point_cell_map.printInfo(std::cout); 
}

template<int dim>
double 
Test_Chemicals<dim>::get_chemical_time_step()
{
	const double min_time_step = 0.01;

	const double maximal_velocity = std::max(stokes_solution.get_maximal_velocity(), 0.1);

	double cfl_time_step = dealii::GridTools::minimal_cell_diameter(triangulation) 
							/ maximal_velocity;

	// cfl_time_step = std::min(min_time_step, cfl_time_step);

	std::cout << "CFL_TIME_STEP: " << cfl_time_step << std::endl;

	return parameters->getTimeStepFactor()*cfl_time_step;
}

// RUN TEST:
// -------------------------------------------------------------------------------
template<int dim>
void 
Test_Chemicals<dim>::run_test()
{
	setup_system();

	// setup initial gaussians:
	// Point<dim> center = ( (dim == 2) ? Point<dim>(3, 2.5) : Point<dim>(2.5, 2.5, 2.5) );
	// const double height = 1.;
	// const double width = 0.5;

	// chemicals[0].project_initial_condition(ExactFunctions::Gaussian<dim>(center,height,width));
	// chemicals[1].project_initial_condition(ExactFunctions::Gaussian<dim>(center,height,width));

	// artificial sources:
	const std::vector<Point<dim> > locations = parameters->getInitialLocations(); // ({center});
	const std::vector<double> amounts(locations.size(), 1000);

	// save period:
	const unsigned int modsave 
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );

	const unsigned int n_save = static_cast<unsigned int>( std::ceil(run_time / chemical_time_step)) + 2;
	std::vector<double> saved_times, mass_c1, mass_c2, c1_min, c1_max, c2_min, c2_max;

	mass_c1.reserve(n_save);
	mass_c2.reserve(n_save);
	c1_min.reserve(n_save);
	c1_max.reserve(n_save);
	c2_min.reserve(n_save);
	c2_max.reserve(n_save);
	saved_times.reserve(n_save);

	// loop through time:
	do{
		// update chemicals:
		chemicals[0].update(locations, amounts);
		chemicals[1].update(locations, amounts);

		// update time:
		time += chemical_time_step; 
		++time_step_number;

		// mass:
		mass_c1.emplace_back(chemicals[0].getMass());
		mass_c2.emplace_back(chemicals[1].getMass());

		// c1 min max
		c1_min.emplace_back(chemicals[0].getMin());
		c1_max.emplace_back(chemicals[0].getMax());

		// c2 min max
		c2_min.emplace_back(chemicals[1].getMin());
		c2_max.emplace_back(chemicals[1].getMax());

		// record time (for easy plotting):
		saved_times.emplace_back(time);

		// save:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;
			if(parameters->isSavingChemicals())
				SimulationTools::output_chemicals<dim,2>(chemicals,
					output_directory,save_step_number);

			++save_step_number;
		}

	} while( time < run_time );

	// times:
	SimulationTools::output_vector(saved_times, "saved_times", output_directory);

	// masses
	SimulationTools::output_vector(mass_c1, "mass_c1", output_directory);
	SimulationTools::output_vector(mass_c2, "mass_c2", output_directory);

	// chemical 1 max/min
	SimulationTools::output_vector(c1_min, "c1_min", output_directory);
	SimulationTools::output_vector(c1_max, "c1_max", output_directory);

	// chemical 2 max/min
	SimulationTools::output_vector(c2_min, "c2_min", output_directory);
	SimulationTools::output_vector(c2_max, "c2_max", output_directory);
}


}} // close namespaces
#endif