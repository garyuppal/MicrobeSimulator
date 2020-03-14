#ifndef MICROBESIMULATOR_FULL_SIMULATOR_H
#define MICROBESIMULATOR_FULL_SIMULATOR_H

#include <deal.II/grid/tria.h>
using dealii::Triangulation;

// Simulator components:
// #include "../bacteria/bacteria.h"
#include "../bacteria/bacteria_handler2.h" // temp name
#include "../bacteria/bacteria_fitness.h"


#include "../chemicals/chemical_handler.h"
#include "../advection/advection_handler.h"
#include "../geometry/geometry.h"

// TOOLS:
// #include "../utility/argparser.h"
#include "./grid_generation_tools.h" // can include  a parameter declaration function here
#include "./simulation_tools.h"

#include "../utility/parameter_handler.h"
#include "../utility/command_line_parser.h"

#include <array>
#include <algorithm>
#include <deal.II/base/function.h>

/** @file
* \brief Simulator file
* Simulation class defined here
*/


namespace MicrobeSimulator{

/** \brief Simulator class
* This class handles reading in parameters, constructing, and executing the simulation.
*/
template<int dim>
class FullSimulator{
public:
	FullSimulator(const CommandLineParameters& cmd_prm);
	void run();
	void function_test();
	void test_chemical_flow();
private:
	ParameterHandler prm;

	Triangulation<dim>									triangulation;
	Geometry<dim> 										geometry;
	Velocity::AdvectionHandler<dim>						velocity_function;
	Chemicals::ChemicalHandler<dim> 					chemicals;

	Chemicals::Controls<dim>							control_functions;

	BacteriaNew::BacteriaHandler<dim>					bacteria;
	BacteriaNew::OR_Fitness<dim>						fitness_function;
		// type of fitness given by parameter as well***

	// std::vector<double>									pg_rates;

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

	double 												mutation_rate;
	double 												mutation_strength;

	// for adding in new groups:
	unsigned int 										number_reintroduced;

	// output:
	void output_bacteria() const;
	void output_chemicals() const;

	// update:
	void update_bacteria();

	// setup:
	void setup_system(); // can even move this outside of this class,
	void assign_local_parameters();
	void setup_time_steps();
	double get_chemical_time_step();
	void setup_fitness();
	// then can get dimension and simulation type in advance

	// Simulators:
	void run_microbes();

};

namespace TestFunction{
	template<int dim>
	class Gaussian : public dealii::Function<dim>{
	public:
		Gaussian(const Point<dim>& c, double a, double w);
		double value(const Point<dim>& p,
			const unsigned int component = 0) const override;
	private:
		Point<dim> center;
		double amplitude;
		double width;
	};

	template<int dim>
	Gaussian<dim>::Gaussian(const Point<dim>& c, double a, double w)
		:
		center(c), amplitude(a), width(w)
	{}

	template<int dim>
	double
	Gaussian<dim>::value(const Point<dim>& p,
			const unsigned int /* component */) const
	{
		return amplitude*std::exp( -(p-center)*(p-center)/width );
	}
} // close namespace

// IMPL
// -----------------------------------------------------------------
template<int dim>
FullSimulator<dim>::FullSimulator(const CommandLineParameters& cmd_prm)
	:
	prm(cmd_prm.getParameterFile(), cmd_prm.getJobID()),
	triangulation(Triangulation<dim>::maximum_smoothing),
	// chemical_fe_base(triangulation), // move into handler, can construct later in program
	// chemical_handler(4),
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

template<int dim>
void
FullSimulator<dim>::run()
{
	setup_system();
	run_microbes();
}

template<int dim>
void
FullSimulator<dim>::function_test()
{
	Point<dim> c = (dim==2) ? Point<dim>(2,2) : Point<dim>(2,2,2);
	double a = 1.0;
	double w = 1.0;

	TestFunction::Gaussian<dim> myfun(c,a,w);

	chemicals.project_function(0, myfun);
	std::vector<Point<dim> > qpoints = geometry.getQuerryPoints();

	for(unsigned int t = 0; t < 100; ++t)
	{
		// update field:
		chemicals.update();

		// get values:
		std::vector<double> values = chemicals[0].value_list(qpoints);

		//output:
		std::ofstream out(output_directory + "/chem_test_" + std::to_string(t) + ".txt");

		for(unsigned int i = 0; i < qpoints.size(); ++i)
			out << qpoints[i] << " " << values[i] << std::endl;

		std::cout << "done t = " << t << std::endl;
	}
}

template<int dim>
void
FullSimulator<dim>::run_microbes()
{
	std::cout << std::endl << std::endl
		<< "Starting microbe simulation" << std::endl
		<<"--------------------------------------------------------------------"
		<< std::endl << std::endl;

	// save period:
	const unsigned int modsave
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );
	const bool isSavingChemicals = prm.get_bool("Chemicals","Save chemicals");

	// spread out intial bacteria...
	const unsigned int intial_spread = 15;
	for(unsigned int i = 0; i < intial_spread; ++i)
		bacteria.randomWalk(chemical_time_step, geometry, velocity_function);

	// check after intial spread:
	output_bacteria();
	++save_step_number;

	do{
		// if(parameters->isAddingNewGroups())
		// 	reintro_bacteria();

		// update time:
		time += chemical_time_step;
		++time_step_number;

		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;
			if(isSavingChemicals)
				output_chemicals(); // *** try saving projected version to triple check***
			output_bacteria();
			++save_step_number;
		}

		// chemicals.update(bacteria.getAllLocations(), bacteria.getAllRates());

		control_functions.update_time(chemical_time_step); // increment internal clock
		chemicals.update(bacteria.getAllLocations(), bacteria.getAllRates(), control_functions);

		if(time_step_number % bacteria_time_step_multiplier == 0)
			update_bacteria();
			// std::cout << "number bacteria: " << bacteria.getTotalNumber() << std::endl;

	   	if( !bacteria.isAlive() )
	   		std::cout << "\n\nEverybody died!" << std::endl;
	}while( (time < run_time) && bacteria.isAlive() );

	// SimulationTools::output_vector(pg_rates,"pg_rates",output_directory);
}

template<int dim>
void
FullSimulator<dim>::test_chemical_flow()
{
	setup_system();

	std::cout << std::endl << std::endl
		<< "Starting chemical flow test" << std::endl
		<<"--------------------------------------------------------------------"
		<< std::endl << std::endl;

	// save period:
	const unsigned int modsave
		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );
	// const bool isSavingChemicals = prm.get_bool("Chemicals","Save chemicals");

	// spread out intial bacteria...
	// const unsigned int intial_spread = 15;
	// for(unsigned int i = 0; i < intial_spread; ++i)
	// 	bacteria.randomWalk(chemical_time_step, geometry, velocity_function);

	// check after intial spread:
	// output_bacteria();
	// ++save_step_number;

	Point<dim> c = (dim==2) ? Point<dim>(2,2) : Point<dim>(2,2,2);
	double a = 1.0;
	double w = 1.0;

	TestFunction::Gaussian<dim> myfun(c,a,w);

	chemicals.project_function(0, myfun);
	std::vector<Point<dim> > qpoints = geometry.getQuerryPoints();

	do{
		// if(parameters->isAddingNewGroups())
		// 	reintro_bacteria();

		// update time:
		time += chemical_time_step;
		++time_step_number;

		// output:
		if(time_step_number % modsave == 0)
		{
			std::cout << "saving at time: " << time << std::endl;
			// if(isSavingChemicals)
			output_chemicals(); // *** try saving projected version to triple check***
			// output_bacteria();
			++save_step_number;
		}

		// chemicals.update(bacteria.getAllLocations(), bacteria.getAllRates());
		chemicals.update();

		// if(time_step_number % bacteria_time_step_multiplier == 0)
		// {
		// 	update_bacteria();
		// 	std::cout << "number bacteria: " << bacteria.getTotalNumber() << std::endl;
		// }

	   	// if( !bacteria.isAlive() )
	   		// std::cout << "\n\nEverybody died!" << std::endl;
	}while( (time < run_time) );

}


// UPDATE:
// --------------------------------------------------------------------------
template<int dim>
void
FullSimulator<dim>::update_bacteria()
{
	bacteria.randomWalk(bacteria_time_step, geometry, velocity_function);

	bacteria.reproduce(bacteria_time_step, fitness_function);

	bacteria.mutate_simple(bacteria_time_step, mutation_rate, mutation_strength);
}

// OUTPUT:
// --------------------------------------------------------------------------

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

template<int dim>
void
FullSimulator<dim>::output_chemicals() const
{
	chemicals.output(output_directory, save_step_number);
}


// SETUP:
// ------------------------------------------------------------------------

template<int dim>
void
FullSimulator<dim>::setup_system()
{
	std::cout << "\n\nSETTING UP SYSTEM\n" << std::endl;

	std::cout << "...setting up parameters" << std::endl;
	SimulationTools::declare_common_parameters(prm); // @ todo: maybe make a local function
	Velocity::AdvectionHandler<dim>::declare_parameters(prm);
	Geometry<dim>::declare_parameters(prm);
	GridGenerationTools::declare_parameters(prm);
	Chemicals::ChemicalHandler<dim>::declare_parameters(prm);
	Chemicals::Controls<dim>::declare_parameters(prm);
	BacteriaNew::BacteriaHandler<dim>::declare_parameters(prm);
	BacteriaNew::Fitness::declare_parameters(prm);

	prm.parse_parameter_file();
	prm.print(std::cout);
	std::ofstream out(output_directory + "/parameters.dat");
	prm.print(out);

	prm.printLoopedParameterGrid(std::cout);

	assign_local_parameters();

	std::cout << "...setting up geometry" << std::endl;
		geometry.init(prm);
		geometry.printInfo(std::cout);
		SimulationTools::output_geometry(output_directory, geometry);

	std::cout << "...setting up grid" << std::endl;
		SimulationTools::setup_grid(prm, geometry, triangulation);
		SimulationTools::output_grid(output_directory,"before_stokes_grid",triangulation);

	std::cout << "...setting up velocity" << std::endl;
		velocity_function.init(prm,geometry,triangulation,output_directory);

	SimulationTools::output_grid(output_directory,"after_stokes_grid",triangulation);

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

template<int dim>
void
FullSimulator<dim>::assign_local_parameters()
{
	run_time = prm.get_double("Run time");
	save_period = prm.get_double("Save period");
	mutation_rate = prm.get_double("Bacteria", "Mutation rate");
	mutation_strength = prm.get_double("Bacteria", "Mutation strength");
}

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

template<int dim>
double
FullSimulator<dim>::get_chemical_time_step()
{
	const double min_time_step = 0.001;
	double maximal_velocity = velocity_function.get_maximum_velocity(0); // SimulationTools::getMaxVelocity(geometry, velocity_function);
	maximal_velocity = std::max(maximal_velocity, 0.1);

	// const double maximal_velocity = std::max(velocity_function.get_maximal_velocity(), 0.1);
	double cfl_time_step = dealii::GridTools::minimal_cell_diameter(triangulation)
							/ maximal_velocity;
	std::cout << "\tCFL_TIME_STEP: " << cfl_time_step << std::endl;

	cfl_time_step = std::min( min_time_step,
		(prm.get_double("Chemicals", "Time step factor"))*cfl_time_step );

	return cfl_time_step;
}

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



} // close namespace
#endif
