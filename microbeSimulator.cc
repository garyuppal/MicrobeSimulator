#include "src/utility/command_line_parser.h"
#include "src/simulators/full_simulator.h"

#include "src/simulators/simulator.h"

#include <list>
#include <fstream>
#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <chrono>
#include <functional>

// *************************************************************************************
// todo:
//               Assert(false, ExcNotImplemented());

// print info for velocity

// - add easy chemical debug ...

//*** parameter handler fails if no matching end in subsection for parameter files
// maybe invoke a precheck that gives a clearer error

// @todo add parameter handler pattern checking
// also should give error if file given doesn't exist!!	

// MASS SEEMS TO CHANGE AS WE CROSS DIFFERENT REFINEMENT AREAS
// not sure if this is due to wrong integration or actual bug (downward spikes)
// *** check also by getting chemical values at location of microbe over time
// also test without flow, but perhaps small flow for microbes only...


// - add circle geometries for vortex and cylindrical pipe...
// - add a DG implementation??? // refactor chemicals	


// build tests to be sure of the following:
// 1) figure out if sticking is an issue (chemical mass conservation)
// 2) groups streching, not splitting?  with pipe shear -- but seem ok in filter

// also try simulation continuous version to check pattern types...

// ** first establish proof of concept, try binary mutation with lower rates

// switch to new fitness class

// *** error using middle refined and periodic bc


// ...
// eventually might be nice to be able to define new variables in
// configuration file, can then have some more flexibility with looping and such
// - can probably speed up reproduction/death as well as perhaps other things...

// - later: generalize to 3d ... can modify geomety parameter patterns accordingly
	// - get_point_list of parameter_handler.h, generalize...

// - eventually generalize bacteria types....

// if we want to incorporate command line parameters into the same class
// as parameter handler, could implement a static declare parameter method
// for simulator that we call in main before running simulation
// *************************************************************************************

int main(int argc, char** argv)
{
	dealii::MultithreadInfo::set_thread_limit(1);
	try
    {
		using namespace std::chrono;
		using namespace MicrobeSimulator;

		// start time (for timing)
		auto start = high_resolution_clock::now();

		// parse command line (also creates output directory):
		const CommandLineParameters cmd_prm(argc,argv);
		cmd_prm.print(std::cout);
		cmd_prm.output();

		unsigned int seed = 0;
		if(cmd_prm.isSeed())
			seed = cmd_prm.getSeed(); // fix seed
		else
			seed =  time(0) + cmd_prm.getJobID(); // initialze random seed

		std::cout << "Using Random Seed: " << seed << std::endl;
		srand(seed);


		// run simulation:
		if(cmd_prm.getDimension() == 2)
		{
			Simulation::Simulator<2> sim(cmd_prm);
			sim.run();
			// sim.run_convergence_check();			
				// FullSimulator<2> fsim(cmd_prm);
				// fsim.run();
		}
		else if(cmd_prm.getDimension() == 3)
		{
			std::cout << "Still need to implement some things" << std::endl;
			// 	FullSimulator<3> sim(cmd_prm);
			// 	sim.run();
		}

		// end time, to get run time:
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);

		std::cout << "\n\n\nTotal program run time: "
		<< duration.count() << " seconds\n\n" << std::endl;

    }
	catch (std::exception &exc)
	{
		std::cerr << std::endl << std::endl
		    << "----------------------------------------------------"
		    << std::endl;
		std::cerr << "Exception on processing: " << std::endl << exc.what()
		    << std::endl << "Aborting!" << std::endl
		    << "----------------------------------------------------"
		    << std::endl << std::endl;

		return 1;
	}
	catch (...)
    {
		std::cerr << std::endl << std::endl
		    << "----------------------------------------------------"
		    << std::endl;
		std::cerr << "Unknown exception!" << std::endl << "Aborting!"
		    << std::endl
		    << "----------------------------------------------------"
		    << std::endl << std::endl;
		return 1;
    } // try

  return 0;

} // main()
