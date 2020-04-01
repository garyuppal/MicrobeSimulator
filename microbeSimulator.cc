#include "src/utility/command_line_parser.h"
#include "src/simulators/full_simulator.h"

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

// - add easy chemical debug ...
// change debug option in cmdprm to just a fixed seed option

// - add a DG implementation???

// build tests to be sure of the following:
// 1) figure out if sticking is an issue
// 2) groups streching, not splitting?  with pipe shear -- but seem ok in filter

// also try simulation continuous version to check pattern types...

// ** first establish proof of concept, try binary mutation with lower rates

// switch to new fitness class

// - add circle geometries for vortex and cylindrical pipe...


// ...
// eventually might be nice to be able to define new variables in
// configuration file, can then have some more flexibility with looping and such
// - can probably speed up reproduction/death as well as perhaps other things...

// - later: generalize to 3d ... can modify geomety parameter patterns accordingly
	// - get_point_list of parameter_handler.h, generalize...

// - eventually generalize bacteria types....
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
		CommandLineParameters cmd_prm(argc,argv);
		cmd_prm.print(std::cout);

		unsigned int seed;
		if(cmd_prm.isDebug())
		{
			// fix seed:	
			seed = 101;
			std::cout << "Using Fixed Seed: " << seed << std::endl;
		}
		else
		{
			// initialze random seed
			const unsigned int seed =  time(0) + cmd_prm.getJobID(); 
			std::cout << "Using Seed: " << seed << std::endl;
		}
		srand(seed);


		// run simulation:
		if(cmd_prm.getDimension() == 2)
		{
				FullSimulator<2> sim(cmd_prm);
				sim.run();
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
