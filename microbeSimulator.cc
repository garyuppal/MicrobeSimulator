#include "src/utility/command_line_parser.h"
// #include "src/simulators/full_simulator.h"

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
// - static assertions    Assert(false, ExcNotImplemented());
// - filter continuous mutation...

// - can already test random walk for cylinder geometry,
// just need to implement grid update and check flow

// - need to make sure bacteria time step does not cross channel ***or modify geometry class to better handle channel walls

// - add circle geometries for vortex and cylindrical pipe...
// - add a DG implementation??? // refactor chemicals	

// - can probably speed up reproduction/death as well as perhaps other things...
// - finish 3D implementations

// - eventually generalize bacteria types

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
			// sim.test_mesh();
		}
		else if(cmd_prm.getDimension() == 3)
		{
			std::cout << "Still need to implement some things for 3D" << std::endl;

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
