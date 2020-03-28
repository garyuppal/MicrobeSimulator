#include "src/simulators/full_simulator.h"
#include "src/utility/command_line_parser.h"

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

// ** clean up geometry class

// increase # of -'s in print outs, can even call line a constant from utility
// short_line, long_line for example
// print out if no controls used ("none" type)
// print out velocity information

// try adding buffer for microbes, to eliminate sticking, done --- left buffer turned off
// also try simulation continuous version to check pattern types...

// ** first establish proof of concept, try binary mutation with lower rates

// ***bug in fitler trial run 8? also for run 37
// also add check to make sure time step is suitable for microbes with given diffusion constant
// and channel thickness, don't want them to go through...

// add ability to save chemicals as grid with fixed resolution 0.2 X done

// essential:
// 1. implement looping within list parameters // X should be good
// 2. have microbe reintroduction implemented if not done so already X
	// - save fallen and allow reintro even if extinct				X // done

// somewhere between 2 and 3, maybe less long channels (shoud be able to scale all with
// time and mutaiton rate). need many runs and look at output...

// eventually might be nice to be able to define variables in
// configuration file, can then have some more flexibility with looping and such

// *** can add a small slip velocity (can even make this a parameter)
// - run base tests:
	// initial group locations?
		// - check parsing, still bug with initial locations comment on same line
// - run filter and mixer tests (need initial location)
	// need reintroduction
	// *** for testing filters, initialize bacteria at point in long channel and compare
	// to pipe of same dimensions!!!
// - add chemical control (left off chemical_handler.h line 200)
// - problem with paramter handler needed jobID > 0 even if not looping? (set default to 1 for now)
// - can probably speed up reproduction/death as well as perhaps other things...
// - check geometry setups and clean up geometry class (ok for now 1/22/2020) and comments within file
// - add easy chemical debug ...
// - later: generalize to 3d ... can modify geomety parameter patterns accordingly
	// - get_point_list of parameter_handler.h, generalize...
// - cleaner print out of parameters ...
// - eventually generalize bacteria types....
// - todo's within "full_simulator.h" ...
// - add circle geometries for vortex and cylindrical pipe...
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

		// iniitialze random seed
		const unsigned int seed =  time(0) + cmd_prm.getJobID(); // + parameters.getJobID(); // 100; //
		std::cout << "Using Seed: " << seed << std::endl;
		srand(seed);

		// run simulation:
		if(cmd_prm.getDimension() == 2){
			FullSimulator<2> sim(cmd_prm);
			sim.run();
			// sim.test_chemical_flow();
		}
		// else if(cmd_prm.getDimension() == 3){
		// 	FullSimulator<3> sim(cmd_prm);
		// 	sim.run();

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
