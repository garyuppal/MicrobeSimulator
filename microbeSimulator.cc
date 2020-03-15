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
// @ todo:
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

		/// start time (for timing)
		auto start = high_resolution_clock::now();

		// parse command line:
		CommandLineParameters cmd_prm(argc,argv);
		// also creates output directory

		cmd_prm.print(std::cout);

		/// iniitialze random seed
		// const unsigned int seed =  time(0) + cmd_prm.getJobID(); // + parameters.getJobID(); // 100; //
		// std::cout << "Using Seed: " << seed << std::endl;
		// srand(seed);

		if(cmd_prm.getDimension() == 2){
			FullSimulator<2> sim(cmd_prm);
			sim.run();
			// sim.test_chemical_flow();
		}
		// else if(cmd_prm.getDimension() == 3){
		// 	FullSimulator<3> sim(cmd_prm);
		// 	sim.run();
/* ADD:
* get querry points() geometry issue with 3d?
*
*/

		/// end time, to get run time:
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
