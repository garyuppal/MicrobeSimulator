// #ifndef MICROBESIMULATOR_TESTS_H
// #define MICROBESIMULATOR_TESTS_H

// #include <deal.II/grid/tria.h>
// using dealii::Triangulation;

// #include "../bacteria/bacteria_handler.h" 
// #include "../bacteria/bacteria_fitness.h"
// #include "../chemicals/chemical_handler.h"
// #include "../advection/advection_handler.h"
// #include "../geometry/geometry.h"
// #include "../geometry/geometry_builder.h" 
// #include "../utility/parameter_handler.h"
// #include "../utility/command_line_parser.h"

// #include <array>
// #include <algorithm>
// #include <deal.II/base/function.h>

// // ideally we dont want to repeat too much code from simulator,
// // the point of this is to easily debug integrated problems
// // would be nicer to have individual unit tests for different components too however
// // and then test each of those here, maybe eventually...


// namespace MicrobeSimulator{ 
// 	/** \brief Namespace to implement some tests */
// 	namespace Tests{
// template<int dim>
// class Debugger{
// public:
// 	Debugger(const CommandLineParameters& cmd_prm);

// 	void run() {std::cout << "running debugger..." << std::endl;}
// private:
// 	ParameterHandler prm;

// };

// template<int dim>
// Debugger<dim>::Debugger(const CommandLineParameters& cmd_prm)
// 	:
// 	prm(cmd_prm.getParameterFile(), cmd_prm.getJobID())
// {}

// }} // CLOSE NAMESPACES

// /* TESTS:
// * conservation of mass when we have a steady source and decay rate
// * conservation of mass when there is only diffusion and advection, 
// 	- no decay or source, just an initial condition
// * motion of microbes with reproduction vs without
// * convergence of solutions as we refine the grid
// 	- known solutions to pdes?

// * fix seed and initial location of group

// * test framework:
// 	- may go top down at this point, can build unit tests later on (probably won't)
// 	- establish integrated tests to run, think about unit tests later
// 	- define set of debugging parameters... limit need to recompile
// */

// 	/** \brief Test function to debug chemicals */
// 	/** @todo Move this into a separate class to test simulator */
// // 	namespace TestFunction{
// // 		template<int dim>
// // 		class Gaussian : public dealii::Function<dim>{
// // 		public:
// // 			Gaussian(const Point<dim>& c, double a, double w);
// // 			double value(const Point<dim>& p,
// // 				const unsigned int component = 0) const override;
// // 		private:
// // 			Point<dim> center;
// // 			double amplitude;
// // 			double width;
// // 		};

// // 		template<int dim>
// // 		Gaussian<dim>::Gaussian(const Point<dim>& c, double a, double w)
// // 			:
// // 			center(c), amplitude(a), width(w)
// // 		{}

// // 		template<int dim>
// // 		double
// // 		Gaussian<dim>::value(const Point<dim>& p,
// // 				const unsigned int /* component */) const
// // 		{
// // 			return amplitude*std::exp( -(p-center)*(p-center)/width );
// // 		}
// // 	} // close namespace



// // /** \brief Test field with gaussian function */
// // /** @todo Move this into a separate class to test simulator */
// // template<int dim>
// // void
// // FullSimulator<dim>::function_test()
// // {
// // 	Point<dim> c = (dim==2) ? Point<dim>(2,2) : Point<dim>(2,2,2);
// // 	double a = 1.0;
// // 	double w = 1.0;

// // 	TestFunction::Gaussian<dim> myfun(c,a,w);

// // 	chemicals.project_function(0, myfun);
// // 	std::vector<Point<dim> > qpoints = geometry.getQuerryPoints();

// // 	for(unsigned int t = 0; t < 100; ++t)
// // 	{
// // 		// update field:
// // 		chemicals.update();

// // 		// get values:
// // 		std::vector<double> values = chemicals[0].value_list(qpoints);

// // 		//output:
// // 		std::ofstream out(output_directory + "/chem_test_" + std::to_string(t) + ".txt");

// // 		for(unsigned int i = 0; i < qpoints.size(); ++i)
// // 			out << qpoints[i] << " " << values[i] << std::endl;

// // 		std::cout << "done t = " << t << std::endl;
// // 	}
// // }


// // /** \brief Test chemical evolution */
// // template<int dim>
// // void
// // FullSimulator<dim>::test_chemical_flow()
// // {
// // 	setup_system();

// // 	std::cout << std::endl << std::endl
// // 		<< "Starting chemical flow test" << std::endl
// // 		<<"--------------------------------------------------------------------"
// // 		<< std::endl << std::endl;

// // 	// save period:
// // 	const unsigned int modsave
// // 		= static_cast<unsigned int>( std::ceil(save_period / chemical_time_step) );

// // 	Point<dim> c = (dim==2) ? Point<dim>(2,2) : Point<dim>(2,2,2);
// // 	double a = 1.0;
// // 	double w = 1.0;

// // 	TestFunction::Gaussian<dim> myfun(c,a,w);

// // 	chemicals.project_function(0, myfun);
// // 	std::vector<Point<dim> > qpoints = geometry.getQuerryPoints();

// // 	do{
// // 		// update time:
// // 		time += chemical_time_step;
// // 		++time_step_number;

// // 		// output:
// // 		if(time_step_number % modsave == 0)
// // 		{
// // 			std::cout << "saving at time: " << time << std::endl;
// // 			output_chemicals(); 
// // 			++save_step_number;
// // 		}

// // 		// chemicals.update(bacteria.getAllLocations(), bacteria.getAllRates());
// // 		chemicals.update();
// // 	}while( (time < run_time) );

// // }


// #endif