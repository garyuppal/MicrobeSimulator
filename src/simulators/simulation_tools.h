#ifndef MICROBESIMULATOR_SIMULATION_TOOLS_H
#define MICROBESIMULATOR_SIMULATION_TOOLS_H

// #include "../bacteria/bacteria.h"
#include "../chemicals/chemical_fe_base.h"
#include "../chemicals/fe_chemical.h"
#include "../geometry/geometry.h"
#include "../advection/velocity_interface.h"
#include "../advection/stokes_solver.h"
#include "../advection/velocity_functions.h"
#include "../advection/advection_handler.h"

// #include "../utility/argparser.h"
#include "../utility/cell_iterator_map.h"
#include "../utility/enum_types.h"

// NEW:
#include "../utility/parameter_handler.h"

#include "./grid_generation_tools.h"

namespace MicrobeSimulator{ namespace SimulationTools{ 

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// FUNCTION DECLARATIONS:
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------

// PARAMETER DECLARATIONS:
// -------------------------------------------------------------------------------
// void declare_simulation_parameters(ParameterHandler& prm);
void declare_common_parameters(ParameterHandler& prm);
// void declare_flow_parameters(ParameterHandler& prm); // moved to velocity_functions.h
// void declare_geometry_parameters(ParameterHandler& prm);
	// void declare_filter_parameters(ParameterHandler& prm);
	// void declare_mixer_parameters(ParameterHandler& prm);
	// void declare_splitter_parameters(ParameterHandler& prm);
void declare_mesh_parameters(ParameterHandler& prm);
// void declare_chemical_parameters(ParameterHandler& prm);
void declare_bacteria_parameters(ParameterHandler& prm);
void declare_fitness_parameters(ParameterHandler& prm);

// SETUP FUNCTIONS:
// -------------------------------------------------------------------------------
// template<int dim>
// void setup_geometry(const ParameterHandler& prm, Geometry<dim> geo);

template<int dim>
void setup_grid(const ParameterHandler& prm, const Geometry<dim>& geo,
					Triangulation<dim>& tria); // done

// template<int dim>
// void setup_velocity(const ParameterHandler& prm,
// 					const Geometry<dim>& geometry,
// 					Velocity::AdvectionHandler<dim>& velocity_function,
// 					Triangulation<dim>& triangulation,
// 					const std::string& output_directory);

// OUTPUT FUNCTIONS:
// -------------------------------------------------------------------------------
template<int dim>
void output_geometry(const std::string& output_directory, 
			const Geometry<dim>& geometry); // done

template<int dim>
void output_grid(const std::string& output_directory, 
			const std::string& file_name, 
			const Triangulation<dim>& triangulation); // done


void output_vector(const std::vector<double>& vec, 
			const std::string& data_name,
			const std::string& output_directory); // done

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// DEFINITIONS:
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------

// PARAMETER DECLARATIONS:

// void declare_simulation_parameters(ParameterHandler& prm)
// {
// 	// declare_common_parameters(prm);
// 	// Velocity::declare_flow_parameters(prm);
// 	// declare_geometry_parameters(prm);
// 	// declare_mesh_parameters(prm);
// 	// declare_chemical_parameters(prm);
// 	// declare_bacteria_parameters(prm);
// 	// declare_fitness_parameters(prm);
// }

void declare_common_parameters(ParameterHandler& prm)
{
	// *** set job id and dimension in command line
	prm.declare_entry("Simulator type",
						"Bacteria",
						Patterns::Selection("Bacteria|IC_AS|Aging|Chemotaxis"),
						"Simulator type.");
	prm.declare_entry("Time step", "1", Patterns::Double());
	prm.declare_entry("Run time","0",Patterns::Double());
	prm.declare_entry("Save period","1",Patterns::Double());
	prm.declare_entry("Output directory","./",Patterns::Anything());
	// prm.declare_entry("Run cycles") // member of specific simulator ...
}

// void declare_flow_parameters(ParameterHandler& prm)
// {
// 	prm.enter_subsection("Advection"); 
// 		prm.declare_entry("Velocity type",
// 						"No flow",
// 						Patterns::Selection("No flow|Stokes|Square pipe|Cylindrical pipe|"
// 										"Vortex|Constant|Linear"),
// 						"Velocity flow type.");

// 		prm.declare_entry("Maximum velocity", "0", Patterns::Double());
// 		prm.declare_entry("Vortex radius", "0", Patterns::Double());
// 		prm.declare_entry("Vortex circulation", "0", Patterns::Double());
// 		prm.declare_entry("Stokes refinement","0",Patterns::Double());
// 	prm.leave_subsection();
// }

// void declare_geometry_parameters(ParameterHandler& prm)
// {
// 	prm.enter_subsection("Geometry");
// 		prm.declare_entry("Geometry type",
// 							"Box",
// 							Patterns::Selection("Box|Filter|Mixer|Splitter|File"),
// 							"Geometry type. Options are Box, Filter, Mixer "
// 							"Splitter, of File. For File, name of file must also be"
// 							"provided in the \"Geometry file\" parameter. ");
// 		prm.declare_entry("Geometry file",
// 							"",
// 							Patterns::Anything(),
// 							"Name of geometry file to initialize geometry from. "
// 							"Mesh type of geometry file should match mesh type "
// 							"given for grid construction.");
// 		prm.declare_entry("Bottom left","{0,0}",Patterns::List(Patterns::Double()));
// 		prm.declare_entry("Top right","{5,5}",Patterns::List(Patterns::Double()));
// 		prm.declare_entry("Boundary conditions",
// 							"{WRAP,WRAP}",
// 							Patterns::List(Patterns::Selection("WRAP|REFLECT|OPEN")),
// 							"Boundary conditions.");
// 		declare_filter_parameters(prm);
// 		declare_mixer_parameters(prm);
// 		declare_splitter_parameters(prm);
// 	prm.leave_subsection();
// }

// void declare_filter_parameters(ParameterHandler& prm)
// {
// 	prm.enter_subsection("Filter");
// 		prm.declare_entry("Number channels","1",Patterns::Double());
// 		prm.declare_entry("Channel thickness","1",Patterns::Double());
// 		prm.declare_entry("Wall thickness","1",Patterns::Double());
// 		prm.declare_entry("Left length","1",Patterns::Double());
// 		prm.declare_entry("Center length","1",Patterns::Double());
// 		prm.declare_entry("Right length","1",Patterns::Double());
// 	prm.leave_subsection();
// }
// void declare_mixer_parameters(ParameterHandler& prm)
// {
// 	prm.enter_subsection("Mixer");
// 		prm.declare_entry("Left length","1",Patterns::Double());
// 		prm.declare_entry("Right length","1",Patterns::Double());
// 		prm.declare_entry("Height","1",Patterns::Double());
// 		prm.declare_entry("Radius","1",Patterns::Double());
// 	prm.leave_subsection();
// }

// void declare_splitter_parameters(ParameterHandler& prm)
// {
// 	prm.enter_subsection("Splitter");
// 		prm.declare_entry("Left length","1",Patterns::Double());
// 		prm.declare_entry("Right length","1",Patterns::Double());
// 		prm.declare_entry("Height","1",Patterns::Double());
// 		prm.declare_entry("Radius","1",Patterns::Double());
// 	prm.leave_subsection();
// }

void declare_mesh_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Mesh");
		prm.declare_entry("Global refinement","0",Patterns::Unsigned());
		prm.declare_entry("Obstacle refinement","0",Patterns::Unsigned());
		prm.declare_entry("Boundary refinement","0",Patterns::Unsigned());
		prm.declare_entry("Mesh type",
							"Box",
							Patterns::Selection("Box|Filter|Mixer|File"));
		prm.declare_entry("Mesh file","",Patterns::Anything());
	prm.leave_subsection();
}

// void declare_chemical_parameters(ParameterHandler& prm)
// {
// 	prm.enter_subsection("Chemicals");
// 		prm.declare_entry("Number chemicals","2",Patterns::Unsigned());
// 		prm.declare_entry("Implementation","{}",
// 						Patterns::List(Patterns::Selection("FE|DG|FDM")),
// 						"Implementation for chemicals. ...");
// 		prm.declare_entry("Diffusion",
// 							"{5,15}",
// 							Patterns::List(Patterns::Double()));
// 		prm.declare_entry("Decay rate",
// 							"{50,15}",
// 							Patterns::List(Patterns::Double()));
// 		prm.declare_entry("Save chemicals", "False",Patterns::Bool());
// 		prm.declare_entry("Time step factor", "1", Patterns::Double());
// 		prm.declare_entry("Viscosity beta","1",Patterns::Double());
// 	prm.leave_subsection();
// }

void declare_bacteria_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Bacteria");
		prm.declare_entry("Number bacteria","100",Patterns::Unsigned());
		prm.declare_entry("Number groups","1",Patterns::Unsigned());
		prm.declare_entry("Diffusion","0.1",Patterns::Double());
		prm.declare_entry("Secretion rate",
							"{100,100}",
							Patterns::List(Patterns::Double()));
		prm.declare_entry("Mutation rate","0",Patterns::Double());
		prm.declare_entry("Mutation strength", "0", Patterns::Double());
		prm.declare_entry("Initial number cheaters","0",Patterns::Unsigned());
		prm.declare_entry("Deterministic number mutate","0",Patterns::Unsigned());
		prm.declare_entry("Deterministic mutate time","0",Patterns::Double());
		prm.declare_entry("Initial locations",
							"{{0,0}}",
							Patterns::List(Patterns::List(Patterns::Double())));
		prm.declare_entry("Reintroducing","False",Patterns::Bool());
		prm.declare_entry("Reintroduction period","0",Patterns::Double());
		// left length, for reintoduction
		prm.declare_entry("Left subdomain length","-1",Patterns::Double());
	prm.leave_subsection();
}

void declare_fitness_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Fitness");
		prm.declare_entry("Fitness type",
							"OR",
							Patterns::Selection("OR|AND|SUM"));
		prm.declare_entry("Chemical fitness","{0,0}",Patterns::List(Patterns::Double()));
		// prm.declare_entry("Chemical harm","{0}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Secretion cost","{0,0}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Chemical saturation","{1,1}",Patterns::List(Patterns::Double()));
	prm.leave_subsection();
}


// SETUP FUNCTIONS:
// -----------------------------------------------------------------------------
// template<int dim>
// void setup_geometry(const ParameterHandler& prm, Geometry<dim> geo)
// {
// 	std::string mesh_type_str = prm.get_string("Mesh type");
// 	MeshType mesh_type = stringToMeshType(mesh_type_str);

// 	// if(mesh_type == MeshType::FILTER)
// 	// 	geometry.initialize(parameters.getFilter());
// 	// else if(mesh_type == MeshType::MIXER)
// 	// 	geometry.initialize(parameters.getMixer());
// 	// else if(mesh_type == MeshType::BOX_MESH)
// 	// 	geometry.initialize(parameters.getPipe()); // need to implement enum type
// 	// else if(mesh_type == MeshType::FILE_MESH)
// 	// 	geometry.initialize(parameters.getGeometryFile(), 
// 	// 					parameters.getMeshFile());
// 	// else if(mesh_type == MeshType::SPLITTER)
// 	// 	geometry.initialize(parameters.getSplitter());
// 	// else
// 	// 	throw std::invalid_argument("Could not initialize; Geometry type not implemented");

// 	// if(dim == 2)
// 	// {
// 	// 	std::array<BoundaryCondition, dim> bcs = {parameters.getXBoundaryBC(), parameters.getYBoundaryBC()};
// 	// 	geometry.setBoundaryConditions(bcs); // want this to be dimension independent...
// 	// }
	
// 	// geometry.printInfo(std::cout);
// }


template<int dim>
void setup_grid(const ParameterHandler& prm, 
					const Geometry<dim>& geometry,
					Triangulation<dim>& triangulation)
{
	const std::string section = "Mesh"; 
	GridGenerationTools::build_mesh_from_geometry(geometry, triangulation);

	triangulation.refine_global(prm.get_unsigned(section, "Global refinement")); //parameters.getGlobalRefinement());

	GridGenerationTools::refine_obstacles(prm.get_unsigned(section, "Obstacle refinement"), //parameters.getSphereRefinement(), 
						geometry, 
						triangulation);

	GridGenerationTools::refine_boundary(prm.get_unsigned(section, "Boundary refinement"), //parameters.getBoundaryRefinement(),
						geometry,
						triangulation);
}

// template<int dim>
// void setup_velocity(const ParameterHandler& prm,
// 					const Geometry<dim>& geometry,
// 					Velocity::AdvectionHandler<dim>& velocity_function,
// 					Triangulation<dim>& triangulation,
// 					const std::string& output_directory)
// {
// 	std::string section = "Advection";

// 	std::string velocity_type_str = prm.get_string(section,"Velocity type");
// 	VelocityType velocity_type = stringToVelocityType(velocity_type_str);//parameters.getVelocityType();

// 	if(velocity_type == VelocityType::STOKES)
// 	{
// 		std::cout << "...setting up stokes velocity" << std::endl;

// 		velocity_function.setup_stokes(triangulation,
// 								prm.get_unsigned(section, "Stokes refinement"),
// 								prm.get_double(section, "Maximum velocity"),
// 								geometry.getNumberSpheres(), 
// 								geometry.getNumberRectangles(),
// 								output_directory); //parameters.getOutputDirectory());
// 	}
// 	else if(velocity_type == VelocityType::NO_FLOW)
// 	{
// 		std::cout << "...setting up no flow " << std::endl;
		
// 		velocity_function.setup_constant(0); // ''constant flow'' with value 0
// 	}
// 	else if(velocity_type == VelocityType::CONSTANT)
// 	{
// 		std::cout << "...setting up constant flow " << std::endl << std::endl;
		
// 		velocity_function.setup_constant(prm.get_double(section, "Maximum velocity")); //parameters.getMaximumVelocity()); // ''constant flow'' with value 0
// 	}
// 	else
// 	{
// 		throw std::runtime_error("Need to implement other flow types...");
// 	}
// }


// OUTPUT FUNCTIONS:
// -----------------------------------------------------------------------------

template<int dim>
void output_geometry(const std::string& output_directory, const Geometry<dim>& geometry)
{
	std::string geo_out_file = output_directory + "/geometryInfo.dat";
	std::ofstream geo_out(geo_out_file);
	geometry.printInfo(geo_out);

	// output boundary and obstacles for easier postprocessing:
	geometry.outputGeometry(output_directory); 
}


template<int dim>
void output_grid(const std::string& output_directory, 
				const std::string& file_name,
				const Triangulation<dim>& triangulation)
{
	std::string grid_out_file = output_directory + "/" + file_name + ".eps";

	std::ofstream out (grid_out_file);
	dealii::GridOut grid_out;
	grid_out.write_eps (triangulation, out);
	std::cout << "...Grid written to " << grid_out_file << std::endl;
}

void output_vector(const std::vector<double>& vec, 
			const std::string& data_name,
			const std::string& output_directory)
{
	std::string outFile = output_directory + "/" + data_name + ".dat";
	std::ofstream out(outFile);

	for(unsigned int i = 0; i < vec.size(); ++i)	
		out << vec[i] << std::endl;
}














// -------------------------------------------------------------------------------
// LEGACY TOOLS:
// -------------------------------------------------------------------------------

// FUNCTION DECLARATIONS:
// template<int dim>
// void setup_geometry(const ArgParser& parameters, 
// 				Geometry<dim>& geometry); // done

// template<int dim> 
// void setup_grid(const ArgParser& parameters, 
// 			const Geometry<dim>& geometry, 
// 			Triangulation<dim>& triangulation);  // done

// template<int dim>
// void setup_velocity(const ArgParser& parameters,
// 					const Geometry<dim>& geometry,
// 					Velocity::StokesSolver<dim>& velocity_function);

// template<int dim>
// void setup_velocity(const ArgParser& parameters,
// 					const Geometry<dim>& geometry,
// 					Velocity::Constant<dim>& velocity_function);

// template<int dim>
// void setup_velocity(const ArgParser& parameters,
// 					const Geometry<dim>& geometry,
// 					Velocity::AdvectionHandler<dim>& velocity_function);

// template<int dim>
// void setup_velocity(const ArgParser& parameters,
// 					const Geometry<dim>& geometry,
// 					Velocity::AdvectionHandler<dim>& velocity_function,
// 					Triangulation<dim>& triangulation);

template<int dim>
void setup_stokes_velocity(unsigned int stokes_refinement,
							double inlet_velocity,
							unsigned int number_spheres,
							unsigned int number_rectangles,
							Velocity::StokesSolver<dim>& stokes_solution,
							const std::string& output_directory); // done

// move these to advection handler class
// template<int dim>
// double getMaxVelocity(const Geometry<dim>& geometry, 
// 						const Velocity::StokesSolver<dim>& velocity_function);

// template<int dim>
// double getMaxVelocity(const Geometry<dim>& geometry, 
// 						const Velocity::Constant<dim>& velocity_function);

// template<int dim, int numchem>
// void setup_bacteria(const ArgParser& parameters, 
// 			const Geometry<dim>& geometry, 
// 			Bacteria<dim, numchem>& bacteria); // done

// template<int dim>
// std::vector<Point<dim> >
// get_bacteria_locations(const ArgParser& parameters, const Geometry<dim>& geometry);

template<int dim>
std::vector<Point<dim> >
get_bacteria_locations(const Geometry<dim>& geometry, unsigned int number_groups,
		double left_length = -1); 


// OUTPUT:
template<int dim, int numchem>
void output_chemicals(const std::array<Chemicals::FE_Chemical<dim>, numchem>& chemicals,
			const std::string& output_directory,
			unsigned int save_step_number); // done

// template<int dim, int numchem>
// void output_bacteria(const Bacteria<dim,numchem>& bacteria,
// 			const std::string& output_directory, 
// 			unsigned int save_step_number); // done

template<int dim, int numchem>
void output_chemicals(const std::array<Chemicals::FE_Chemical<dim>, numchem>& chemicals,
			const std::string& output_directory,
			unsigned int save_step_number,
			unsigned int cycle); // done

// template<int dim, int numchem>
// void output_bacteria(const Bacteria<dim,numchem>& bacteria,
// 			const std::string& output_directory, 
// 			unsigned int save_step_number,
// 			unsigned int cycle); // done



// IMPLEMENTATION:
// ---------------------------------------------------------------------------------------
// template<int dim>
// void setup_geometry(const ArgParser& parameters, 
// 	Geometry<dim>& geometry)
// {
// 	std::cout << "...Setting up geometry" << std::endl;
// 	MeshType mesh_type = parameters.getMeshType();

// 	if(mesh_type == MeshType::FILTER)
// 		geometry.initialize(parameters.getFilter());
// 	else if(mesh_type == MeshType::MIXER)
// 		geometry.initialize(parameters.getMixer());
// 	else if(mesh_type == MeshType::BOX_MESH)
// 		geometry.initialize(parameters.getPipe()); // need to implement enum type
// 	else if(mesh_type == MeshType::FILE_MESH)
// 		geometry.initialize(parameters.getGeometryFile(), 
// 						parameters.getMeshFile());
// 	else if(mesh_type == MeshType::SPLITTER)
// 		geometry.initialize(parameters.getSplitter());
// 	else
// 		throw std::invalid_argument("Could not initialize; Geometry type not implemented");

// 	if(dim == 2)
// 	{
// 		std::array<BoundaryCondition, dim> bcs = {parameters.getXBoundaryBC(), parameters.getYBoundaryBC()};
// 		geometry.setBoundaryConditions(bcs); // want this to be dimension independent...
// 	}
	
// 	geometry.printInfo(std::cout);
// }

// template<int dim> 
// void setup_grid(const ArgParser& parameters, 
// 	const Geometry<dim>& geometry, 
// 	Triangulation<dim>& triangulation)
// {
// 	std::cout << "...Setting up grid" << std::endl;

// 	GridGenerationTools::build_mesh_from_geometry(geometry, triangulation);

// 	// const double max_diameter = GridTools::maximal_cell_diameter(triangulation);
// 	// std::cout << "\n\t ...initial max cell diameter: " << max_diameter << std::endl;

// 	triangulation.refine_global(parameters.getGlobalRefinement());

// 	GridGenerationTools::refine_obstacles(parameters.getSphereRefinement(), 
// 						geometry, 
// 						triangulation);

// 	GridGenerationTools::refine_boundary(parameters.getBoundaryRefinement(),
// 						geometry,
// 						triangulation);
// }


// // POLYMORPHIC FUNCTIONS:
// template<int dim>
// void setup_velocity(const ArgParser& parameters,
// 					const Geometry<dim>& geometry,
// 					Velocity::StokesSolver<dim>& velocity_function)
// {
// 	setup_stokes_velocity(parameters.getStokesRefinement(),
// 							parameters.getMaximumVelocity(),
// 							geometry.getNumberSpheres(),
// 							geometry.getNumberRectangles(),
// 							velocity_function,
// 							parameters.getOutputDirectory());
// }

// template<int dim>
// void setup_velocity(const ArgParser& parameters,
// 					const Geometry<dim>& geometry,
// 					Velocity::Constant<dim>& velocity_function)
// {
// 	velocity_function.setFlowRate(parameters.getMaximumVelocity());
// }

// template<int dim>
// void setup_velocity(const ArgParser& parameters,
// 					const Geometry<dim>& geometry,
// 					Velocity::AdvectionHandler<dim>& velocity_function)
// {
// 	velocity_function.setup_constant(parameters.getMaximumVelocity()); 
// }

// // use this to set up advection handler ... get type to setup from parameters...
// template<int dim>
// void setup_velocity(const ArgParser& parameters,
// 					const Geometry<dim>& geometry,
// 					Velocity::AdvectionHandler<dim>& velocity_function,
// 					Triangulation<dim>& triangulation)
// {
// 	VelocityType velocity_type = parameters.getVelocityType();

// 	if(velocity_type == VelocityType::STOKES)
// 	{
// 		std::cout << "...setting up stokes velocity" << std::endl;

// 		velocity_function.setup_stokes(triangulation,
// 								parameters.getStokesRefinement(),
// 								parameters.getMaximumVelocity(),
// 								geometry.getNumberSpheres(),
// 								geometry.getNumberRectangles(),
// 								parameters.getOutputDirectory());
// 	}
// 	else if(velocity_type == VelocityType::NO_FLOW)
// 	{
// 		std::cout << "...setting up no flow " << std::endl;
		
// 		velocity_function.setup_constant(0); // ''constant flow'' with value 0
// 	}
// 	else if(velocity_type == VelocityType::CONSTANT)
// 	{
// 		std::cout << "...setting up constant flow " << std::endl << std::endl;
		
// 		velocity_function.setup_constant(parameters.getMaximumVelocity()); // ''constant flow'' with value 0
// 	}
// 	else
// 	{
// 		throw std::runtime_error("Need to implement other flow types...");
// 	}

// }

// SPECIFIC IMPLEMENTATIONS:

template<int dim>
void setup_stokes_velocity(unsigned int stokes_refinement,
							double inlet_velocity,
							unsigned int number_spheres,
							unsigned int number_rectangles,
							Velocity::StokesSolver<dim>& stokes_solution,
							const std::string& output_directory)
{
	if(inlet_velocity < 1e-6)
		return;
	
	std::cout << "...Setting up stokes velocity" << std::endl;

	// get boundary conditions:
	std::vector<unsigned int> no_slip_ids = {2,3}; // always no slip top and bottom walls
	if(dim == 3)
	{
		no_slip_ids.emplace_back(4);
		no_slip_ids.emplace_back(5);
	}

	unsigned int sphere_id = GridGenerationTools::id_sphere_begin;
	for(unsigned int i = 0; i < number_spheres; ++i, ++sphere_id)
		no_slip_ids.emplace_back(sphere_id);

	unsigned int rect_id = GridGenerationTools::id_rectangle_begin;
	for(unsigned int i = 0; i < number_rectangles; ++i, ++rect_id)
		no_slip_ids.emplace_back(rect_id);

	// set boundaries (taking 0 to be inlet, and 1 to be outlet):
	stokes_solution.setNoSlipIDs(no_slip_ids);

	// set number of refinement cycles:
	stokes_solution.setNumberRefinementSteps(stokes_refinement);

	// solve and output:
	stokes_solution.solve(inlet_velocity, output_directory);
}

// template<int dim>
// double getMaxVelocity(const Geometry<dim>& geometry, 
// 						const Velocity::StokesSolver<dim>& velocity_function)
// {
// 	return velocity_function.get_maximal_velocity();
// }

// template<int dim>
// double getMaxVelocity(const Geometry<dim>& geometry, 
// 						const Velocity::Constant<dim>& velocity_function)
// {
// 	return velocity_function.getFlowRate();
// }


// BACTERIA SETUP:
// ----------------------------------------------------------------------
// template<int dim, int numchem>
// void setup_bacteria(const ArgParser& parameters, 
// 	const Geometry<dim>& geometry, 
// 	Bacteria<dim, numchem>& bacteria)
// {
// 	std::cout << "... setting up bacteria" << std::endl;

// 	if(numchem != 2)
// 		std::cout << "Setup bacteria function for 2 chem only for now" << std::endl;

// 	std::vector<Point<2> > initial_locations = parameters.getInitialLocations();

// 	if(initial_locations.empty())
// 	{
// 		unsigned int number_groups = parameters.getNumberGroups();
// 		if(number_groups == 0)
// 			number_groups = parameters.getNumberBacteria();
// 		const double left_length = parameters.getLeftSubdomainLength();
// 		initial_locations = get_bacteria_locations(geometry, number_groups, left_length);
// 	}

// 	bacteria.init(parameters.getBacteriaDiffusionConstant(),
// 		parameters.getNumberBacteria(),
// 		std::array<double, 2>(
// 		{parameters.getGoodSecretionRate(), parameters.getWasteSecretionRate()}),
// 		initial_locations); 

// 	bacteria.printInfo(std::cout);
// }

// template<int dim>
// std::vector<Point<dim> >
// get_bacteria_locations(const ArgParser& parameters, const Geometry<dim>& geometry)
// {
// 	const unsigned int number_bacteria = parameters.getNumberBacteria();
// 	unsigned int number_groups = parameters.getNumberGroups();

// 	std::cout << "...Finding " << number_groups
// 		<< " group positions" << std::endl;
	
// 	if(number_groups == 0)
// 		number_groups = number_bacteria; // no groups same as NB ``groups''

// 	std::vector<Point<dim> > group_locations;
// 	group_locations.reserve(number_groups);

// 	Point<dim> temp_point;
// 	for(unsigned int i = 0; i < number_groups; ++i)
// 	{
// 		bool found = false;

// 		while(!found)
// 		{
// 		  for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
// 		  {
// 		  	double width = geometry.getWidth(dim_itr);

// 		  	/// possibly initialize in subdomain:
// 		  	// if( (dim_itr == 0) && (left_length > 0) )
// 		  	// 	width = left_length; 

// 		    temp_point[dim_itr] = (width)*((double)rand() / RAND_MAX) 
// 		      + geometry.getBottomLeftPoint()[dim_itr];
// 		  }

// 		  if( geometry.isInDomain(temp_point) )
// 		  {
// 		    group_locations.emplace_back(temp_point);
// 		    found = true;
// 		  }
// 		} // while not found
// 	} // for group locations
// 	std::cout << "...Group positions found." << std::endl;

// 	return group_locations;
// }

template<int dim>
std::vector<Point<dim> >
get_bacteria_locations(const Geometry<dim>& geometry, unsigned int number_groups,
		double left_length = -1)
{
	std::cout << "...Finding " << number_groups
		<< " group positions" << std::endl;

	std::vector<Point<dim> > group_locations;
	group_locations.reserve(number_groups);

	Point<dim> temp_point;
	for(unsigned int i = 0; i < number_groups; ++i)
	{
		bool found = false;

		while(!found)
		{
		  for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		  {
		  	double width = geometry.getWidth(dim_itr);

		  	// possibly initialize in subdomain:
		  	if( (dim_itr == 0) && (left_length > 0) )
		  		width = left_length; 

		    temp_point[dim_itr] = (width)*((double)rand() / RAND_MAX) 
		      + geometry.getBottomLeftPoint()[dim_itr];
		  }

		  if( geometry.isInDomain(temp_point) )
		  {
		    group_locations.emplace_back(temp_point);
		    found = true;
		  }
		} // while not found
	} // for group locations
	std::cout << "...Group positions found." << std::endl;

	return group_locations;
}


// OUTPUT:
// ------------------------------------------------------------------------------------


template<int dim, int numchem>
void output_chemicals(const std::array<Chemicals::FE_Chemical<dim>, numchem>& chemicals,
	const std::string& output_directory,
	unsigned int save_step_number)
{
	for(unsigned int c = 0; c < numchem; ++c)
		chemicals[c].output_solution(output_directory, c, save_step_number);
}

// template<int dim, int numchem>
// void output_bacteria(const Bacteria<dim,numchem>& bacteria,
// 	const std::string& output_directory, 
// 	unsigned int save_step_number)
// {
// 	std::string outfile = output_directory
// 							+ "/bacteria_" 
// 							+ dealii::Utilities::int_to_string(save_step_number,4)
// 							+ ".dat";
// 	std::ofstream out(outfile);
// 	bacteria.print(out);
// }

template<int dim, int numchem>
void output_chemicals(const std::array<Chemicals::FE_Chemical<dim>, numchem>& chemicals,
	const std::string& output_directory,
	unsigned int save_step_number,
	unsigned int cycle)
{
	for(unsigned int c = 0; c < numchem; ++c)
		chemicals[c].output_solution(output_directory, c, save_step_number, cycle);
}

// template<int dim, int numchem>
// void output_bacteria(const Bacteria<dim,numchem>& bacteria,
// 	const std::string& output_directory, 
// 	unsigned int save_step_number,
// 	unsigned int cycle)
// {
// 	std::string outfile = output_directory
// 							+ "/bacteria_" 
// 							+ dealii::Utilities::int_to_string(cycle,2)
// 	              			+ "-"
// 							+ dealii::Utilities::int_to_string(save_step_number,4) // cycle
// 							+ ".dat";
// 	std::ofstream out(outfile);
// 	bacteria.print(out);
// }



}} // close namespace
#endif