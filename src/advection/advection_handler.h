#ifndef MICROBESIMULATOR_ADVECTION_HANDLER_H
#define MICROBESIMULATOR_ADVECTION_HANDLER_H

#include <deal.II/grid/tria.h>

// #include "../utility/enum_types.h"
#include "./velocity_functions.h"
#include "./stokes_solver.h"
#include "../geometry/geometry_builder.h"

#include <memory>

// boost library:
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

namespace MicrobeSimulator{ namespace Velocity{
	using dealii::Triangulation;

template<int dim>
class AdvectionHandler{
public:
	AdvectionHandler();

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm,
			const Geometry<dim>& geometry,
			Triangulation<dim>& triangulation,
			const std::string& output_directory); // instead of using simulation tools

	// from interface:
	Tensor<1, dim> value(const Point<dim>& location) const;
	double get_maximum_velocity(double max_coordinate) const;
	std::vector<std::vector<dealii::Tensor<1, dim> > >
		get_fe_velocity_values(unsigned int fe_degree, 
			const FiniteElement<dim>& fe,
			const dealii::DoFHandler<dim>& dof) const;

	// setup functions:
	void setup_stokes(Triangulation<dim>& tria,
					unsigned int stokes_refinement,
					double inlet_velocity,
					unsigned int number_spheres,
					unsigned int number_rectangles,
					unsigned int number_lines,
					const std::string& output_directory);
	void setup_constant(double flow_rate);
	void setup_couette(double shear);
	void setup_square_pipe(double h, double v);
	void setup_rankine_vortex(double gamma, double radius);

	void printInfo(std::ostream& out) const;
private:
	std::shared_ptr<VelocityInterface<dim> >	velocity_function;
};

// IMPLEMENTATION:
// ----------------------------------------------------------------------
template<int dim>
AdvectionHandler<dim>::AdvectionHandler()
{}

// SETUP:
template<int dim>
void 
AdvectionHandler<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Advection"); 
		prm.declare_entry("Velocity type",
						"No flow",
						Patterns::Selection("No flow|Stokes|Square pipe|Cylindrical pipe|"
										"Vortex|Constant|Linear"),
						"Velocity flow type.");
		prm.declare_entry("Maximum velocity", "0", Patterns::Double());
		prm.declare_entry("Vortex radius", "0", Patterns::Double());
		prm.declare_entry("Vortex circulation", "0", Patterns::Double());
		prm.declare_entry("Stokes refinement","0",Patterns::Double());
	prm.leave_subsection();
}


template<int dim>
void 
AdvectionHandler<dim>::init(const ParameterHandler& prm,
			const Geometry<dim>& geometry,
			Triangulation<dim>& triangulation,
			const std::string& output_directory)
{
	std::string section = "Advection";

	std::string velocity_type = prm.get_string(section,"Velocity type");

	if( boost::iequals(velocity_type, "NO FLOW") )
	{
		std::cout << "...setting up no flow " << std::endl;
		
		setup_constant(0); // ''constant flow'' with value 0
	}
	else if( boost::iequals(velocity_type, "STOKES") )
	{
		std::cout << "...setting up stokes velocity" << std::endl;

		setup_stokes(triangulation,
					prm.get_unsigned(section, "Stokes refinement"),
					prm.get_double(section, "Maximum velocity"),
					geometry.getNumberSpheres(), 
					geometry.getNumberRectangles(),
					geometry.getNumberLines(),
					output_directory); //parameters.getOutputDirectory());
	}
	else if( boost::iequals(velocity_type, "CONSTANT") )
	{
		std::cout << "...setting up constant flow " << std::endl << std::endl;
		
		setup_constant(prm.get_double(section, "Maximum velocity")); 
	}
	else if( boost::iequals(velocity_type, "COUETTE") )
	{
		std::cout << "...setting up constant flow " << std::endl << std::endl;
		
		const double shear_rate = prm.get_double(section, "Maximum velocity")/
									(geometry.getTopRightPoint()[1]); // shear = velocity/height
		setup_couette(shear_rate); 	
	}
	else if( boost::iequals(velocity_type, "SQUARE PIPE") )
	{
		std::cout << "...setting up square pipe flow " << std::endl << std::endl;
		
		setup_square_pipe( geometry.getTopRightPoint()[1], // y height from geometry
							prm.get_double(section, "Maximum velocity") ); 	
	}
	else
	{
		throw std::runtime_error("Need to implement other flow types...");
	}

	// add....
	// 	void setup_square_pipe(double h, double v);
	// void setup_rankine_vortex(double gamma, double radius);
}


// FROM INTERFACE:
template<int dim>
Tensor<1, dim> 
AdvectionHandler<dim>::value(const Point<dim>& location) const
{
	return velocity_function->value(location);
}

template<int dim>
double 
AdvectionHandler<dim>::get_maximum_velocity(double max_coordinate) const
{
	return velocity_function->get_maximum_velocity(max_coordinate);
}

template<int dim>
std::vector<std::vector<dealii::Tensor<1, dim> > >
AdvectionHandler<dim>::get_fe_velocity_values(unsigned int fe_degree, 
			const FiniteElement<dim>& fe,
			const dealii::DoFHandler<dim>& dof) const
{
	return velocity_function->get_fe_velocity_values(fe_degree, fe, dof);
}		

// SETUP METHODS:
// create new velocity interface function of appropriate type...
template<int dim>
void 
AdvectionHandler<dim>::setup_stokes(Triangulation<dim>& tria,
							unsigned int stokes_refinement,
							double inlet_velocity,
							unsigned int number_spheres,
							unsigned int number_rectangles,
							unsigned int number_lines,
							const std::string& output_directory)
{

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

	unsigned int line_id = GridGenerationTools::id_line_begin;
	for(unsigned int i = 0; i < number_lines; ++i, ++line_id)
		no_slip_ids.emplace_back(line_id);

	// can also give fe_degree to constructor 
		//.. stick with default of linear elements for now
	velocity_function = std::make_shared<StokesSolver<dim> >(tria, 
		no_slip_ids, stokes_refinement, inlet_velocity, output_directory); 
}

template<int dim>
void 
AdvectionHandler<dim>::setup_constant(double flow_rate)
{
	velocity_function = std::make_shared<Constant<dim> >(flow_rate);
}

template<int dim>
void
AdvectionHandler<dim>::setup_couette(double shear)
{
	velocity_function = std::make_shared<Couette<dim> >(shear);
}

template<int dim>
void 
AdvectionHandler<dim>::setup_square_pipe(double h, double v)
{
	velocity_function = std::make_shared<SquarePipe<dim> >(h,v);
}

template<int dim>
void 
AdvectionHandler<dim>::setup_rankine_vortex(double gamma, double radius)
{
	velocity_function = std::make_shared<RankineVortex<dim> >(gamma, radius);
}

template<int dim>
void 
AdvectionHandler<dim>::printInfo(std::ostream& out) const
{
	out << "Velocity Info: " << std::endl;
	velocity_function->printInfo(out);
}


}} // close namespaces
#endif