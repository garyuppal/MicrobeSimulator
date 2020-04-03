#ifndef MICROBESIMULATOR_GEOMETRY_BUILDER_H
#define MICROBESIMULATOR_GEOMETRY_BUILDER_H

#include <deal.II/grid/tria.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>

#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <memory>

// boost library:
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "../utility/parameter_handler.h"
#include "./geometry.h"

/** @file
* @todo add vortex/cylinder
* @todo add 3D versions of filter, mixer, splitter, and box
* @todo add file type
* @todo add swiss cheese
*/


namespace MicrobeSimulator{ 

// NAMESPACE GRIDGENERATIONTOOLS
// --------------------------------------------------------------------------------------
/** \brief Grid generation namespace */
/** Set of functions to aid mesh grid generation */
namespace GridGenerationTools{
using dealii::Triangulation;

  	/** \brief Output grid to eps file */
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

	/** \brief Declares parameters needed for mesh generation 
	* to be read from configuration file
	*/
	void declare_parameters(ParameterHandler& prm)
	{
	  prm.enter_subsection("Mesh");
	    prm.declare_entry("Global refinement","0",Patterns::Unsigned());
	    prm.declare_entry("Obstacle refinement","0",Patterns::Unsigned());
	    prm.declare_entry("Boundary refinement","0",Patterns::Unsigned());
	    prm.declare_entry("Max cell size","-1",Patterns::Double());
	  prm.leave_subsection();
	}

	// CONSTANTS FOR GRID GENERATION AND BOUNDARY IDENTIFICATION
	// --------------------------------------------------------------------------------------
	/** left boundary id */
	static constexpr unsigned int id_left = 0;
	/** right boundary id */
	static constexpr unsigned int id_right = 1;
	/**  top boundary id */
	static constexpr unsigned int id_top = 2;
	/** bottom boundary id */
	static constexpr unsigned int id_bottom = 3;

	/** For 3d, front boundary id */
	static constexpr unsigned int id_front = 4;
	/** For 3d, back boundary id */
	static constexpr unsigned int id_back = 5;

	// static constexpr unsigned int id_other = 7;

	/** id of first sphere */
	static constexpr unsigned int id_sphere_begin = 10;

	/** id of first rectangle */
	static constexpr unsigned int id_rectangle_begin = 100;

	// arrays for dimension independent access:
	/** lower faces ids */
	static constexpr std::array< unsigned int , 3 > lower_ids = {id_left, id_bottom, id_back};
	/** upper faces ids */
	static constexpr std::array< unsigned int , 3 > upper_ids = {id_right, id_top, id_front};
} // NAMESPACE GRIDGENERATIONTOOLS



// NAMESPACE GEOMETRYTOOLS
// --------------------------------------------------------------------------------------
/** \brief Namespace for geometry builder and its classes */
/** Classes to construct various geometries and meshes */
namespace GeometryTools{

/** \brief Base class for building geometry objects and meshes */
/** This is a pure virtual class. Implemented default methods are given for setting boundary
* labels, manifolds, and grid refinement. These are mainly common to all constructions, 
* but are set as virtual if they need to be overridden. 
* Derived classes then specialize in creating specific geometries.
*/
template<int dim>
class BuilderBase{
public:
	BuilderBase(const ParameterHandler& prm);
	virtual ~BuilderBase() {}

	// main methods for each class
	virtual void build_geometry(Geometry<dim>& geo) const = 0; 
	virtual void build_grid_base(const Geometry<dim>& geo, 
								Triangulation<dim>& tria) const = 0; 
	virtual void printInfo(std::ostream& out) const=0;

	// mesh construction support methods:
	// boundary labels:
	virtual void set_edge_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria);
	virtual void set_sphere_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria);
	virtual void set_rectangle_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria);

	// manifolds:
	virtual void attach_mesh_manifolds(const Geometry<dim>& geo, Triangulation<dim>& tria);

	// refinement:
	virtual void refine_obstacles(const Geometry<dim>& geo, Triangulation<dim>& tria);
	virtual void refine_boundary(const Geometry<dim>& geo, Triangulation<dim>& tria);
	virtual void refine_global(Triangulation<dim>& tria);
	virtual void refine_largest_cells(Triangulation<dim>& tria);

	void printMeshInfo(std::ostream& out);

protected:
	unsigned int global_refinement;
	unsigned int obstacle_refinement;
	unsigned int boundary_refinement;

	double sphere_tolerance;
	double max_cell_size;
};

// IMPL
// ---------------------------------------------------------
// mesh construction support methods:
/** \brief Constuctor for base class for builder classes */
template<int dim>
BuilderBase<dim>::BuilderBase(const ParameterHandler& prm)
	:
	sphere_tolerance(-1),
	max_cell_size(-1)
{
	const std::string section = "Mesh";
	global_refinement = prm.get_unsigned(section, "Global refinement");
	obstacle_refinement = prm.get_unsigned(section, "Obstacle refinement");
	boundary_refinement = prm.get_unsigned(section, "Boundary refinement");
	max_cell_size = prm.get_double(section, "Max cell size");
}

/** \brief Set boundary labels for outer faces */
template<int dim>
void 
BuilderBase<dim>::set_edge_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria)
{
	const double edge_tolerance = 1e-8;
	const Point<dim> bottom_left = geo.getBottomLeftPoint();
	const Point<dim> top_right = geo.getTopRightPoint();

	for (typename Triangulation<dim>::active_cell_iterator
			cell = tria.begin_active();
			cell != tria.end();
			++cell)
		for (unsigned int f=0; f<dealii::GeometryInfo<dim>::faces_per_cell; ++f)
			if (cell->face(f)->at_boundary())
				for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
					if ( std::fabs( cell->face(f)->center()[dim_itr]
							- bottom_left[dim_itr] ) < edge_tolerance )
						cell->face(f)->set_boundary_id(GridGenerationTools::lower_ids[dim_itr]);
					else if ( std::fabs( cell->face(f)->center()[dim_itr]
							- top_right[dim_itr] ) < edge_tolerance )
						cell->face(f)->set_boundary_id(GridGenerationTools::upper_ids[dim_itr]);
} // set_edge_boundary_ids()

/** \brief Set boundary labels for interior spheres */
template<int dim>
void 
BuilderBase<dim>::set_sphere_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria)
{
	std::vector<Sphere<dim> > spheres = geo.getSpheres();

	if( sphere_tolerance < 0 )
		sphere_tolerance = 0.1*dealii::GridTools::minimal_cell_diameter(tria);

	std::cout << "SPHERE TOLERANCE IS SET TO: " << sphere_tolerance << std::endl;

	for (typename Triangulation<dim>::active_cell_iterator
			cell = tria.begin_active();
			cell != tria.end();
			++cell)
	{
		for(unsigned int f=0; f < dealii::GeometryInfo<dim>::faces_per_cell; ++f)
		{
			unsigned int sphere_bid = GridGenerationTools::id_sphere_begin; // start

			for(unsigned int sp = 0; sp < spheres.size(); ++sp)
			{
				const Point<dim> sphere_center = spheres[sp].getCenter();
				const double sphere_radius = spheres[sp].getRadius();

				const double distance_from_center = sphere_center.distance(cell->face(f)->center());

				if( std::fabs(distance_from_center - sphere_radius) < sphere_tolerance )
					cell->face(f)->set_boundary_id(sphere_bid);

			++ sphere_bid;
			} // for spheres
		} // for faces
	} // for cells
} // set_sphere_boundary_ids()

/** \brief Set boundary labels for interior rectangles */
template<int dim>
void 
BuilderBase<dim>::set_rectangle_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria)
{
	const std::vector<HyperRectangle<dim> > rects = geo.getRectangles();

	const double edge_tolerance = 1e-8;

	for (typename Triangulation<dim>::active_cell_iterator
			cell = tria.begin_active();
			cell != tria.end();
			++cell)
	{
		for(unsigned int f=0; f < dealii::GeometryInfo<dim>::faces_per_cell; ++f)
		{
			unsigned int bid_rect = GridGenerationTools::id_rectangle_begin;

			for(unsigned int i = 0; i < rects.size(); ++i)
			{

				const double distance_to_rectangle_border =
				rects[i].distance_from_border(cell->face(f)->center());

				if ( distance_to_rectangle_border < edge_tolerance)
					cell->face(f)->set_boundary_id(bid_rect);
			} // for rectangles

			++bid_rect;

		} // for faces
	} // for cells
} // set rectangle boundary ids

/** \brief Attach manifolds to interior spheres */
/** Manifolds help ensure grid refinement adds verticies in correct locations.
* If refining a sphere, rather than adding a new vertex in midpoint of an edge,
* the point should be added to the midpoint of the arc of the sphere. Deal.II manifolds
* help ensure this is done properly. See the Deal.II documentation on 
* "Manifold description for triangulations" for more information.
*/ 
template<int dim>
void 
BuilderBase<dim>::attach_mesh_manifolds(const Geometry<dim>& geo, Triangulation<dim>& tria)
{
	std::vector<Sphere<dim> > spheres = geo.getSpheres();

	tria.reset_all_manifolds();

	unsigned int man_id = GridGenerationTools::id_sphere_begin;
	for(unsigned int i = 0; i < spheres.size(); i++)
	{
		std::cout << "setting manifold: " << man_id << std::endl;
		tria.set_all_manifold_ids_on_boundary (man_id, man_id);
		dealii::SphericalManifold<dim> sphere_manifold( spheres[i].getCenter() );
		tria.set_manifold(man_id, sphere_manifold);
		++man_id;
	}
}


// REFINEMENT:
// -------------------------

/** \brief Refine cells at boundaries of interior spheres and rectangles */
template<int dim>
void 
BuilderBase<dim>::refine_obstacles(const Geometry<dim>& geo, Triangulation<dim>& triangulation)
{
	for(unsigned int step = 0; step < obstacle_refinement; ++step)
	{
		for(auto cell : triangulation.active_cell_iterators())
		{
			for(unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
			{
				// for each circle:
				unsigned int number_spheres = geo.getNumberSpheres();
				for(unsigned int i = 0; i < number_spheres; ++i)
				{
					const Point<dim> center = geo.getSphereAt(i).getCenter();
					const double radius = geo.getSphereAt(i).getRadius();

					const double distance_from_center = center.distance(cell->vertex(v));
					if (std::fabs(distance_from_center - radius) < 1e-10)
					{
						cell->set_refine_flag();
						// break;
					} // if vertex on circle boundary
				} // for each sphere

				// for each rectangle
				unsigned int number_rectangles = geo.getNumberRectangles();
				for(unsigned int i=0; i < number_rectangles; ++i)
				{
					if(geo.getRectangleAt(i).distance_from_border(cell->vertex(v)) < 1e-8)
					{
						cell->set_refine_flag();
					}
				}

			} // for each vertex
		} // for each cell in mesh
		triangulation.execute_coarsening_and_refinement();
	} // for each refinement step

	std::cout << "...refined obstacles " << obstacle_refinement << " times" << std::endl;
}

/** \brief Refine Boundary of bounding domain */
template<int dim>
void 
BuilderBase<dim>::refine_boundary(const Geometry<dim>& geo, Triangulation<dim>& triangulation)
{
	for(unsigned int step = 0; step < boundary_refinement; ++step)
	{
		for(auto cell : triangulation.active_cell_iterators())
		{
			for(unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
			{
				// for each face/edge:
				for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
				{
					// lower edge/face:
					const double lower_face = geo.getBottomLeftPoint()[dim_itr];
					if( std::fabs(cell->vertex(v)[dim_itr] - lower_face) < 1e-8)
						cell->set_refine_flag();

					// upper edge/face:
					const double upper_face = geo.getTopRightPoint()[dim_itr];
					if( std::fabs(cell->vertex(v)[dim_itr] - upper_face) < 1e-8)
						cell->set_refine_flag();
				} // for each dimension
			} // for each vertex
		} // for each cell in mesh
		triangulation.execute_coarsening_and_refinement();
	} // for each refinement step

	std::cout << "...refined boundary " << boundary_refinement << " times" << std::endl;
}

/** \brief Globally refine mesh */
template<int dim>
void 
BuilderBase<dim>::refine_global(Triangulation<dim>& tria)
{
	tria.refine_global(global_refinement);
}

/** \brief Refine cells larger than given max_cell_size */
/** @todo note, may also want option to coarsen very small cells */
template<int dim>
void 
BuilderBase<dim>::refine_largest_cells(Triangulation<dim>& triangulation)
{
	if(max_cell_size <= 0)
		return; 

	unsigned int n_refine;
	unsigned int iter = 0;

	std::cout << "max cell size: " << max_cell_size << std::endl;
	do
	{
		n_refine = 0; // initially no cells to refine

		for(auto cell : triangulation.active_cell_iterators())
		{

			if(cell->diameter() > max_cell_size)
			{
				// std::cout << "cell diameter: " << cell->diameter() << std::endl;
				++ n_refine;
				cell->set_refine_flag();
			}

		} // for each cell in mesh
		++iter;
		std::cout << "refining: " << n_refine << "cells on iteration: " << iter << std::endl;
		triangulation.execute_coarsening_and_refinement();
	}while(n_refine != 0); // while there are still cells to refine

	std::cout << "...refined boundary " << boundary_refinement << " times" << std::endl;	
}


// print mesh info:

/** \brief Display mesh refinement info */
template<int dim>
void 
BuilderBase<dim>::printMeshInfo(std::ostream& out)
{
	out << Utility::short_line << std::endl
		<< "\tMESH:" << std::endl
		<< Utility::short_line << std::endl
		<< "Global refinement: " << global_refinement << std::endl
		<< "Obstacle refinement: " << obstacle_refinement << std::endl
		<< "Boundary refinement: " << boundary_refinement << std::endl
		<< Utility::short_line << std::endl << std::endl;
}







// -------------------------------------------------------------------------------
// 		RECTANGLE BOX:
// -------------------------------------------------------------------------------
/** \brief Class to help constuct simple box geometry */
template<int dim>
class Box : public BuilderBase<dim>{
public:
	Box(const ParameterHandler& prm);

	// class parameters:
	static void declare_parameters(ParameterHandler& prm);

	// override virtual methods:
	void build_geometry(Geometry<dim>& geo) const override; 
	void build_grid_base(const Geometry<dim>& geo, 
						Triangulation<dim>& tria) const override; 
	void printInfo(std::ostream& out) const override;
private:
	Point<dim> lower;
	Point<dim> upper;

	std::array<BoundaryCondition, dim> boundary_conditions; 

	int n_refineCenter;

	void refine_center(Triangulation<dim>& tria) const;
};

// IMPL
// ---------------------------------------------
/** \brief Constructor for simple box type geometry builder */
template<int dim>
Box<dim>::Box(const ParameterHandler& prm)
	:
	BuilderBase<dim>(prm),
	n_refineCenter(-1) // for debugging hanging nodes
{
	const std::string section = "Geometry.Box";
	std::vector<double> b = prm.get_double_vector(section, "Bottom left");
	std::vector<double> u = prm.get_double_vector(section, "Top right");
	std::vector<std::string> bcs = prm.get_list(section, "Boundary conditions");

	assert(b.size() == dim);
	assert(u.size() == dim);

	for(unsigned int i = 0; i < dim; ++i)
	{
		lower[i] = b[i];
		upper[i] = u[i];

		boundary_conditions[i] = stringToBoundaryCondition(bcs[i]);
	}

	n_refineCenter = prm.get_int("Debug", "RefineBoxCenter");
}

// class parameters:
/** \brief Two dimensional declaration for box class parameters */
template<>
void 
Box<2>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.enter_subsection("Box");
			prm.declare_entry("Bottom left","{0,0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Top right","{0,0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Boundary conditions",
			          "{WRAP,WRAP}",
			          Patterns::List(Patterns::Selection("WRAP|REFLECT|OPEN")),
			          "Boundary conditions.");
		prm.leave_subsection();
	prm.leave_subsection();

	prm.enter_subsection("Debug");
		prm.declare_entry("RefineBoxCenter", "-1", Patterns::Integer());
	prm.leave_subsection();
}

/** \brief Three dimensional declaration for box class parameters */
template<>
void 
Box<3>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.enter_subsection("Box");
			prm.declare_entry("Bottom left","{0,0,0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Top right","{0,0,0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Boundary conditions",
			          "{WRAP,WRAP}",
			          Patterns::List(Patterns::Selection("WRAP|REFLECT|OPEN")),
			          "Boundary conditions.");
		prm.leave_subsection();
	prm.leave_subsection();
}

// override virtual methods:

/** \brief Build geometry for box type*/
template<int dim>
void 
Box<dim>::build_geometry(Geometry<dim>& geo) const
{
	geo.setBottomLeftPoint(lower);
	geo.setTopRightPoint(upper);
	geo.setBoundaryConditions(this->boundary_conditions);
}

/** \brief Build grid for box geometry */
template<int dim>
void 
Box<dim>::build_grid_base(const Geometry<dim>& /* geo */, Triangulation<dim>& tria) const
{
	tria.clear();
	std::vector<unsigned int> repetitions;
	repetitions.reserve(dim);
	std::vector<double> widths;
	widths.reserve(dim);

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		widths.emplace_back( upper[dim_itr] - lower[dim_itr] );

	// find min length:
	const double min_width = *std::min_element(widths.begin(),widths.end());

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
	{
		unsigned int rep = ceil(widths[dim_itr]/min_width);
		repetitions.emplace_back(rep);
	}

	std::cout << "...generating hyper rectangle" << std::endl;
	dealii::GridGenerator::subdivided_hyper_rectangle(tria, repetitions, lower, upper,
			                                        /*colorize*/ false); // relabeling anyways

	// for debugging:
	if(n_refineCenter > 0)
		refine_center(tria);
}

/** \brief Refine center third portion of mesh as many times as given.
* Mainly use to debug effects of hanging nodes or transitions in mesh refinement
*/
template<int dim>
void
Box<dim>::refine_center(Triangulation<dim>& triangulation) const
{
	const double third = (upper[0] - lower[0])/3.0;
	const double left = lower[0] + third;
	const double right = lower[0] + 2.0*third;

	// refine center in x direction
	for(unsigned int i = 0; i < n_refineCenter; ++i)
	{
		for(auto cell : triangulation.active_cell_iterators())
		{
			// mark cells in center third:
			if( (cell->center()[0] >= left) &&
					(cell->center()[0] <= right) )
				cell->set_refine_flag();
		} // for each cell in mesh

		triangulation.execute_coarsening_and_refinement();
	} // for n_refineCenter steps

	std::cout << "...refined center " << n_refineCenter << " times" << std::endl;		
}

/** \brief Print info for box builder */
template<int dim>
void 
Box<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\tBOX:" << std::endl
		<< Utility::short_line << std::endl
		<< "Bottom left: " << lower << std::endl
		<< "Top right: " << upper << std::endl
		<< Utility::short_line << std::endl << std::endl;
}



// -------------------------------------------------------------------------------
// 		FILTER:
// -------------------------------------------------------------------------------
/** \brief Class to help construct filter geometry */
/** @todo Option to extrude to 3D */
template<int dim>
class Filter : public BuilderBase<dim>{
public:
	// constructor
	Filter(const ParameterHandler& prm);

	// class parameters:
	static void declare_parameters(ParameterHandler& prm);

	// override virtual methods:
	void build_geometry(Geometry<dim>& geo) const override; 
	void build_grid_base(const Geometry<dim>& geo, Triangulation<dim>& tria) const override; 
	void printInfo(std::ostream& out) const override;
private:
	unsigned int number_channels;
	double channel_thickness;
	double wall_thickness;
	double left_length;
	double center_length;
	double right_length;

	double fixed_height;

	/** @todo option to fix total height */

	// support methods:
	void construct_filter_side(const double width,
                          const double height,
                          const std::vector< std::vector< double > > &  step_sizes,
                          Triangulation<2>& filter_side) const;
		// THIS DOESNT NEED TO BE A MEMBER FUNCTION...

	void attach_filter_channels(Triangulation<dim>& left_side) const;
	void reassign_channel_thickness(double buffer);
	// extrusion to 3D:
	// void extrude(Triangulation<dim>& filter_twodim);
};

// IMPL
// -----------------------------------------------------------------

/** \brief Parameter declaration for filter class */
template<int dim>
void 
Filter<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.enter_subsection("Filter");
			prm.declare_entry("Number channels","1",Patterns::Double());
			prm.declare_entry("Channel thickness","1",Patterns::Double());
			prm.declare_entry("Wall thickness","1",Patterns::Double());
			prm.declare_entry("Fixed height","-1",Patterns::Double());
			prm.declare_entry("Left length","1",Patterns::Double());
			prm.declare_entry("Center length","1",Patterns::Double());
			prm.declare_entry("Right length","1",Patterns::Double());
		prm.leave_subsection();
	prm.leave_subsection();
} // maybe separate declarations for 2 and 3 (option to extrude)

/** \brief Constructor for filter class */
template<int dim>
Filter<dim>::Filter(const ParameterHandler& prm)
	:
	BuilderBase<dim>(prm),
	fixed_height(-1)
{
	const std::string subsection = "Geometry.Filter";
	number_channels = prm.get_unsigned(subsection, "Number channels");
    channel_thickness = prm.get_double(subsection, "Channel thickness");
	wall_thickness = prm.get_double(subsection, "Wall thickness");
	left_length = prm.get_double(subsection, "Left length");
	center_length = prm.get_double(subsection, "Center length");
	right_length = prm.get_double(subsection, "Right length");

	fixed_height = prm.get_double(subsection, "Fixed height");

	double edge_buffer = prm.get_double("Bacteria", "Edge buffer");

	if(fixed_height > 0)
		reassign_channel_thickness(edge_buffer);
}

/** \brief Calculate needed channel thickness if using a fixed height */
template<int dim>
void
Filter<dim>::reassign_channel_thickness(double buffer)
{
	// check if assigned fixed height is sufficiently large:
	double gap = fixed_height - (number_channels-1.0)*(wall_thickness + 2.0*buffer) - 2.0*buffer;
	if(gap <= 0)
		throw std::runtime_error("Fixed height is too small for given wall thickness, buffer, and number of channels");

	// if ok, find new channel thickness:
	channel_thickness = ( fixed_height - (number_channels-1.0)*wall_thickness )/number_channels;
}

/** \brief Build geometry object for filter geometry */
template<int dim>
void 
Filter<dim>::build_geometry(Geometry<dim>& geo) const
{
	if(dim != 2)
		throw std::runtime_error("filter geometry currently only implemented for dim == 2");

	Point<dim> bottom_left, top_right;

	// check parameters...
	const unsigned int num_rect = number_channels - 1;

	const double width = left_length + center_length + right_length;
	const double height = number_channels*channel_thickness + num_rect*wall_thickness;

	
	// reassign anyways:
	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		bottom_left[dim_itr] = 0;
	
	top_right[0] = width;
	top_right[1] = height;

	geo.setBottomLeftPoint(bottom_left);
	geo.setTopRightPoint(top_right);

	// add rectangles:
	const double x_left = left_length;
	const double x_right = left_length + center_length;

	double y_bottom = channel_thickness;
	double y_top = channel_thickness + wall_thickness;

	for(unsigned int i = 0; i < num_rect; ++i)
	{
		/** @todo add ability to extrude to 3D or include z as a parameter */
		geo.addRectangle(
			HyperRectangle<2>( Point<2>(x_left,y_bottom),Point<2>(x_right,y_top) ));
		y_bottom += (channel_thickness + wall_thickness);
		y_top += (channel_thickness + wall_thickness);
	}

	/** @todo or should we use base?
	* 	geo.setBoundaryConditions(this->boundary_conditions);
	*/

	// SET PROPER BOUNDARY CONDITIONS FOR FILTER:
	geo.setBoundaryCondition(0, BoundaryCondition::OPEN); // OPEN IN x DIRECTION
	for(unsigned int i = 1; i < dim; ++i)
		geo.setBoundaryCondition(i, BoundaryCondition::REFLECT); // REFLECT REST
}

/** \brief Build mesh for filter geometry */
/** @todo make this for 2D only, the 3D version then builds a 2D grid base calling
* this method and extrudes to the final 3D triangulation
*/
template<int dim>
void 
Filter<dim>::build_grid_base(const Geometry<dim>& /* geo */, Triangulation<dim>& tria) const
{
	Triangulation<2> filter_side;

	const double height = ((double)number_channels - 1.)*wall_thickness
	              + ((double)number_channels)*channel_thickness;

	const double max_thickness = std::max(wall_thickness, channel_thickness);
	const unsigned int number_x_divisions = std::max( (int)std::floor(left_length/max_thickness), 1);
	const double x_step = left_length/((double)number_x_divisions);

	// x divisions:
	std::vector<double> left_x_divisions(number_x_divisions, x_step);

	// y divisions:
	std::vector<double> y_divisions;

	for(unsigned int i = 0; i < number_channels-1; ++i)
	{
		y_divisions.push_back(channel_thickness);
		y_divisions.push_back(wall_thickness);
	}
	y_divisions.push_back(channel_thickness);

	const std::vector<std::vector<double> > left_step_sizes = {left_x_divisions, y_divisions}; /// starting from bottom left

	construct_filter_side(left_length, height, left_step_sizes, filter_side);
	attach_filter_channels(filter_side);

	tria.copy_triangulation(filter_side);
	filter_side.clear();

	const unsigned int number_right_x_divisions = std::max( (int)std::floor(right_length/max_thickness), 1);
	const double right_x_step = right_length/((double)number_right_x_divisions);

	// right x divisions:
	std::vector<double> right_x_divisions(number_right_x_divisions, right_x_step);

	const std::vector<std::vector<double> > right_step_sizes = {right_x_divisions, y_divisions}; /// starting from bottom left

	construct_filter_side(right_length, height, right_step_sizes, filter_side);

	// shift right side:
	Tensor<1,2>  shift_vector;
	shift_vector[0] = left_length + center_length;
	shift_vector[1] = 0;
	dealii::GridTools::shift(shift_vector, filter_side);

	// merge:
	dealii::GridGenerator::merge_triangulations(filter_side, tria, tria);
} // build_mesh_base -- FILTER

/** \brief Supports filter grid constuction. Builds left and right sides to which
* channels are attached
*/
template<>
void 
Filter<2>::construct_filter_side(const double width,
                          const double height,
                          const std::vector< std::vector< double > > &  step_sizes,
                          Triangulation<2>& filter_side) const
{
	dealii::GridGenerator::subdivided_hyper_rectangle(filter_side,
                                                    step_sizes,
                                                    Point<2>(0,0),
                                                    Point<2>(width,height));
}

/** \brief Supports filter grid construction. Attaches channels to left side
* constucted by constuct_filter_side()
*/
template<>
void 
Filter<2>::attach_filter_channels(Triangulation<2>& left_side) const
{
	double bottom_y = 0.;
	double top_y = bottom_y + channel_thickness;

	std::vector<unsigned int> repetitions = {1, 1};
	repetitions[0] = std::max( (int)std::floor(center_length/channel_thickness), 1);

	for(unsigned int i = 0; i < number_channels; ++i)
	{
		Triangulation<2> channel;

		dealii::GridGenerator::subdivided_hyper_rectangle(channel,
		                                                  repetitions,
		                                                  Point<2>(left_length, bottom_y),
		                                                  Point<2>(left_length + center_length, top_y));

		dealii::GridGenerator::merge_triangulations(left_side, channel, left_side);

		bottom_y += (channel_thickness+wall_thickness);
		top_y += (channel_thickness+wall_thickness);
	}
}


/** \brief Output information for filter class */
template<int dim>
void 
Filter<dim>::printInfo(std::ostream& out) const
{		
	out << Utility::short_line << std::endl
		<< "\tFILTER:" << std::endl
		<< Utility::short_line << std::endl
		<< "Number channels: " << number_channels << std::endl
		<< "Channel thickness: " << channel_thickness << std::endl
		<< "Wall thickness: " << wall_thickness << std::endl
		<< "Left length: " << left_length << std::endl
		<< "Center length: " << center_length << std::endl
		<< "Right length: " << right_length << std::endl
		<< Utility::short_line << std::endl << std::endl;
}


// -------------------------------------------------------------------------------
// 		MIXER:
// -------------------------------------------------------------------------------
/** \brief Class to construct Mixer type geometry */
/** @todo Add ability to extrude to 3D. Will require geometry bc handling
* against cyllinders for microbes...
*/
template<int dim>
class Mixer : public BuilderBase<dim>{
public:
	// constructor
	Mixer(const ParameterHandler& prm);

	// class parameters:
	static void declare_parameters(ParameterHandler& prm);

	// override virtual methods:
	void build_geometry(Geometry<dim>& geo) const override; 
	void build_grid_base(const Geometry<dim>& geo, Triangulation<dim>& tria) const override; 
	void printInfo(std::ostream& out) const override;
private:
	double left_length;
	double right_length;
	double height;
	double radius;

	// support methods:
	void construct_mixer_center(Triangulation<dim>& tria) const;
	void generate_half_circle_hole_tile(Triangulation<dim>& tile) const;
	void add_mixer_ends(Triangulation<dim>& center, Triangulation<dim>& tria) const;
};

// IMPL
// -----------------------------------------------------------------

// constructor
/** \brief Constuctor for mixer class */
template<int dim>
Mixer<dim>::Mixer(const ParameterHandler& prm)
	:
	BuilderBase<dim>(prm)
{
	const std::string subsection = "Geometry.Mixer";
	left_length = prm.get_double(subsection, "Left length");
	right_length = prm.get_double(subsection, "Right length");
	height = prm.get_double(subsection, "Height");
	radius = prm.get_double(subsection, "Radius");

	const double buffer = 0.5*height - radius;
	this->sphere_tolerance = 0.5*buffer;
}

// class parameters:
/** \brief Declare parameters needed for mixer */
template<int dim>
void 
Mixer<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.enter_subsection("Mixer");
			prm.declare_entry("Left length","1",Patterns::Double());
			prm.declare_entry("Right length","1",Patterns::Double());
			prm.declare_entry("Height","1",Patterns::Double());
			prm.declare_entry("Radius","1",Patterns::Double());
		prm.leave_subsection();
	prm.leave_subsection();
}

// override virtual methods:

/** \brief Override build geometry for Mixer */
/** @todo Add a depth parameter for possible 3D 
*  @todo will also probably need to add cylinder type reflection in geometry class
*/
template<int dim>
void 
Mixer<dim>::build_geometry(Geometry<dim>& geo) const
{
	if(height < 2.*radius)
		throw std::invalid_argument("Mixer height must be greater than sphere diameter");

	const double width = left_length + 2.*radius + right_length;

	Point<dim> bottom_left, top_right;

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		bottom_left[dim_itr] = 0;

	top_right[0] = width;
	top_right[1] = height;

	// set bounding domain:
	geo.setBottomLeftPoint(bottom_left);
	geo.setTopRightPoint(top_right);

	// add spheres:
	/** @ todo would need to add cylinders for 3D case */
	const double center_x = left_length + radius;

	geo.addSphere( Sphere<2>(Point<2>(center_x, 0.) ,radius) );
	geo.addSphere( Sphere<2>(Point<2>(center_x, height) ,radius) );

	// set boundary conditions:
	// SET PROPER BOUNDARY CONDITIONS FOR Mixer:
	geo.setBoundaryCondition(0, BoundaryCondition::OPEN); // OPEN IN x DIRECTION
	for(unsigned int i = 1; i < dim; ++i)
		geo.setBoundaryCondition(i, BoundaryCondition::REFLECT); // REFLECT REST
}

/** \brief Override build grid for Mixer */
/** @todo add option to extrude to 3D */
template<int dim>
void 
Mixer<dim>::build_grid_base(const Geometry<dim>& /* geo */, Triangulation<dim>& tria) const
{
	// /** @todo Change if...throw 's to assertions */
	// if(geo.getNumberSpheres() != 2)
	// 	throw std::runtime_error("Geometry should have two spheres for mixer mesh");
	// if(geo.getSphereAt(0).getCenter()[0] != geo.getSphereAt(1).getCenter()[0])
	// 	throw std::runtime_error("Two spheres should have same x-coordinate for mixer mesh");

	// // bottom left point should be at origin
	// if(geo.getBottomLeftPoint()[0] != 0 || geo.getBottomLeftPoint()[1] != 0)
	// 	throw std::runtime_error("Bottom left corner for mixer mesh should be located at the origin");

	// if(right_length < 0)
	// 	std::runtime_error("Top right point should be located to the right of sphere for mixer mesh");

	// // const double buffer = build_mixer_mesh(left_length, right_length, height, radius, tria);
	// // :::::
	// 	if(height < (2.*radius) )
	// 		throw std::invalid_argument("Mixer height must be larger than diameter");

		Triangulation<2> center;
		construct_mixer_center(center); 
		add_mixer_ends(center, tria);

		/** @todo add here extrusion for 3D implementation */
}  

/** \brief Output mixer info */
template<int dim>
void 
Mixer<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\tMixer:" << std::endl
		<< Utility::short_line << std::endl
		<< "Left length: " << left_length << std::endl
		<< "Right length: " << right_length << std::endl
		<< "Height: " << height << std::endl
		<< "Radius: " << radius << std::endl
		<< Utility::short_line << std::endl << std::endl;
}


// Mesh constuction support methods:

/** \brief Constuct center portion of mixer mesh */
template<int dim>
void 
Mixer<dim>::construct_mixer_center(Triangulation<dim>& center) const
{
	generate_half_circle_hole_tile(center); // generated with corner at origin

	Triangulation<dim> auxillary;
	auxillary.copy_triangulation(center);

	const double angle = 0.5*dealii::numbers::PI;
	dealii::GridTools::rotate(angle,auxillary); // bottom
	dealii::GridTools::rotate(-angle,center);  // top

	Tensor<1,dim> shift_vector;
	shift_vector[0] = height;
	shift_vector[1] = 0;
	// move bottom left corner back to origin
	dealii::GridTools::shift(shift_vector,auxillary);

	shift_vector[0] = 0;
	shift_vector[1] = height;
	// move above bottom tile
	dealii::GridTools::shift(shift_vector,center);

	// merge:
	dealii::GridGenerator::merge_triangulations(auxillary,center,center);
}

/** \brief Constuct tile with semicircle removed to use in cetner portion of 
* mixer mesh 
*/
template<int dim>
void 
Mixer<dim>::generate_half_circle_hole_tile(Triangulation<dim>& hole_tile) const
{
	if(dim != 2)
		throw std::runtime_error("Function not implemented for dim != 2");

	hole_tile.clear();

	const double outer_radius = 0.5*height;

	dealii::GridGenerator::hyper_cube_with_cylindrical_hole(hole_tile,
					                                        radius,
					                                        outer_radius);

	// shift so corner at origin:
	Tensor<1,dim> shift_vector;
	shift_vector[0] = outer_radius;
	shift_vector[1] = outer_radius;
	dealii::GridTools::shift(shift_vector,hole_tile);

	// find cells to remove:
	std::set< typename Triangulation<dim>::active_cell_iterator > cells_to_remove;

	for (const auto cell : hole_tile.active_cell_iterators())
		if (cell->center()[0] < outer_radius)
			cells_to_remove.insert(cell);

	// remove LHS:
	dealii::GridGenerator::create_triangulation_with_removed_cells(hole_tile,
					                                              cells_to_remove,
					                                              hole_tile);
	// shift new corner to origin
	shift_vector[0] = -outer_radius;
	shift_vector[1] = 0;
	dealii::GridTools::shift(shift_vector,hole_tile);
}

/** \brief Add left and right sides to mixer center */
template<int dim>
void 
Mixer<dim>::add_mixer_ends(Triangulation<dim>& center,
		                Triangulation<dim>& tria) const
{
	// center tile will be a `buffer' length wider than radius
	// we adjust left and right lengths accordingly to get desired lengths
	const double buffer = 0.5*height - radius;
	const double left_width = left_length - buffer;
	const double right_width = right_length - buffer;

	// assert modified lengths are positive:
	assert(left_width>0);
	assert(right_width>0);

	Triangulation<dim> side_tile;

	std::vector<unsigned int> repetitions(dim,2);

	repetitions[0] = 1 + std::floor(left_width/height);

	// for dim == 2 -- can extrude this tile for dim == 3 :
	dealii::GridGenerator::subdivided_hyper_rectangle(side_tile,
	                                                repetitions,
	                                                Point<2>(0,0),
	                                                Point<2>(left_width,height));

	// shift center and merge:
	Tensor<1,2> shift_vector;
	shift_vector[0] = left_width;
	shift_vector[1] = 0.;

	dealii::GridTools::shift(shift_vector, center);
	dealii::GridGenerator::merge_triangulations(side_tile,center,tria);

	// generate RIGHT side:
	side_tile.clear();

	repetitions[0] = 1 + std::floor(right_width/height);

	// for dim == 2 -- can extrude this tile for dim == 3 :
	dealii::GridGenerator::subdivided_hyper_rectangle(side_tile,
	                                                repetitions,
	                                                Point<2>(0,0),
	                                                Point<2>(right_width,height));

	shift_vector[0] = left_width + height; // center width = height
	dealii::GridTools::shift(shift_vector, side_tile);

	// final merge:
	dealii::GridGenerator::merge_triangulations(side_tile,tria,tria);
}




// -------------------------------------------------------------------------------
// 		SPLITTER:
// -------------------------------------------------------------------------------
/** \brief Splitter class to construct splitter geometry and grid */
template<int dim>
class Splitter : public BuilderBase<dim>{
public:
	// constructor
	Splitter(const ParameterHandler& prm);

	// class parameters:
	static void declare_parameters(ParameterHandler& prm);

	// override virtual methods:
	void build_geometry(Geometry<dim>& geo) const override; 
	void build_grid_base(const Geometry<dim>& geo, Triangulation<dim>& tria) const override; 
	void printInfo(std::ostream& out) const override;
private:
	double height;
	double left;
	double right;
	double radius;

	// support methods:
	void construct_splitter_center(Triangulation<dim>& tria) const;
	void add_splitter_ends(Triangulation<dim>& center, Triangulation<dim>& tria) const;
};

// IMPL
// -------------------------------------------------------------------------

// constructor
/** \brief Constuctor for Splitter class */
template<int dim>
Splitter<dim>::Splitter(const ParameterHandler& prm)
	:
	BuilderBase<dim>(prm) // gets mesh refinement parameters
{
	const std::string section = "Geometry.Splitter";

	left = prm.get_double(section, "Left length");
	right = prm.get_double(section, "Right length");
	height = prm.get_double(section, "Height");
	radius = prm.get_double(section, "Radius");
}

// class parameters:
/** \brief Declare parameters for Splitter class */
template<int dim>
void 
Splitter<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.enter_subsection("Splitter");
			prm.declare_entry("Left length","1",Patterns::Double());
			prm.declare_entry("Right length","1",Patterns::Double());
			prm.declare_entry("Height","1",Patterns::Double());
			prm.declare_entry("Radius","1",Patterns::Double());
		prm.leave_subsection();
	prm.leave_subsection();
}

// override virtual methods:

/** \brief Override build geometry for Splitter */
template<int dim>
void 
Splitter<dim>::build_geometry(Geometry<dim>& geo) const
{
	if(dim != 2)
		throw std::invalid_argument("Splitter only implemented for dim == 2");

	Point<dim> bottom_left, top_right; 

	const double width = left + 2.*radius + right;

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		bottom_left[dim_itr] = 0;

	top_right[0] = width;
	top_right[1] = height;

	// set bounding domain:
	geo.setBottomLeftPoint(bottom_left);
	geo.setTopRightPoint(top_right);

	// add spheres:
	const double center_x = left + radius;
	const double center_y = 0.5*height;
	geo.addSphere(Sphere<2>(Point<2>(center_x, center_y) ,radius));

	// SET PROPER BOUNDARY CONDITIONS FOR Splitter:
	geo.setBoundaryCondition(0, BoundaryCondition::OPEN); // OPEN IN x DIRECTION
	for(unsigned int i = 1; i < dim; ++i)
		geo.setBoundaryCondition(i, BoundaryCondition::REFLECT); // REFLECT REST
}

/** \brief Override build mesh for splitter */
/** Splitter mesh consists of a rectangle with a hole cut out in the center 
*/
template<int dim>
void 
Splitter<dim>::build_grid_base(const Geometry<dim>& /* geo */, Triangulation<dim>& tria) const
{
  Triangulation<2> center;
  construct_splitter_center(center);
  add_splitter_ends(center, tria);
}

/** \brief Display info for Splitter object */
template<int dim>
void 
Splitter<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\tSplitter:" << std::endl
		<< Utility::short_line << std::endl
		<< "Height: " << height << std::endl
		<< "Left length: " << left << std::endl
		<< "Right length: " << right << std::endl
		<< "Radius: " << radius << std::endl
		<< Utility::short_line << std::endl << std::endl;
}


// support methods:

/** \brief Construct center portion of splitter grid */
template<int dim>
void 
Splitter<dim>::construct_splitter_center(Triangulation<dim>& tria) const
{
    const double outer_radius = 0.5*height;
    const double inner_radius = radius;

    dealii::GridGenerator::hyper_cube_with_cylindrical_hole (tria,
						                                    inner_radius,
						                                    outer_radius);
    // shift corner to origin:
    Tensor<1,2> shift_vector;
    shift_vector[0] = outer_radius;
    shift_vector[1] = outer_radius;

    dealii::GridTools::shift(shift_vector, tria);
}

/** \brief Add left and right ends to splitter center */
template<int dim>
void 
Splitter<dim>::add_splitter_ends(Triangulation<dim>& center, 
									Triangulation<dim>& tria) const
{
	// center tile will be a `buffer' length wider than radius
	// we adjust left and right lengths accordingly to get desired lengths
	const double buffer = 0.5*height - radius;
	const double left_width = left - buffer;
	const double right_width = right - buffer;

	// assert modified lengths are positive:
	assert(left_width>0);
	assert(right_width>0);

	Triangulation<dim> side_tile;

	std::vector<unsigned int> repetitions(dim,2);

	repetitions[0] = 1 + std::floor(left_width/height);

	// for dim == 2 -- can extrude this tile for dim == 3 :
	dealii::GridGenerator::subdivided_hyper_rectangle(side_tile,
	                                                repetitions,
	                                                Point<2>(0,0),
	                                                Point<2>(left_width,height));

	// shift center and merge:
	Tensor<1,2> shift_vector;
	shift_vector[0] = left_width;
	shift_vector[1] = 0.;

	dealii::GridTools::shift(shift_vector, center);
	dealii::GridGenerator::merge_triangulations(side_tile,center,tria);

	// generate RIGHT side:
	side_tile.clear();

	repetitions[0] = 1 + std::floor(right_width/height);

	// for dim == 2 -- can extrude this tile for dim == 3 :
	dealii::GridGenerator::subdivided_hyper_rectangle(side_tile,
	                                                repetitions,
	                                                Point<2>(0,0),
	                                                Point<2>(right_width,height));

	shift_vector[0] = left_width + height; // center width = height
	dealii::GridTools::shift(shift_vector, side_tile);

	// final merge:
	dealii::GridGenerator::merge_triangulations(side_tile,tria,tria);
}




// -------------------------------------------------------------------------------
// 		VORTEX / CYLINDRICAL PIPE (extruded vortex):
// -------------------------------------------------------------------------------
/** \brief Cylinder geometry builder. 2D version gives a circle,
* and 3D a cylindrical pipe */
template<int dim>
class Cylinder : public BuilderBase<dim>{
public:
	// constructor
	Cylinder(const ParameterHandler& prm);

	// class parameters:
	static void declare_parameters(ParameterHandler& prm);

	// override virtual methods:
	void build_geometry(Geometry<dim>& geo) const override; 
	void build_grid_base(const Geometry<dim>& geo, Triangulation<dim>& tria) const override; 
	void printInfo(std::ostream& out) const override;

	// boundary labels:
	// virtual void set_edge_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria);
	// virtual void set_sphere_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria);
	// virtual void set_rectangle_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria);

	// manifolds:
	void attach_mesh_manifolds(const Geometry<dim>& geo, Triangulation<dim>& tria) override;

private:
	double radius;

	// for 3D:
	double length;
};

// IMPL
// --------------------------------------------------------------------

/** \brief Constuctor for cylinder geometry builder */
template<int dim>
Cylinder<dim>::Cylinder(const ParameterHandler& prm)
	:
	BuilderBase<dim>(prm)
{
	const std::string subsection = "Geometry.Cylinder";
	radius = prm.get_double(subsection, "Radius");
	length = prm.get_double(subsection, "Length");
}

// class parameters:
/** \brief Declare parameters for cylinder geometry */
template<int dim>
void 
Cylinder<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.enter_subsection("Cylinder");
			prm.declare_entry("Radius","1",Patterns::Double());
			prm.declare_entry("Length","1",Patterns::Double());
		prm.leave_subsection();
	prm.leave_subsection();
}

// override virtual methods:

/** \brief Build cylinder geometry */
template<int dim>
void 
Cylinder<dim>::build_geometry(Geometry<dim>& geo) const
{
	Point<dim> origin, bl, tr;
	for(unsigned int i = 0; i < dim; ++i)
	{
		origin[i] = 0;
		bl[i] = -radius;
		tr[i] = radius;
		geo.setBoundaryCondition(i, BoundaryCondition::REFLECT);
	}
	if(dim==3)
	{
		bl[0] = -0.5*length;
		tr[0] = 0.5*length; // extruded in x-direction, circle in y-z plane
					// following deal.ii cylinder builder, putting origin in middle
	}

	geo.setBottomLeftPoint(bl);
	geo.setTopRightPoint(tr);

	// add exterior sphere, centered at origin:
	geo.addSphere(Sphere<dim>(origin, radius, ObstacleType::EXTERIOR));
}

/** \brief Build cylinder grid */
/** @todo still need to implement, also may need to override boundary labeling and manifolds */
template<int dim>
void 
Cylinder<dim>::build_grid_base(const Geometry<dim>& /* geo */, 
									Triangulation<dim>& tria) const
{
	dealii::GridGenerator::cylinder(tria, radius, 0.5*length/* half_length */);
}

/** \brief Manifold assignment for cylinder */
/* Override manifold assignment since deal.ii grid builder for cylinder already
* gives proper manifold assignment. 
*/
template<int dim>
void 
Cylinder<dim>::attach_mesh_manifolds(const Geometry<dim>& /* geo */, 
										Triangulation<dim>& /* tria */ )
{}

/** \brief Print info for cylinder geometry builder */
template<int dim>
void 
Cylinder<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\tCYLINDER:" << std::endl
		<< Utility::short_line << std::endl
		<< "Radius: " << radius << std::endl
		<< "Length: " << length << std::endl
		<< Utility::short_line << std::endl << std::endl;
}


// -------------------------------------------------------------------------------
// 		FROM FILE(S):
// -------------------------------------------------------------------------------
// template<int dim>
// class FileGeometry : public BuilderBase<dim>{

// };







// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------


// -------------------------------------------------------------------------------
// 		HANDLING CLASS FOR GEOMETRY AND GRID CONSTRUCTION:
// -------------------------------------------------------------------------------

/** \brief Main class to build geometry and mesh */
template<int dim>
class GeometryBuilder{
public:
	GeometryBuilder(const ParameterHandler& prm);
	static void declare_parameters(ParameterHandler& prm);

	// main methods:
	void build_geometry(Geometry<dim>& geo) const; 
	void build_grid(const Geometry<dim>& geo, Triangulation<dim>& tria) const; 
	void printInfo(std::ostream& out) const;
private:
	std::shared_ptr<BuilderBase<dim> >		builder;
};

// IMPL
// --------------------------------------------------------------------
/** \brief Constuctor: Construct appropriate building class */
template<int dim>
GeometryBuilder<dim>::GeometryBuilder(const ParameterHandler& prm)
{
	const std::string section = "Geometry";
	std::string geometry_type = prm.get_string(section, "Geometry type");

	if( boost::iequals(geometry_type, "Filter"))
		builder = std::make_shared<Filter<dim> >(prm);
	else if( boost::iequals(geometry_type, "Mixer") )
		builder = std::make_shared<Mixer<dim> >(prm);
	else if( boost::iequals(geometry_type, "Box") )
		builder = std::make_shared<Box<dim> >(prm);
	else if( boost::iequals(geometry_type, "Splitter") )
		builder = std::make_shared<Splitter<dim> >(prm);
	else if( boost::iequals(geometry_type, "Cylinder") ) // 2d gives vortex
		builder = std::make_shared<Cylinder<dim> >(prm);
	// else if( boost::iequals(geometry_type, "File") )
	// 	std::cout << "Need to implement" << std::endl;
	else
		throw std::runtime_error("Invalid geometry type: <" + geometry_type + ">");
}

/** \brief Declare all parameters for all builder classes */
template<int dim>
void 
GeometryBuilder<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.declare_entry("Geometry type",
		          "Box",
		          Patterns::Selection("Box|Filter|Mixer|Splitter|File"),
		          "Geometry type. Options are Box, Filter, Mixer "
		          "Splitter, of File. For File, name of file must also be"
		          "provided in the \"Geometry file\" parameter. ");
	prm.leave_subsection();

	GridGenerationTools::declare_parameters(prm); 
	Box<dim>::declare_parameters(prm);  
	Filter<dim>::declare_parameters(prm);
	Mixer<dim>::declare_parameters(prm);
	Splitter<dim>::declare_parameters(prm);
	Cylinder<dim>::declare_parameters(prm);
}

// Main Methods:
/** \brief Build geometry */
template<int dim>
void 
GeometryBuilder<dim>::build_geometry(Geometry<dim>& geo) const
{
	builder->build_geometry(geo);
}

/** \brief Build grid */
template<int dim>
void 
GeometryBuilder<dim>::build_grid(const Geometry<dim>& geo, Triangulation<dim>& tria) const
{
	// const int sphere_tolerance = -1; // may need to vary based on geometry type

	// build grid:
	builder->build_grid_base(geo, tria);

	// set boundary labels:
	builder->set_edge_boundary_ids(geo, tria);
	builder->set_sphere_boundary_ids(geo, tria); 
	builder->set_rectangle_boundary_ids(geo, tria); 

	// attach manifolds:
	builder->attach_mesh_manifolds(geo, tria);

	// refine:
	builder->refine_global(tria);
	builder->refine_obstacles(geo, tria);
	builder->refine_boundary(geo, tria);
	builder->refine_largest_cells(tria);
}

/** \brief Display builder (geometry type) info */
template<int dim>
void 
GeometryBuilder<dim>::printInfo(std::ostream& out) const
{
	builder->printInfo(out); // can do geometry and mesh info
	out << std::endl << std::endl;
	builder->printMeshInfo(out); // inherited from base 
}

}} // close namespace
#endif