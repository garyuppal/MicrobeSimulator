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

// STILL TO ADD:
// -----------------------------------------------
// SPLITTER
// VORTEX AND CYLINDER
// SWISS CHEESE

namespace MicrobeSimulator{ namespace GeometryTools{

// mesh tools: ...
// declare mesh parameters ???
// template<int dim>
// void output_grid(const std::string& output_directory, 
//         const std::string& file_name,
//         const Triangulation<dim>& triangulation)
// {
//   std::string grid_out_file = output_directory + "/" + file_name + ".eps";

//   std::ofstream out (grid_out_file);
//   dealii::GridOut grid_out;
//   grid_out.write_eps (triangulation, out);
//   std::cout << "...Grid written to " << grid_out_file << std::endl;
// }

// /** \brief Declares parameters needed for mesh generation to be read from configuration file
// */
// void declare_parameters(ParameterHandler& prm)
// {
//   prm.enter_subsection("Mesh");
//     prm.declare_entry("Global refinement","0",Patterns::Unsigned());
//     prm.declare_entry("Obstacle refinement","0",Patterns::Unsigned());
//     prm.declare_entry("Boundary refinement","0",Patterns::Unsigned());
//     prm.declare_entry("Mesh type",
//               "Box",
//               Patterns::Selection("Box|Filter|Mixer|File"));
//     prm.declare_entry("Mesh file","",Patterns::Anything());
//   prm.leave_subsection();
// }


/** \brief Base class for building geometry objects and meshes */
/** Creates a local copy of geometry to create grid. 
* Maybe this should also be intialized in constructor ...
* also 
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
	virtual void set_sphere_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria,
                            double sphere_tolerance);
	virtual void set_rectangle_boundary_ids(const Geometry<dim>& geo, Triangulation<dim>& tria);

	// manifolds:
	virtual void attach_mesh_manifolds(const Geometry<dim>& geo, Triangulation<dim>& tria);

	// refinement:
	virtual void refine_obstacles(const Geometry<dim>& geo, Triangulation<dim>& tria);
	virtual void refine_boundary(const Geometry<dim>& geo, Triangulation<dim>& tria);
	virtual void refine_global(Triangulation<dim>& tria);

	void printMeshInfo(std::ostream& out);

protected:
	// mesh parameters*** & boundary conditions
	unsigned int global_refinement;
	unsigned int obstacle_refinement;
	unsigned int boundary_refinement;

	std::array<BoundaryCondition, dim> boundary_conditions; // read in from prm
			// but not always, eg filter is fixed
};

// IMPL
// ---------------------------------------------------------
// mesh construction support methods:
template<int dim>
BuilderBase<dim>::BuilderBase(const ParameterHandler& prm)
{
	const std::string section = "Mesh";
	global_refinement = prm.get_unsigned(section, "Global refinement");
	obstacle_refinement = prm.get_unsigned(section, "Obstacle refinement");
	boundary_refinement = prm.get_unsigned(section, "Boundary refinement");

	// default to periodic boundaries // will be overridden based on geometry type?
	for(unsigned int i = 0; i < dim; ++i)
		boundary_conditions[i] = BoundaryCondition::WRAP;
}

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

template<int dim>
void 
BuilderBase<dim>::set_sphere_boundary_ids(const Geometry<dim>& geo, 
		Triangulation<dim>& tria, double sphere_tolerance)
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

template<int dim>
void 
BuilderBase<dim>::refine_global(Triangulation<dim>& tria)
{
	tria.refine_global(global_refinement);
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
};

// IMPL
// ---------------------------------------------
template<int dim>
Box<dim>::Box(const ParameterHandler& prm)
	:
	BuilderBase<dim>(prm)
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

		this->boundary_conditions[i] = stringToBoundaryCondition(bcs[i]);
	}
}

// class parameters:
template<>
void 
Box<2>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Box");
		prm.declare_entry("Bottom left","{0,0}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Top right","{0,0}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Boundary conditions",
		          "{WRAP,WRAP}",
		          Patterns::List(Patterns::Selection("WRAP|REFLECT|OPEN")),
		          "Boundary conditions.");
	prm.leave_subsection();
}

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
template<int dim>
void 
Box<dim>::build_geometry(Geometry<dim>& geo) const
{
	// maybe make builder a friend? if possible
	geo.setBottomLeftPoint(lower);
	geo.setTopRightPoint(upper);
	geo.setBoundaryConditions(this->boundary_conditions);
}

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
	/** @todo need to make sure grids are labeled and manifolds set */
}

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

	// support methods:
	void construct_filter_side(const double width,
                          const double height,
                          const std::vector< std::vector< double > > &  step_sizes,
                          Triangulation<2>& filter_side) const;
		// THIS DOESNT NEED TO BE A MEMBER FUNCTION...

	void attach_filter_channels(Triangulation<dim>& left_side) const;

	// extrusion to 3D:
	// void extrude(Triangulation<dim>& filter_twodim);
};

// IMPL
// -----------------------------------------------------------------

/** Parameter declaration for filter class */
template<int dim>
void 
Filter<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.enter_subsection("Filter");
			prm.declare_entry("Number channels","1",Patterns::Double());
			prm.declare_entry("Channel thickness","1",Patterns::Double());
			prm.declare_entry("Wall thickness","1",Patterns::Double());
			prm.declare_entry("Left length","1",Patterns::Double());
			prm.declare_entry("Center length","1",Patterns::Double());
			prm.declare_entry("Right length","1",Patterns::Double());
		prm.leave_subsection();
	prm.leave_subsection();
} // maybe separate declarations for 2 and 3 (option to extrude)

/** Constructor for filter class */
template<int dim>
Filter<dim>::Filter(const ParameterHandler& prm)
	:
	BuilderBase<dim>(prm)
{
	const std::string subsection = "Geometry.Filter";
	number_channels = prm.get_unsigned(subsection, "Number channels");
    channel_thickness = prm.get_double(subsection, "Channel thickness");
	wall_thickness = prm.get_double(subsection, "Wall thickness");
	left_length = prm.get_double(subsection, "Left length");
	center_length = prm.get_double(subsection, "Center length");
	right_length = prm.get_double(subsection, "Right length");
}

/** Build geometry object for filter geometry */
template<int dim>
void 
Filter<dim>::build_geometry(Geometry<dim>& geo) const
{
	error, check boundaries;
	
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
}

/** Build mesh for filter geometry */
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


/** Output information for filter class */
template<int dim>
void 
Filter<dim>::printInfo(std::ostream& out) const
{		out << Utility::short_line << std::endl
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
	// void constuct_mixer_center(Triangulation<dim>& tria);
	// void generate_half_circle_hole_tile(Triangulation<dim>& tile);
	// void add_mixer_ends(Triangulation<dim>& center, Triangulation<dim>& tria);
};

// IMPL
// -----------------------------------------------------------------

// constructor
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
}

// class parameters:
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
template<int dim>
void 
Mixer<dim>::build_geometry(Geometry<dim>& geo) const
{

}

template<int dim>
void 
Mixer<dim>::build_grid_base(const Geometry<dim>& geo, Triangulation<dim>& tria) const
{

}  // add option to extrude to 3D

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


// -------------------------------------------------------------------------------
// 		SPLITTER:
// -------------------------------------------------------------------------------


// -------------------------------------------------------------------------------
// 		VORTEX / CYLINDRICAL PIPE (extruded vortex):
// -------------------------------------------------------------------------------
// template<int dim>
// class Cylinder : public BuilderBase<dim>{
// 	// constructor
// 	Cylinder(const ParameterHandler& prm);

// 	// class parameters:
// 	static void declare_parameters(ParameterHandler& prm);

// 	// override virtual methods:
// 	void build_geometry(Geometry<dim>& geo) const override; 
// 	void build_grid_base(Triangulation<dim>& tria) const override; 
// 	void printInfo(std::ostream& out) const override;
// private:
// 	double radius;
// 	// for 3D:
// 	double length;
// };

// -------------------------------------------------------------------------------
// 		FROM FILE(S):
// -------------------------------------------------------------------------------








// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------


// -------------------------------------------------------------------------------
// 		HANDLING CLASS FOR GEOMETRY AND GRID CONSTRUCTION:
// -------------------------------------------------------------------------------
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
/** Construct appropriate building class */
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
		std::cout << "Need to implement" << std::endl;
	else if( boost::iequals(geometry_type, "Cylinder") ) // 2d gives vortex
		std::cout << "Need to implement" << std::endl;
	else if( boost::iequals(geometry_type, "File") )
		std::cout << "Need to implement" << std::endl;
	else
		throw std::runtime_error("Invalid geometry type: <" + geometry_type + ">");
}

/** Declare all parameters for all builder classes*/
template<int dim>
void 
GeometryBuilder<dim>::declare_parameters(ParameterHandler& prm)
{
	// declare mesh parameters ??? GridGenerationTools::declare_parameters(prm);

	prm.enter_subsection("Geometry");
		prm.declare_entry("Geometry type",
		          "Box",
		          Patterns::Selection("Box|Filter|Mixer|Splitter|File"),
		          "Geometry type. Options are Box, Filter, Mixer "
		          "Splitter, of File. For File, name of file must also be"
		          "provided in the \"Geometry file\" parameter. ");
	prm.leave_subsection();

	Box<dim>::declare_parameters(prm);  // thess are not subsubsections ...***
	Filter<dim>::declare_parameters(prm);
	Mixer<dim>::declare_parameters(prm);
	// Cylinder<dim>::declare_parameters(prm);
}

// Main Methods:
/** Build geometry */
template<int dim>
void 
GeometryBuilder<dim>::build_geometry(Geometry<dim>& geo) const
{
	builder->build_geometry(geo);
}

/** Build grid */
template<int dim>
void 
GeometryBuilder<dim>::build_grid(const Geometry<dim>& geo, Triangulation<dim>& tria) const
{
	const int sphere_tolerance = -1; // may need to vary based on geometry type

	// build grid:
	builder->build_grid_base(geo, tria);

	// set boundary labels:
	builder->set_edge_boundary_ids(geo, tria);
	builder->set_sphere_boundary_ids(geo, tria, sphere_tolerance); 
	builder->set_rectangle_boundary_ids(geo, tria); 

	// attach manifolds:
	builder->attach_mesh_manifolds(geo, tria);

	// refine:
	builder->refine_global(tria);
	builder->refine_obstacles(geo, tria);
	builder->refine_boundary(geo, tria);
}

/** Display builder (geometry type) info */
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