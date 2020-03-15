#ifndef MICROBE_SIMULATOR_GRID_GENERATION_TOOLS
#define MICROBE_SIMULATOR_GRID_GENERATION_TOOLS

#include <deal.II/grid/tria.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/grid/manifold_lib.h>

#include "../geometry/geometry.h"
#include "../geometry/sphere.h"

#include <vector>
#include <functional>
#include <algorithm>

namespace MicrobeSimulator{
    /** \brief Grid generation namespace
    */
    /**
    *  Set of functions to aid mesh grid generation
    */
    namespace GridGenerationTools{
  using dealii::Triangulation;

/** \brief Declares parameters needed for mesh generation to be read from configuration file
*/
void declare_parameters(ParameterHandler& prm)
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

// CONSTANTS FOR GRID GENERATION AND BOUNDARY IDENTIFICATION
// --------------------------------------------------------------------------------------
/// left boundary id
  const unsigned int id_left = 0;
/// right boundary id
  const unsigned int id_right = 1;
/// top boundary id
  const unsigned int id_top = 2;
/// bottom boundary id
  const unsigned int id_bottom = 3;

/// For 3d, front boundary id
  const unsigned int id_front = 4;
/// For 3d, back boundary id
  const unsigned int id_back = 5;

  const unsigned int id_other = 7;

  /// first sphere:
  const unsigned int id_sphere_begin = 10;

  /// first rectangle:
  const unsigned int id_rectangle_begin = 100;

  const std::array< unsigned int , 3 > lower_ids = {id_left, id_bottom, id_back};
  const std::array< unsigned int , 3 > upper_ids = {id_right, id_top, id_front};


/** Squish Transform
* used to help create hex_cheese meshes
*/
Point<2> squish_transform (const Point<2> &in, double tile_width)
{

  double difference = (1.0 - 0.5*std::sqrt(3))*tile_width;
  double shift = 0.5*difference;

  double y_value = in(1);
  if(in(1) == 0)
    y_value = 0 + shift;
  else if(in(1) == tile_width)
    y_value = tile_width - shift;

  return Point<2>(in(0),
                  y_value);
}

// PROTOTYPES:
// ========================================================================

/** \brief Mesh generation from geometry
*/
/**
* In general, the geometry class will store the type of mesh desired
* including a possible mesh file to read in
* grid generation can then be given by calling the build_mesh_from_geometry
* function. This will then call the necessary functions to construct the mesh
*/
template<int dim>
void build_mesh_from_geometry(const Geometry<dim>& geo,
                              Triangulation<dim>& tria);

/** \brief Mesh from file
*/
/**
* Import mesh from file, add boundary labels and manifolds if needed
*/
template<int dim>
void build_mesh_from_file(const Geometry<dim>& geo,
                          Triangulation<dim>& tria);


/** \brief Simple pipe mesh
*/
/**
* Construct hyper rectangle from geometry end points
*/
template<int dim>
void build_hyper_rectangle(const Geometry<dim>& geo,
                            Triangulation<dim>& tria);


/** \brief Mixer
*/
/**
* This function takes four arguments to specify the lengths of the
* left end, the right end, the height, and the radius of the hemispheres used.
* The current implementation requires left and right lenghts be greater than the height.
* We also require the radius to be less than height.
* For a more uniform mesh, probably best to choose
* left = right = height + (height - radius).
* It is up to the user to make sure the geometry discription for bacteria match
* this mesh.
* @todo Could perhaps test to see if geometry agrees with mesh and warn
* user if not. Implement 3D using extrusion.
*/
template<int dim>
double build_mixer_mesh(const double left_length,
                      const double right_length,
                      const double height,
                      const double radius,
                      Triangulation<dim>& tria,
                      const bool relabel = false);

    // helpers:
    /** \brief Helper function to constuct center of mixer mesh.
    */
    /** Constucts the center portion of the mixer mesh, using the height
    * and radius of the hemispheres used.
    */
    template<int dim>
    void construct_mixer_center(const double height,
                                const double radius,
                                Triangulation<dim>& tria);
    /** \brief Helper function to generate a tile with half a sphere cut out.
    */
    /** This method is used by construct_mixer_center() to stack to hemispheres
    * to build the center portion of the mixer mesh.
    */
    template<int dim>
    void generate_half_circle_hole_tile(const double height,
                                        const double radius,
                                        Triangulation<dim>& tile);
    // width is calculated in build_mixer_mesh, width = left - (height - radius)

    /** \brief Helper function to attach left and right ends to mixer center
    */
    template<int dim>
    void add_mixer_ends(const double left_width,
                        const double right_width,
                        const double height,
                        Triangulation<dim>& center,
                        Triangulation<dim>& tria);

    /** \brief Geometry wrapper for mixer mesh
    */
    template<int dim>
    void build_mixer_mesh(const Geometry<dim>& geo,
                          Triangulation<dim>& tria);



// FILTERS:
// -----------------------------------------------------------------------------------

/** \brief Function to build filter type mesh.
*/
/**
* Build filter mesh takes 6 parameters.
* A wrapper is also available to use with geometry class.
* Here, the height of the system is computed as
* (number of channels - 1)x(wall_thickness) +
* (number channels)xchannel_thickness.
* The mesh is constructed by attaching channels to left side and then adding
* the right side.
* @todo currently implemented only for 2d, again, can extrude for 3d
*/
template<int dim>
void build_filter_mesh(const double left_length,
                      const double center_filter_length,
                      const double right_length,
                      const unsigned int number_channels,
                      const double wall_thickness,
                      const double channel_thickness,
                      Triangulation<dim>& tria,
                      const bool relabel = false);

    // helper functions:
    /** \brief Builds sides of filter mesh
    */
    /**
    * Builds side so that vertices of mesh match placement of filter channels
    */
    void construct_filter_side(const double width,
                              const double height,
                              const std::vector< std::vector< double > > &  step_sizes,
                              Triangulation<2>& filter_side);

    /** \brief Attaches filter channels to filter side tile.
    */
    void attach_fitler_channels(const double left_length,
                                const double center_filter_length,
                                const unsigned int number_channels,
                                const double wall_thickness,
                                const double channel_thickness,
                                Triangulation<2>& left_side);

    /** \brief geometry wrapper for filter construction
    */
    template<int dim>
    void build_filter_mesh(const Geometry<dim>& geo,
                          Triangulation<dim>& tria);

// @todo look at subdivided hyper rectangle with option to remove cells,
// perhaps we can use this and attach spheres to get porous mesh...


// Boundary labeling and manifold attachment from geometry:
// -----------------------------------------------------------------------

  /** \brief Wrapper to call with geometry object
  */
  /**
  * Assigns IDs to boundary faces and interior manifolds to allow
  * further refinement of proper shape and to deal with boundary conditions.
  * We attach spherical manifolds to the mesh so that refinements add vertices to
  * the proper locations. See deal.ii documentation on manifolds for more details.
  */
  template<int dim>
  void set_boundary_and_manifolds(const Geometry<dim>& geo,
                                  Triangulation<dim>& tria,
                                  double sphere_tolerance = -1);

// SPLITTERS:
// ---------------------------------------------------------------------------
/** \brief Splitter mesh with hole in center
*/
/**
* Built similar to mixer mesh, consists of left end, a center
* with a spherical obstacle and a right end.
*/
template<int dim>
void build_splitter_mesh(const double left_length,
                        const double right_length,
                        const double height,
                        const double radius,
                        Triangulation<dim>& tria);

    // helpers:
    /** \brief Helper function to constuct center portion of splitter mesh
    */
    template<int dim>
    double construct_splitter_center(const double height,
                                const double radius,
                                Triangulation<dim>& tria);

    /** \brief Wrapper to call with geometry object
    */
    template<int dim>
    void build_splitter_mesh(const Geometry<dim>& geo,
                          Triangulation<dim>& tria);

    // General functions to set boundaries and manifolds:
    // --------------------------------------------------------------------------------

    /** \brief Sets boundary IDs for outer faces of mesh
    */
    template<int dim>
    void set_edge_boundary_ids(const Point<dim>& bottom_left,
                        const Point<dim>& top_right,
                        Triangulation<dim>& tria);

    /** \brief Set IDs to inner spherical boundaries
    */
    template<int dim>
    void set_sphere_boundary_ids(const std::vector<Sphere<dim> >& spheres,
                        Triangulation<dim>& tria,
                        double sphere_tolerance = -1);

    /** \brief Set IDs to inner rectangular bounaries
    */
    template<int dim>
    void set_rectangle_boundary_ids(const std::vector<HyperRectangle<dim> >& rects,
                        Triangulation<dim>& tria);

    /** \brief Set any manifolds to mesh
    */
    /**
    * Manifolds help to deal with further possible refinement. See deal.ii documentation
    * on manifolds for further details.
    */
    template<int dim>
    void attach_mesh_manifolds(const std::vector<Sphere<dim> >& spheres,
                              Triangulation<dim>& tria);


// Refinement methods:
// -----------------------------------------------------------------------------------
/** \brief Refine boundaries of interior obstacles
*/
template<int dim>
void refine_obstacles(unsigned int sphere_refinement,
            const Geometry<dim>& geometry,
            Triangulation<dim>& triangulation);

/** \brief Refine outer boundary of mesh
*/
template<int dim>
void refine_boundary(unsigned int boundary_refinement,
            const Geometry<dim>& geometry,
            Triangulation<dim>& triangulation);





// LEGACY FUNCTIONS ...
// ----------------------------------------------------------------------------------------

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void grid_hole_cutouts(const Geometry<dim>& geo,
                      Triangulation<dim>& tria);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void refine_along_spheres(const Geometry<dim> & geo,
                        Triangulation<dim>& tria);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void cutout_spheres(const Geometry<dim> & geo,
                        Triangulation<dim>& tria);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void generateGrid(const Geometry<dim>& geo,
                  Triangulation<dim>& triangulation,
                  unsigned int sphere_refinement = 0);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void generateRectangleWithHole(double height, Triangulation<2,2>& triangulation);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void generateGridWithHoles(const Geometry<dim>& geo,
                          Triangulation<dim>& triangulation);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void square_pipe_with_spherical_hole(Triangulation<dim>& triangulation, double scale_factor = 1);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void set_pipe_edges(Triangulation<dim>& triangulation, double scale_factor = 1);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void set_pipe_manifolds(Triangulation<dim>& triangulation, double scale_factor = 1);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void cube_with_spherical_hole(Triangulation<dim>& triangulation);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void importMeshFromFile(const Geometry<dim>& geo, Triangulation<dim>& triangulation);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void simpleRectangle(Triangulation<dim> &tria,
                     const Geometry<dim>& geo,
                     const bool colorize=false);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void swissCheese(const Geometry<dim>& geo, Triangulation<dim> &tria);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void squished_square_with_hole(double hole_radius, double long_side,
	Triangulation<dim>& triangulation);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void swissRowsToFullMesh(double tile_height, double full_height,
    Triangulation<dim>& bottom_row, Triangulation<dim>& second_row,
    Triangulation<dim>& tria);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void squished_half_circle_hole(double hole_radius, double tile_width,
	Triangulation<dim>& triangulation);

template<int dim>
void generateHalfCircleHole(double hole_radius, double tile_width,
	Triangulation<dim>& triangulation, unsigned int sphere_refinement = 0);

template<int dim>
void constructSwissBottomRow(double tile_width, double full_width,
	 Triangulation<dim>& hole_tile, Triangulation<dim>& bottom_row);

template<int dim>
void constructSwissSecondRow(double tile_width, double full_width,
  	Triangulation<dim>& half_circle, Triangulation<dim>&hole_tile,
  	Triangulation<dim>& second_row);

template<int dim>
void buildMixerMesh(const Geometry<dim>& geo,
    Triangulation<dim>& triangulation, unsigned int sphere_refinement);

template<int dim>
void constructMixerCenter(double tile_width,
  Triangulation<dim>& hole_tile,
  Triangulation<dim>& triangulation);

template<int dim>
void addMixerEnds(double tile_width, double full_width,
  Triangulation<dim>& triangulation,
  unsigned int sphere_refinement);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void querryBoundaryFaces(const Triangulation<dim>& triangulation);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void setEdgeBoundaries(double left, double right,
	double top, double bottom, Triangulation<dim>& tria);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void setEdgeBoundaries(const Geometry<dim>& geo, Triangulation<dim>& triangulation);

/** \brief Legacy function */
/** @ todo: remove dependencies and delete function */
template<int dim>
void attach_manifolds(const Geometry<dim>& geo, Triangulation<dim>& triangulation);










// IMPLEMENTATION:
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

/** MESH CONSTRUCTION FROM GEOMETRY:
*/

// template<int dim>
// void build_mesh_from_geometry(const Geometry<dim>& geo,
//                               Triangulation<dim>& tria)
// {
//   if(dim == 2)
//   {
//     if(geo.getMeshType() == MeshType::FILE_MESH)
//     {
//       build_mesh_from_file(geo, tria);
//       attach_mesh_manifolds(geo.getSpheres(), tria); //  ids already set in gmsh
//     }
//     // else if(geo.getMeshType() == MeshType::BOX_MESH)
//     //   simpleRectangle(triangulation,geo,true);
//     // else if(geo.getMeshType() == MeshType::SQUARE_CHEESE)
//     //   std::cout << "need to implement" << std::endl;
//     // else if(geo.getMeshType() == MeshType::HEX_CHEESE)
//     //   swissCheese(geo, triangulation);
//     else if(geo.getMeshType() == MeshType::MIXER)
//       build_mixer_mesh(geo, tria);
//     else if(geo.getMeshType() == MeshType::FILTER)
//       build_filter_mesh(geo, tria);
//     else if(geo.getMeshType() == MeshType::BOX_MESH)
//       build_hyper_rectangle(geo, tria);
//     else
//       throw std::runtime_error("Desired mesh type not valid or not implemented.");
//   }
//   else if(dim ==3)
//   {
//     if(geo.getMeshType() == MeshType::FILE_MESH)
//       build_mesh_from_file(geo, tria);
//     else if(geo.getMeshType() == MeshType::BOX_MESH)
//       simpleRectangle(tria,geo,true);
//     else
//       throw std::runtime_error("Desired mesh type not valid or not implemented.");
//   }
//   else
//     throw std::runtime_error("Mesh generation not implemented for desired dimension.");

//   if(geo.getMeshType() != MeshType::FILE_MESH &&
//       geo.getMeshType() != MeshType::MIXER)
//     set_boundary_and_manifolds(geo, tria); /// this then should not be called in other methods...
// }

template<>
void build_mesh_from_geometry<2>(const Geometry<2>& geo,
                              Triangulation<2>& tria)
{
  if(geo.getMeshType() == MeshType::FILE_MESH)
  {
    build_mesh_from_file(geo, tria);
    attach_mesh_manifolds(geo.getSpheres(), tria); //  ids already set in gmsh
  }
  // else if(geo.getMeshType() == MeshType::BOX_MESH)
  //   simpleRectangle(triangulation,geo,true);
  // else if(geo.getMeshType() == MeshType::SQUARE_CHEESE)
  //   std::cout << "need to implement" << std::endl;
  // else if(geo.getMeshType() == MeshType::HEX_CHEESE)
  //   swissCheese(geo, triangulation);
  else if(geo.getMeshType() == MeshType::MIXER)
    build_mixer_mesh(geo, tria);
  else if(geo.getMeshType() == MeshType::FILTER)
    build_filter_mesh(geo, tria);
  else if(geo.getMeshType() == MeshType::BOX_MESH)
    build_hyper_rectangle(geo, tria);
  else if(geo.getMeshType() == MeshType::SPLITTER)
    build_splitter_mesh(geo, tria);
  else
    throw std::runtime_error("Desired mesh type not valid or not implemented.");

  if(geo.getMeshType() != MeshType::FILE_MESH &&
      geo.getMeshType() != MeshType::MIXER)
    set_boundary_and_manifolds(geo, tria); /// this then should not be called in other methods...
}

template<>
void build_mesh_from_geometry<3>(const Geometry<3>& geo,
                              Triangulation<3>& tria)
{
    if(geo.getMeshType() == MeshType::FILE_MESH)
      build_mesh_from_file(geo, tria);
    else if(geo.getMeshType() == MeshType::BOX_MESH)
      simpleRectangle(tria,geo,true);
    else
      throw std::runtime_error("Desired mesh type not valid or not implemented.");

  if(geo.getMeshType() != MeshType::FILE_MESH &&
      geo.getMeshType() != MeshType::MIXER)
    set_boundary_and_manifolds(geo, tria); /// this then should not be called in other methods...
}

/** HYPER RECTANGLE MESH
*/

template<int dim>
void build_hyper_rectangle(const Geometry<dim>& geo,
                            Triangulation<dim>& tria)
{
  tria.clear();

  std::vector<unsigned int> repetitions;
  std::vector<double> widths;

  for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
    widths.push_back( geo.getWidth(dim_itr) );

  // find min length:
  const double min_width = *std::min_element(widths.begin(),widths.end());
  // std::cout << "min width is " << min_width << std::endl;

  for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
  {
    unsigned int rep = ceil(widths[dim_itr]/min_width);
    // std::cout << "iterations along dim = " << dim_itr << ": " << rep << std::endl;
    repetitions.push_back(rep);
  }

  std::cout << "...generating hyper rectangle" << std::endl;
  dealii::GridGenerator::subdivided_hyper_rectangle(tria,
                                                repetitions,
                                                geo.getBottomLeftPoint(),
                                                geo.getTopRightPoint(),
                                                /*colorize*/ false); // relabeling anyways
}

/* MESH FROM FILE:
* Import mesh from file, add boundary labels and manifolds if needed
*/
template<int dim>
void build_mesh_from_file(const Geometry<dim>& geo,
                          Triangulation<dim>& tria)
{
  dealii::GridIn<dim> grid_in;
  grid_in.attach_triangulation(tria);

  std::cout << "...Using mesh from file: " << geo.getMeshFile() << std::endl;
  std::ifstream input_file(geo.getMeshFile());
  grid_in.read_msh(input_file);
  std::cout << "...mesh read in successfully" << std::endl;
}



/** MIXER FUNCTION IMPLEMENTATIONS
*/

template<>
double build_mixer_mesh<2>(const double left_length,
                      const double right_length,
                      const double height,
                      const double radius,
                      Triangulation<2>& tria,
                      const bool relabel)
{
  if(height < (2.*radius) )
    throw std::invalid_argument("Mixer height must be larger than diameter");

  const double buffer = 0.5*height - radius; // center tile will add to left and right lengths

  const double left_width = left_length - buffer;
  const double right_width = right_length - buffer;

  // std::cout << "buffer = " << buffer << std::endl;
  // std::cout << "left width = " << left_width << " right width = " << right_width << std::endl;

  if(left_width <= 0)
    throw std::invalid_argument("Left width is too small, length already contained in center");

  if(right_width <= 0)
    throw std::invalid_argument("Right width is too small, length already contained in center");

  // need to choose left and right lengths at least of size buffer ...

  Triangulation<2> center;

  construct_mixer_center(height, radius, center);

  add_mixer_ends(left_width, right_width, height, center, tria);

  if(relabel == true)
  {
      const double total_width = left_length + right_length + 2.*radius;
      std::vector<Sphere<2> > spheres;

      // bottom sphere:
      spheres.push_back(Sphere<2>( Point<2>( left_length + radius, 0), radius));
      // top sphere:
      spheres.push_back(Sphere<2>( Point<2>(  left_length + radius, height ), radius));

      // set_boundary_ids(Point<dim>(0,0), Point<dim>(total_width,height), spheres, tria);
      set_edge_boundary_ids(Point<2>(0,0), Point<2>(total_width,height), tria);
      set_sphere_boundary_ids(spheres, tria);

      attach_mesh_manifolds(spheres, tria);
  }

  return buffer;
}

template<>
double build_mixer_mesh<3>(const double /* left_length */,
                      const double /* right_length */,
                      const double /* height */,
                      const double /* radius */,
                      Triangulation<3>& /* tria */,
                      const bool /* relabel */)
{
  throw std::runtime_error("MIXER not yet implemented for 3 dimensions");
}


    /// helpers:
    template<int dim>
    void construct_mixer_center(const double height,
                                const double radius,
                                Triangulation<dim>& center)
    {
      generate_half_circle_hole_tile(height, radius, center); // generated with corner at origin

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


    template<int dim>
    void generate_half_circle_hole_tile(const double height,
                                        const double radius,
                                        Triangulation<dim>& hole_tile)
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
        std::set< typename Triangulation<dim>::active_cell_iterator >
        cells_to_remove;

        // const double mesh_center = 0.5*tile_width;

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

        /// *** testing
        // std::ofstream out("gridcenter.eps");
        // dealii::GridOut grid_out;
        // grid_out.write_eps(hole_tile, out);
        // std::cout << "hole tile written to gridcenter.eps" << std::endl;
    }


    template<int dim>
    void add_mixer_ends(const double left_width,
                        const double right_width,
                        const double height,
                        Triangulation<dim>& center,
                        Triangulation<dim>& tria)
    {
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


    /** geometry wrapper
    */

    template<int dim>
    void build_mixer_mesh(const Geometry<dim>& geo,
                          Triangulation<dim>& tria)
    {
      if(geo.getNumberSpheres() != 2)
        throw std::runtime_error("Geometry should have two spheres for mixer mesh");
      if(geo.getSphereAt(0).getCenter()[0] != geo.getSphereAt(1).getCenter()[0])
        throw std::runtime_error("Two spheres should have same x-coordinate for mixer mesh");

      // /// also check y-coordinates of spheres match y coordinates of corner points:
      // if( (geo.getSphereAt(0).getCenter()[1] != 0) || (geo.getSphereAt(0).getCenter()[]) )

      // // or if spheres are equal...

      /// Assuming geometry is suitable constructed for a mixer...
      const double height = geo.getWidth(1); ///geo.getTopRightPoint()[1] - geo.getBottomLeftPoint()[1];

      // Assuming spheres are aligned along x
      const double radius = geo.getSphereAt(0).getRadius();
      const double center_x = geo.getSphereAt(0).getCenter()[0];

      /// bottom left point should be at origin
      if(geo.getBottomLeftPoint()[0] != 0 || geo.getBottomLeftPoint()[1] != 0)
        throw std::runtime_error("Bottom left corner for mixer mesh should be located at the origin");

      const double left_length = (center_x - radius); // - geo.getBottomLeftPoint()[0];

      const double right_length = geo.getTopRightPoint()[0] - center_x - radius;

      if(right_length < 0)
        std::runtime_error("Top right point should be located to the right of sphere for mixer mesh");

      const double buffer = build_mixer_mesh(left_length, right_length, height, radius, tria);

      set_boundary_and_manifolds(geo, tria, 0.5*buffer);
    }


// ---------------------------------------------------------------------
/** SPLITTER:
*/
template<>
void build_splitter_mesh<2>(const double left_length,
                        const double right_length,
                        const double height,
                        const double radius,
                        Triangulation<2>& tria)
{
  Triangulation<2> center;

  const double buffer = construct_splitter_center(height, radius, center);
  // with bottom left corner of center at origin

  const double left = left_length - buffer;
  const double right = right_length - buffer;

  add_mixer_ends(left, right, height, center, tria);
}

  /// helpers:
  template<int dim>
  double construct_splitter_center(const double height,
                                const double radius,
                                Triangulation<dim>& tria)
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

    // return buffer:
    return outer_radius - inner_radius;
  }

  /** Wrapper for call with geometry object:
  */
  template<int dim>
  void build_splitter_mesh(const Geometry<dim>& geo,
                          Triangulation<dim>& tria)
  {
    if(geo.getNumberSpheres() != 1)
      throw std::runtime_error("Geometry should have one sphere for splitter mesh");

    const double height = geo.getWidth(1); // y -direction
    const double radius = geo.getSphereAt(0).getRadius();
    const double center_x = geo.getSphereAt(0).getCenter()[0];

    const double left_length = center_x - radius;
    const double right_length = geo.getTopRightPoint()[0] - center_x - radius;

    build_splitter_mesh(left_length, right_length, height, radius, tria);
  }


// ---------------------------------------------------------------------
/** FILTER
*/

template<int dim>
void build_filter_mesh(const double left_length,
                      const double center_filter_length,
                      const double right_length,
                      const unsigned int number_channels,
                      const double wall_thickness,
                      const double channel_thickness,
                      Triangulation<dim>& tria,
                      const bool relabel)
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

  attach_fitler_channels(left_length,
                        center_filter_length,
                        number_channels,
                        wall_thickness,
                        channel_thickness,
                        filter_side);

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
  shift_vector[0] = left_length + center_filter_length;
  shift_vector[1] = 0;
  dealii::GridTools::shift(shift_vector, filter_side);

  // merge:
  dealii::GridGenerator::merge_triangulations(filter_side, tria, tria);

  /// @todo option to extrude if 3d

  // label boundaries:
  if(relabel == true)
  {
      const double total_width = left_length + center_filter_length + right_length;
      set_edge_boundary_ids (Point<dim>(0,0), Point<dim>(total_width, height) , tria);

      std::vector<HyperRectangle<dim> > rectangles; /// can also pass in through geometry...

      double bottom_y = channel_thickness;
      double top_y = channel_thickness + wall_thickness;
      for(unsigned int i = 0; i < number_channels-1; i++)
      {
        Point<dim> lower(left_length, bottom_y);
        Point<dim> upper(left_length + center_filter_length, top_y);
        rectangles.push_back(HyperRectangle<dim>(lower, upper));

        bottom_y += (channel_thickness + wall_thickness);
        top_y += (channel_thickness + wall_thickness);

        std::cout << "lower: " << lower << std::endl;
        std::cout << "upper: " << upper << std::endl;
      }

      set_rectangle_boundary_ids(rectangles, tria);

      // attach manifolds ... -- not needed
  }
}

    /// filter helper functions:
    void construct_filter_side(const double width,
                              const double height,
                              const std::vector< std::vector< double > > &  step_sizes,
                              Triangulation<2>& filter_side)
    {
      dealii::GridGenerator::subdivided_hyper_rectangle(filter_side,
                                                        step_sizes,
                                                        Point<2>(0,0),
                                                        Point<2>(width,height));
    }


    void attach_fitler_channels(const double left_length,
                                const double center_filter_length,
                                const unsigned int number_channels,
                                const double wall_thickness,
                                const double channel_thickness,
                                Triangulation<2>& left_side)
    {
      double bottom_y = 0.;
      double top_y = bottom_y + channel_thickness;

      std::vector<unsigned int> repetitions = {1, 1};
      repetitions[0] = std::max( (int)std::floor(center_filter_length/channel_thickness), 1);

      for(unsigned int i = 0; i < number_channels; ++i)
      {
        Triangulation<2> channel;

        dealii::GridGenerator::subdivided_hyper_rectangle(channel,
                                                          repetitions,
                                                          Point<2>(left_length, bottom_y),
                                                          Point<2>(left_length + center_filter_length, top_y));

        dealii::GridGenerator::merge_triangulations(left_side, channel, left_side);

        bottom_y += (channel_thickness+wall_thickness);
        top_y += (channel_thickness+wall_thickness);
      }

    }

    /** geometry wrapper:
    */
    template<int dim>
    void build_filter_mesh(const Geometry<dim>& geo,
                          Triangulation<dim>& tria)
    {
      // assuming geometry class has checked and assembled parameters properly

      const unsigned int number_channels = geo.getNumberRectangles() + 1;
      const double left_length = geo.getRectangleAt(0).getBottomLeft()[0]; // assuming corner at origin
      const double center_filter_length = geo.getRectangleAt(0).getTopRight()[0]
          - geo.getRectangleAt(0).getBottomLeft()[0];
      const double right_length = geo.getTopRightPoint()[0]
          - geo.getRectangleAt(0).getTopRight()[0];
      const double wall_thickness = geo.getRectangleAt(0).getTopRight()[1]
          - geo.getRectangleAt(0).getBottomLeft()[1];

      // thickness = (height - number_rectangles*wall_thickness)/number_channels
      const double channel_thickness = ( geo.getWidth(1)
             - (double)geo.getNumberRectangles()*wall_thickness ) /( (double)number_channels );

      build_filter_mesh(left_length, center_filter_length, right_length,
                        number_channels, wall_thickness, channel_thickness,
                        tria);
    }








/** BOUNDARY LABELS AND MANIFOLD ATTACHMENTS:
*/

  /** Wrapper for call with geometry object:
  */
  template<int dim>
  void set_boundary_and_manifolds(const Geometry<dim>& geo,
                                  Triangulation<dim>& tria,
                                  double sphere_tolerance)
  {
    set_edge_boundary_ids(geo.getBottomLeftPoint(), geo.getTopRightPoint(), tria);
    set_sphere_boundary_ids(geo.getSpheres(), tria, sphere_tolerance);
    set_rectangle_boundary_ids(geo.getRectangles(), tria);
    attach_mesh_manifolds(geo.getSpheres(), tria);
  }

/** GENERAL BOUNDARY AND MANIFOLD SETTING FUNCTIONS
*/

template<int dim>
void set_edge_boundary_ids(const Point<dim>& bottom_left,
                      const Point<dim>& top_right,
                      Triangulation<dim>& tria)
{
  const double edge_tolerance = 1e-8;

  for (typename Triangulation<dim>::active_cell_iterator
         cell = tria.begin_active();
         cell != tria.end();
         ++cell)
       for (unsigned int f=0; f<dealii::GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
                if ( std::fabs( cell->face(f)->center()[dim_itr]
                    - bottom_left[dim_itr] ) < edge_tolerance )
                  cell->face(f)->set_boundary_id(lower_ids[dim_itr]);
                else if ( std::fabs( cell->face(f)->center()[dim_itr]
                    - top_right[dim_itr] ) < edge_tolerance )
                  cell->face(f)->set_boundary_id(upper_ids[dim_itr]);
} // set_edge_boundary_ids()





template<int dim>
void set_sphere_boundary_ids(const std::vector<Sphere<dim> >& spheres,
                            Triangulation<dim>& tria,
                            double sphere_tolerance)
{
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
      unsigned int sphere_bid = id_sphere_begin; // start

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
void set_rectangle_boundary_ids(const std::vector<HyperRectangle<dim> >& rects,
                               Triangulation<dim>& tria)
{
  const double edge_tolerance = 1e-8;

  for (typename Triangulation<dim>::active_cell_iterator
       cell = tria.begin_active();
       cell != tria.end();
       ++cell)
  {
    for(unsigned int f=0; f < dealii::GeometryInfo<dim>::faces_per_cell; ++f)
    {
      unsigned int bid_rect = id_rectangle_begin;

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


}



template<int dim>
void attach_mesh_manifolds(const std::vector<Sphere<dim> >& spheres,
                            Triangulation<dim>& tria)
{
  tria.reset_all_manifolds();

  const unsigned int sphere_id = 10;

  unsigned int man_id = sphere_id;

  for(unsigned int i = 0; i < spheres.size(); i++)
  {
    std::cout << "setting manifold: " << man_id << std::endl;
    tria.set_all_manifold_ids_on_boundary (man_id, man_id);
    dealii::SphericalManifold<dim> sphere_manifold( spheres[i].getCenter() );
    tria.set_manifold(man_id, sphere_manifold);
    ++man_id;
  }
}


/**
* Refinement methods:
*/
template<int dim>
void refine_obstacles(unsigned int obstacle_refinement,
            const Geometry<dim>& geometry,
            Triangulation<dim>& triangulation)
{
  for(unsigned int step = 0; step < obstacle_refinement; ++step)
  {
    for(auto cell : triangulation.active_cell_iterators())
    {
      for(unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
      {

        // for each circle:
      unsigned int number_spheres = geometry.getNumberSpheres();
      for(unsigned int i = 0; i < number_spheres; ++i)
      {
        const Point<dim> center = geometry.getSphereAt(i).getCenter();
        const double radius = geometry.getSphereAt(i).getRadius();

        const double distance_from_center = center.distance(cell->vertex(v));
              if (std::fabs(distance_from_center - radius) < 1e-10)
              {
                  cell->set_refine_flag();
                  // break;
              } // if vertex on circle boundary
          } // for each sphere

          // for each rectangle
          unsigned int number_rectangles = geometry.getNumberRectangles();
          for(unsigned int i=0; i < number_rectangles; ++i)
          {
            if(geometry.getRectangleAt(i).distance_from_border(cell->vertex(v)) < 1e-8)
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


template<int dim>
void refine_boundary(unsigned int boundary_refinement,
            const Geometry<dim>& geometry,
            Triangulation<dim>& triangulation)
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
          const double lower_face = geometry.getBottomLeftPoint()[dim_itr];
          if( std::fabs(cell->vertex(v)[dim_itr] - lower_face) < 1e-8)
            cell->set_refine_flag();

          // upper edge/face:
          const double upper_face = geometry.getTopRightPoint()[dim_itr];
          if( std::fabs(cell->vertex(v)[dim_itr] - upper_face) < 1e-8)
            cell->set_refine_flag();
        } // for each dimension
      } // for each vertex
    } // for each cell in mesh
    triangulation.execute_coarsening_and_refinement();
  } // for each refinement step

  std::cout << "...refined boundary " << boundary_refinement << " times" << std::endl;
}
















/** *******************************************************************************************************
* LEGACY FUNCTIONS:
*/




// template<int dim>
// MyGridGenerator()
// {}

template<int dim>
void
generateGridWithHoles(const Geometry<dim>& geo,
  Triangulation<dim>& triangulation)
{
    simpleRectangle(triangulation,geo,false);

    triangulation.refine_global(3);

    setEdgeBoundaries(geo,triangulation);

    attach_manifolds(geo,triangulation);
}


/*
  generateGrid first querries geometry to see which function to call,
  from there, it passes geometry to the right helper function to generate
  the mesh
*/
template<int dim>
void generateGrid(const Geometry<dim>& geo,
  Triangulation<dim>& triangulation, unsigned int sphere_refinement)
{
	if(dim == 2)
	{
		if(geo.getMeshType() == MeshType::FILE_MESH)
			importMeshFromFile(geo, triangulation);
		else if(geo.getMeshType() == MeshType::BOX_MESH)
			simpleRectangle(triangulation,geo,true);
		else if(geo.getMeshType() == MeshType::SQUARE_CHEESE)
			std::cout << "need to implement" << std::endl;
		else if(geo.getMeshType() == MeshType::HEX_CHEESE)
			swissCheese(geo, triangulation);
    else if(geo.getMeshType() == MeshType::MIXER)
      buildMixerMesh(geo,triangulation,sphere_refinement);
		else
			throw std::runtime_error("Desired mesh type not valid or not implemented.");
	}
	else if(dim ==3)
	{
		if(geo.getMeshType() == MeshType::FILE_MESH)
			importMeshFromFile(geo, triangulation);
		else if(geo.getMeshType() == MeshType::BOX_MESH)
			simpleRectangle(triangulation,geo,true);
		else
			throw std::runtime_error("Desired mesh type not valid or not implemented.");
	}
	else
		throw std::runtime_error("Mesh generation not implemented for desired dimension.");
}


template<int dim>
void
cube_with_spherical_hole(Triangulation<dim>& triangulation)
{
  // put corner at origin, unit radius sphere in 2x2 cube...

  const double inner_radius = 0.25;
  const double outer_radius = (dim == 2 ?
                              0.5*std::sqrt(2)
                              : 0.5*std::sqrt(3));

  const double c = 0.5*std::sqrt(dim);

  Point<dim> center = (dim == 2 ?
                        Point<dim>(c,c) :
                        Point<dim>(0,0,0) );

  const unsigned int n_cells = (dim == 2 ? 4 : 0);

  dealii:: GridGenerator::hyper_shell ( triangulation,
                                        center,
                                        inner_radius,
                                        outer_radius,
              /* const unsigned int  n_cells = */ n_cells,
                        /* bool  colorize = */ true); // outside = 1, inside = 0

  // shift/rotate to have corner at origin

  if(dim == 2)
  {
    // dealii::GridTools::rotate(0.25*dealii::numbers::PI,
    //                         // 2, // z axis
    //                         triangulation);
    dealii::GridTools::shift(Point<dim>(0.5,-0.5),
                            triangulation);
  }
  else if(dim == 3)
  {
        dealii::GridTools::shift(Point<dim>(0.5,0.5,0.5),
                            triangulation);
  }
  // else if dim == 3...

  triangulation.set_all_manifold_ids_on_boundary (0,0);
  triangulation.set_all_manifold_ids_on_boundary (1,1);
     // set manifold_id to mimic boundary_id


  if(dim == 2)
  {
    triangulation.reset_all_manifolds(); //(1); // reset outside to flat
    dealii::SphericalManifold<dim> sphere_manifold(Point<dim>(0.5,0.5));
    triangulation.set_manifold(0, sphere_manifold);
  }
  else if(dim == 3)
  {
    triangulation.reset_all_manifolds(); //(1); // reset outside to flat
    dealii::SphericalManifold<dim> sphere_manifold(Point<dim>(0.5,0.5,0.5));
    triangulation.set_manifold(0, sphere_manifold);

	// triangulation.refine_global(1);
  }
  // triangulation.refine_global(1);

  // transfinite manifold:
  // dealii::TransfiniteInterpolationManifold<dim> transfinite_manifold;
  // transfinite_manifold.initialize(triangulation);
  // triangulation.set_manifold(1, transfinite_manifold);

  // display vertices:
  // for( const auto cell : triangulation.active_cell_iterators() )
  // {
  //   std::cout << "cell center: " << cell->center() << std::endl;
  //   for(unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
  //     std::cout << "cell vertex: " << cell->vertex(v) << std::endl;
  // }


} // cube with hole()


template<int dim>
void square_pipe_with_spherical_hole(Triangulation<dim>& triangulation, double scale_factor)
{
	Triangulation<dim> auxillary;

	dealii::GridGenerator::hyper_cube(auxillary); // 0x1 by default

	cube_with_spherical_hole(triangulation);

	Tensor<1,dim> shift_vector;
	shift_vector[0] = 1.0;
	for(unsigned int dim_itr = 1; dim_itr < dim; dim_itr++)
		shift_vector[dim_itr] = 0.0;

	dealii::GridTools::shift(shift_vector, triangulation); // shift right in x

	dealii::GridGenerator::merge_triangulations(auxillary,triangulation,triangulation);

	shift_vector[0] = 2.0;
	dealii::GridTools::shift(shift_vector, auxillary);
	dealii::GridGenerator::merge_triangulations(triangulation,auxillary,triangulation);

	dealii::GridTools::scale(scale_factor, triangulation);

	// need to reattach manifolds... and assign boundary ids ...
	set_pipe_edges(triangulation, scale_factor);
	set_pipe_manifolds(triangulation, scale_factor);
}

template<int dim>
void
set_pipe_edges(Triangulation<dim>& triangulation, double scale_factor)
{

	// edges at x = 0, x = 3
	// y = 0, y = 1; z = 0, z = 1

	const double sphere_tolerance = 0.2; // may need to make this better...
	// can maybe make coarser, since coarse mesh shouldn't have a problem

	const double edge_tolerance = 1e-8;

	const unsigned int id_left = 0;
	const unsigned int id_right = 1;
	const unsigned int id_top = 2;
	const unsigned int id_bottom = 3;

	// 3d: ( z direction )
	const unsigned int id_front = 4;
	const unsigned int id_back = 5;

	const unsigned int id_other = 7;

	// triangulation.set_all_boundary_ids(id_other);

	const unsigned int id_sphere = 10;

	const double radius = 0.25;
	const Point<dim> sphere_center = (dim == 2 ?
                        Point<dim>(scale_factor*1.5, scale_factor*0.5) :
                        Point<dim>(scale_factor*1.5, scale_factor*0.5, scale_factor*0.5) );

	std::array< unsigned int , 3 > lower_ids = {id_left, id_bottom, id_back};
	std::array< unsigned int , 3 > upper_ids = {id_right, id_top, id_front};

	Point<dim> lower;
	Point<dim> upper;
	upper[0] = scale_factor*3.;

	for(unsigned int dim_itr = 1; dim_itr < dim; ++dim_itr)
		upper[dim_itr] = scale_factor*1.;


	for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
	     cell != triangulation.end();
	     ++cell)
	{
 		 for (unsigned int f=0; f<dealii::GeometryInfo<dim>::faces_per_cell; ++f)
 		 {
		    if (cell->face(f)->at_boundary())
		    {
		    	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		    	{
			        if ( std::fabs( cell->face(f)->center()[dim_itr]
	  		          - lower[dim_itr] ) < edge_tolerance )
			        {
				          cell->face(f)->set_boundary_id(lower_ids[dim_itr]);
				          std::cout << lower_ids[dim_itr] << ": set at " <<
				          	cell->face(f)->center() << std::endl;
			        }
			        else if ( std::fabs( cell->face(f)->center()[dim_itr]
			          - upper[dim_itr] ) < edge_tolerance )
			        {
				          cell->face(f)->set_boundary_id(upper_ids[dim_itr]);
  				          std::cout << upper_ids[dim_itr] << ": set at " <<
				          	cell->face(f)->center() << std::endl;
			        }
					// else
					// {
					// 	cell->face(f)->set_boundary_id(id_other);
			  //         std::cout << id_other << ": set at " <<
				 //          	cell->face(f)->center() << std::endl;
					// }
		    	} // for dimensions
		    } // if at boundary
		}

		// sphere in center:
		for(unsigned int f=0; f < dealii::GeometryInfo<dim>::faces_per_cell; ++f)
		{
			const double distance_from_center = sphere_center.distance(cell->face(f)->center());

			if( std::fabs(distance_from_center - radius) < sphere_tolerance )
			{
				std::cout << "in sphere: " << cell->face(f)->center() << std::endl;
				cell->face(f)->set_boundary_id(id_sphere);
			}
		} // for faces

    } // for cells

}


template<int dim>
void
set_pipe_manifolds(Triangulation<dim>& triangulation, double scale_factor)
{
  triangulation.set_all_manifold_ids_on_boundary (10,10);

  if(dim == 2)
  {
    triangulation.reset_all_manifolds(); //(1); // reset outside to flat
    dealii::SphericalManifold<dim> sphere_manifold(Point<dim>(scale_factor*1.5,scale_factor*0.5));
    triangulation.set_manifold(10, sphere_manifold);
  }
  else if(dim == 3)
  {
    triangulation.reset_all_manifolds(); //(1); // reset outside to flat
    dealii::SphericalManifold<dim> sphere_manifold(Point<dim>(scale_factor*1.5,scale_factor*0.5,scale_factor*0.5));
    triangulation.set_manifold(10, sphere_manifold);
  }
}

template<int dim>
void
generateRectangleWithHole(double height, Triangulation<2,2>& triangulation)
{
  if(dim != 2)
    throw std::runtime_error("not implemented");

  Triangulation<2,2> temporary;

  // start with first third:
  Point<2> bottomLeft(0,0);
  Point<2> topCorner(height,height);

  std::vector<unsigned int> subdivisions;
  subdivisions.push_back(2);
  subdivisions.push_back(2);

  dealii::GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            subdivisions,
                                            bottomLeft,
                                            topCorner);

  // center third:
  const double outer_radius = 0.5*height;
  const double inner_radius = 0.5*outer_radius;

  dealii::GridGenerator::hyper_cube_with_cylindrical_hole(temporary,
                                                  inner_radius,
                                                  outer_radius
                                                  );

  Tensor<1,2> shift_vector;
  shift_vector[0] = outer_radius + height;
  shift_vector[1] = outer_radius;

  dealii::GridTools::shift(shift_vector,temporary);

  dealii::GridGenerator::merge_triangulations(temporary,triangulation,triangulation);

  // final third:
  temporary.clear();
  dealii::GridGenerator::subdivided_hyper_rectangle(temporary,
                                                    subdivisions,
                                                    bottomLeft,
                                                    topCorner);
  shift_vector[0] = 2.*height;
  shift_vector[1] = 0.;


  dealii::GridTools::shift(shift_vector,temporary);
  dealii::GridGenerator::merge_triangulations(temporary,triangulation,triangulation);

  // assign boundary ids and manifolds:
  // boundary ids:
  const double tolerance = 0.1;
  // can maybe make coarser, since coarse mesh shouldn't have a problem


  const unsigned int id_left = 0;
  const unsigned int id_right = 1;
  const unsigned int id_top = 2;
  const unsigned int id_bottom = 3;
  const unsigned int id_other = 5;

  for (typename Triangulation<2>::active_cell_iterator
       cell = triangulation.begin_active();
     cell != triangulation.end();
     ++cell)
  for (unsigned int f=0; f<dealii::GeometryInfo<2>::faces_per_cell; ++f)
    if (cell->face(f)->at_boundary())
    {
      // std::cout << "yo" << std::endl;
        if (cell->face(f)->center()[0] == 0)
          cell->face(f)->set_boundary_id(id_left);
        else if (cell->face(f)->center()[0] == 3*height)
          cell->face(f)->set_boundary_id(id_right);
        else if (cell->face(f)->center()[1] == height)
          cell->face(f)->set_boundary_id(id_top);
        else if (cell->face(f)->center()[1] == 0)
          cell->face(f)->set_boundary_id(id_bottom);
        else
          cell->face(f)->set_boundary_id(id_other);
    } // if at boundary

  // set sphere boundaryies:
  Point<2> center(1.5*height, 0.5*height);
  const double radius = inner_radius;

  unsigned int sphere_bid = 10; // start at 10

  for(typename Triangulation<2>::active_cell_iterator
        cell = triangulation.begin_active();
        cell != triangulation.end();
        ++cell)
    {
      for(unsigned int f=0; f < dealii::GeometryInfo<2>::faces_per_cell; ++f)
      {
        const double distance_from_center = center.distance(cell->face(f)->center());

        if( std::fabs(distance_from_center - radius) < tolerance )
            cell->face(f)->set_boundary_id(sphere_bid);
      } // for faces
    } // for cells

  // manifolds:
    triangulation.reset_all_manifolds();

    const unsigned int inner_manifold_id = 0;
    triangulation.set_all_manifold_ids(inner_manifold_id);

    // unsigned int sphere_bid = 6; // start at 6 by convention .. may need to change for 3D
    unsigned int sphere_manifold_id = 1; // start at 1

    // for(unsigned int i = 0; i < number_spheres; i++)
    // {
      triangulation.set_all_manifold_ids_on_boundary(sphere_bid, sphere_manifold_id);

      dealii::SphericalManifold<dim> sphere_manifold(center); // polar vs circle?
      triangulation.set_manifold (sphere_manifold_id, sphere_manifold);

      // ++sphere_bid;
      // ++sphere_manifold_id;
    // }


    dealii::TransfiniteInterpolationManifold<dim> inner_manifold;
    inner_manifold.initialize(triangulation);
    triangulation.set_manifold (inner_manifold_id, inner_manifold);
}



template<int dim>
void
importMeshFromFile(const Geometry<dim>& geo,
  Triangulation<dim>& triangulation)
{
  dealii::GridIn<dim> grid_in;
	grid_in.attach_triangulation(triangulation);

  std::cout << "...Using mesh from file: " << geo.getMeshFile() << std::endl;
	std::ifstream input_file(geo.getMeshFile());
	grid_in.read_msh(input_file);
  std::cout << "...mesh read" << std::endl;

  // scale if needed:
  if(geo.getScaleFactor() != 1)
    dealii::GridTools::scale(geo.getScaleFactor(), triangulation);

  // color boundaries...
  // setEdgeBoundaries(geo.getBottomLeftPoint()[0], geo.getTopRightPoint()[0],
  //     geo.getTopRightPoint()[1], geo.getBottomLeftPoint()[1],
  //     triangulation);
  setEdgeBoundaries(geo, triangulation);


  // attach manifolds:
  attach_manifolds(geo,triangulation);

}


template<int dim>
void
attach_manifolds(const Geometry<dim>& geo,
    Triangulation<dim>& triangulation)
{
    triangulation.reset_all_manifolds();

    const unsigned int inner_manifold_id = 0;
    triangulation.set_all_manifold_ids(inner_manifold_id);

    const unsigned int number_spheres = geo.getNumberSpheres();
    unsigned int sphere_bid = 10; // start at 6 by convention .. may need to change for 3D
    unsigned int sphere_manifold_id = 1; // start at 1

    for(unsigned int i = 0; i < number_spheres; i++)
    {
      triangulation.set_all_manifold_ids_on_boundary(sphere_bid, sphere_manifold_id);

      dealii::SphericalManifold<dim> sphere_manifold(geo.getSphereAt(i).getCenter()); // polar vs circle?
      triangulation.set_manifold (sphere_manifold_id, sphere_manifold);

      ++sphere_bid;
      ++sphere_manifold_id;
    }


    dealii::TransfiniteInterpolationManifold<dim> inner_manifold;
    inner_manifold.initialize(triangulation);
    triangulation.set_manifold (inner_manifold_id, inner_manifold);
}


template<int dim>
void
simpleRectangle(Triangulation<dim> &tria,
			const Geometry<dim>& geo, const bool colorize)
{
  // @todo ***LEFT OFF HERE ***
  // subdivide to account for non equal widths ...
  tria.clear();

  std::vector<unsigned int> repetitions;
  std::vector<double> widths;

  for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
    widths.push_back( geo.getWidth(dim_itr) );

  // find min length:
  const double min_width = *std::min_element(widths.begin(),widths.end());
  std::cout << "min width is " << min_width << std::endl;

  for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
  {
    unsigned int rep = ceil(widths[dim_itr]/min_width);
    std::cout << "iterations along dim = " << dim_itr << ": " << rep << std::endl;
    repetitions.push_back(rep);
  }

std::cout << "generating hyper rectangle" << std::endl;
	dealii::GridGenerator::subdivided_hyper_rectangle(tria,
                                                  repetitions,
                                                  geo.getBottomLeftPoint(),
	                                 								geo.getTopRightPoint(),
                                                  colorize);

  std::cout << "rectangle set" << std::endl;
}














template<int dim>
void grid_hole_cutouts(const Geometry<dim>& geo,
                      Triangulation<dim>& tria)
{
  simpleRectangle(tria,geo); // @todo, clean up grid generation, make all functions with same order of parameters

  tria.refine_global(2);
  refine_along_spheres(geo, tria);

  // cutout_spheres(geo, tria); // doesn't work for adaptively refined meshes...
}


template<int dim>
void refine_along_spheres(const Geometry<dim> & geometry,
                        Triangulation<dim>& triangulation)
{

  const unsigned int sphere_refinement = 4; // make adaptive...

  const unsigned int number_spheres = geometry.getNumberSpheres();


    for(unsigned int step = 0; step < sphere_refinement; ++step)
    {
        const double min_diameter = 1.0 *
          dealii::GridTools::minimal_cell_diameter(triangulation); // use as a guide...
        // refine radius +/- cell diameter,
        // until local cell diameter is ...

      for(unsigned int i = 0; i < number_spheres; i++)
      {
        const Point<dim> center = geometry.getSphereAt(i).getCenter();
        const double radius = geometry.getSphereAt(i).getRadius();

          for(auto cell : triangulation.active_cell_iterators())
          {
            for(unsigned int v = 0; v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
            {
              const double distance_from_center = center.distance(cell->vertex(v));
              const bool within_gap = ( distance_from_center <= (radius+min_diameter) )
               && ( distance_from_center > (radius - 0.2*min_diameter) );
                  if (within_gap)
                  {
                    cell->set_refine_flag();
                    // std::cout << "is within gap? ... " << within_gap << std::endl;
                    break;
                  } // if vertex on circle boundary
            } // for each vertex
          } // for each cell in mesh

    } // for each sphere

   triangulation.execute_coarsening_and_refinement();
  } // for each refinement step

}


template<int dim>
void cutout_spheres(const Geometry<dim> & geometry,
                        Triangulation<dim>& triangulation)
{
  const unsigned int number_spheres = geometry.getNumberSpheres();

  // find cells to remove:
  std::set< typename Triangulation<dim>::active_cell_iterator >
  cells_to_remove;

  for(unsigned int i = 0; i < number_spheres; i++)
  {
    const Point<dim> center = geometry.getSphereAt(i).getCenter();
    const double radius = geometry.getSphereAt(i).getRadius();

    for (const auto cell : triangulation.active_cell_iterators())
    {
      const double distance_from_center = center.distance(cell->center());
      if ( std::fabs(distance_from_center - radius) < 1e-8)
      {
        // std::cout <<" cell center at: " << cell->center() << std::endl;
        cells_to_remove.insert(cell);
      }
    }
  } // for spheres

  // remove LHS:
  dealii::GridGenerator::create_triangulation_with_removed_cells(triangulation,
                                                                cells_to_remove,
                                                                triangulation);
}




















template<int dim>
void
swissCheese(const Geometry<dim>& geo, Triangulation<dim> &tria)
{
	if(dim != 2)
		throw std::runtime_error("Function not implemented for dim != 2");

  if(dim == 2)
    throw std::runtime_error("Need to implement swiss cheese based on geometry dependence");

  geo.printInfo(std::cout);

  Triangulation<dim>    hole_tile;
  Triangulation<dim>    half_circle;

  const double tile_width = 5.0;
  const double hole_radius = 0.25;
  const double tile_height = 0.5*std::sqrt(3)*tile_width;
  // need to check here if tilable ... if radius is small enough given spacing
  // radius must be smaller than height

  squished_square_with_hole(hole_radius, tile_width, hole_tile);
  squished_half_circle_hole(hole_radius, tile_width, half_circle); // already shifted to origin

  // need two rows:
  const double full_width = 9*tile_width; //geo.getXMax()-geo.getXMin();
  const double full_height = 8*tile_height; //geo.getYMax()-geo.getYMin();

  Triangulation<dim> bottom_row;
  bottom_row.copy_triangulation(hole_tile);
  constructSwissBottomRow(tile_width,full_width,hole_tile,bottom_row);

  Triangulation<dim> second_row;
  second_row.copy_triangulation(half_circle);
  constructSwissSecondRow(tile_width,full_width,
    half_circle,hole_tile,second_row);

  swissRowsToFullMesh(tile_height,full_height,
                      bottom_row,second_row,tria);

  const double left = 0; // geo.getXMin(); // should be 0
  const double bottom = 0; //get.getYMin();
  const double top = full_height;
  const double right = full_width;
  setEdgeBoundaries(left,right,top,bottom,tria);
}


template<int dim>
void
swissRowsToFullMesh(double tile_height, double full_height,
    Triangulation<dim>& bottom_row, Triangulation<dim>& second_row,
    Triangulation<dim>& tria)
{
	if(dim != 2)
		throw std::runtime_error("Function not implemented for dim != 2");

  // const unsigned int n_x_stackings = full_width/base_width - 1;
  const unsigned int n_y_stackings = full_height/tile_height - 2;

  std::cout << n_y_stackings << " y stackings" << std::endl;

  // alternate rows:
  Tensor<1,dim> shift_vector;
  shift_vector[0] = 0;
  shift_vector[1] = tile_height;

  tria.clear();
  tria.copy_triangulation(bottom_row);
  dealii::GridTools::shift(shift_vector,second_row);
  dealii::GridGenerator::merge_triangulations(bottom_row,second_row,tria);

  shift_vector[1] = 2*tile_height;
  for(unsigned int i = 0; i < n_y_stackings; i+=2)
  {
    dealii::GridTools::shift(shift_vector,bottom_row);
    dealii::GridTools::shift(shift_vector,second_row);
    dealii::GridGenerator::merge_triangulations(bottom_row,tria,tria);
    dealii::GridGenerator::merge_triangulations(second_row,tria,tria);
  }
}

template<int dim>
void
constructSwissBottomRow(double tile_width, double full_width,
   Triangulation<dim>& hole_tile, Triangulation<dim>& bottom_row)
{
	if(dim != 2)
		throw std::runtime_error("Function not implemented for dim != 2");

  const int number_tiles = full_width/tile_width - 1;
  Tensor<1,dim> shift_vector;
  shift_vector[0] = tile_width;
  shift_vector[1] = 0;

  for(unsigned int i = 0; i < number_tiles; i++)
  {
    dealii::GridTools::shift(shift_vector,hole_tile);
    dealii::GridGenerator::merge_triangulations(hole_tile,bottom_row,bottom_row);
  }

  // // shift back
  // std::cout << "number tiles: " << number_tiles << std::endl;
  shift_vector[0] = -number_tiles*tile_width;
  shift_vector[1] = 0;
  dealii::GridTools::shift(shift_vector,hole_tile);
}



template<int dim>
void
squished_square_with_hole(double hole_radius, double tile_width,
  Triangulation<dim>& triangulation)
{
	if(dim != 2)
		throw std::runtime_error("Function not implemented for dim != 2");

  Triangulation<2>    hole_tile;

  dealii::GridGenerator::hyper_cube_with_cylindrical_hole(hole_tile,
                                                hole_radius,
                                                0.5*tile_width);
  Tensor<1,dim> shift_vector;
  shift_vector[0] = 0.5*tile_width;
  shift_vector[1] = 0.5*tile_width;
  dealii::GridTools::shift(shift_vector,hole_tile); // put corner at origin

  auto transformFunction = std::bind(squish_transform,std::placeholders::_1, tile_width);
  dealii::GridTools::transform (transformFunction, hole_tile);

  // shift back down to origin:
  double difference = (1.0 - 0.5*std::sqrt(3))*tile_width;
  double down_shift = -0.5*difference;
  shift_vector[0] = down_shift;
  shift_vector[1] = 0;
  dealii::GridTools::shift(shift_vector,hole_tile);

  triangulation.copy_triangulation(hole_tile);
}


template<int dim>
void
constructSwissSecondRow(double tile_width, double full_width,
    Triangulation<dim>& half_circle, Triangulation<dim>&hole_tile,
    Triangulation<dim>& second_row)
{
	if(dim != 2)
		throw std::runtime_error("Function not implemented for dim != 2");

  const int number_tiles = full_width/tile_width - 1;
  Tensor<1,dim> shift_vector;
  shift_vector[0] = 0.5*tile_width;
  shift_vector[1] = 0;
  dealii::GridTools::shift(shift_vector,hole_tile);

  shift_vector[0] = tile_width;
  for(unsigned int i = 0; i < number_tiles; i++)
  {
    dealii::GridGenerator::merge_triangulations(hole_tile,second_row,second_row);
    dealii::GridTools::shift(shift_vector,hole_tile);
  }

  // rotate, shift, and attach half circle to other end:
  const double angle = dealii::numbers::PI;
  const double tile_height = 0.5*std::sqrt(3)*tile_width;

  dealii::GridTools::rotate(angle,half_circle);
  shift_vector[0] = tile_height + (number_tiles)*tile_width;
  shift_vector[1] = tile_width;
  dealii::GridTools::shift(shift_vector,half_circle);

  dealii::GridGenerator::merge_triangulations(half_circle,second_row,second_row);
}


template<int dim>
void
generateHalfCircleHole(double hole_radius, double tile_width,
  Triangulation<dim>& triangulation, unsigned int sphere_refinement)
{
	if(dim != 2)
		throw std::runtime_error("Function not implemented for dim != 2");

  // start with rectangle with hole:
  Triangulation<dim>    hole_tile;
  // const double hole_radius = 0.25;
  // const double tile_width = 1.0;

  dealii::GridGenerator::hyper_cube_with_cylindrical_hole(hole_tile,
                                                hole_radius,
                                                0.5*tile_width);

  // refine circle:
  Point<2> center(0,0);
  const double radius = hole_radius;

  for(unsigned int step = 0; step < sphere_refinement; ++step)
  {
    std::cout << "sphere refinement step: " << step << std::endl;
    for(auto cell : hole_tile.active_cell_iterators())
    {
      for(unsigned int v = 0; v < dealii::GeometryInfo<2>::vertices_per_cell; ++v)
      {
        const double distance_from_center = center.distance(cell->vertex(v));
        if (std::fabs(distance_from_center - radius) < 1e-10)
        {
          cell->set_refine_flag();
          break;
        } // if vertex on circle boundary
      } // for each vertex
    } // for each cell in mesh
    hole_tile.execute_coarsening_and_refinement();
  } // for each refinement step


  // shift so corner at origin:
  Tensor<1,dim> shift_vector;
  shift_vector[0] = 0.5*tile_width;
  shift_vector[1] = 0.5*tile_width;
  dealii::GridTools::shift(shift_vector,hole_tile);


  // find cells to remove:
  std::set< typename Triangulation<dim>::active_cell_iterator >
  cells_to_remove;

  const double mesh_center = 0.5*tile_width;

  for (const auto cell : hole_tile.active_cell_iterators())
  {
      if (cell->center()[0] < mesh_center)
      {
        // std::cout <<" cell center at: " << cell->center() << std::endl;
        cells_to_remove.insert(cell);
      }
  }

  // remove LHS:
  dealii::GridGenerator::create_triangulation_with_removed_cells(hole_tile,
                                                          cells_to_remove,
                                                          triangulation);

  // shift new corner to origin
  shift_vector[0] = -0.5*tile_width;
  shift_vector[1] = 0;
  dealii::GridTools::shift(shift_vector,triangulation);
}

template<int dim>
void
squished_half_circle_hole(double hole_radius, double tile_width,
  Triangulation<dim>& triangulation)
{
	if(dim != 2)
		throw std::runtime_error("Function not implemented for dim != 2");

  generateHalfCircleHole(hole_radius, tile_width,triangulation);

  auto transformFunction = std::bind(squish_transform,std::placeholders::_1, tile_width);
  dealii::GridTools::transform (transformFunction, triangulation);

  // shift back down to origin:
  double difference = (1.0 - 0.5*std::sqrt(3))*tile_width;
  double down_shift = -0.5*difference;

  Tensor<1,dim> shift_vector;
  shift_vector[0] = down_shift;
  shift_vector[1] = 0;

  dealii::GridTools::shift(shift_vector,triangulation);
}


template<int dim>
void
buildMixerMesh(const Geometry<dim>& geo,
  Triangulation<dim>& triangulation, unsigned int sphere_refinement)
{
  // std::cout << "In buildMixerMesh() ... can perhaps extrude for 3d ..."
  //   << std::endl;

  if(geo.getNumberSpheres() < 1)
    throw std::invalid_argument("Need spheres for mixer mesh.");

  const double hole_radius = geo.getSphereAt(0).getRadius();
  const double tile_width = geo.getWidth(1); // along y
  const double full_width = geo.getWidth(0); // along x

  Triangulation<dim> hole_tile;

  generateHalfCircleHole(hole_radius, tile_width, hole_tile, sphere_refinement);
  generateHalfCircleHole(hole_radius, tile_width, triangulation, sphere_refinement);
  constructMixerCenter(tile_width,hole_tile,triangulation);
  addMixerEnds(tile_width,full_width,triangulation, 0); //sphere_refinement);

  dealii::GridTools::remove_hanging_nodes(triangulation);

  setEdgeBoundaries(geo, triangulation);

  querryBoundaryFaces(triangulation); // check

  attach_manifolds(geo, triangulation);

}

template<int dim>
void
constructMixerCenter(double tile_width,
  Triangulation<dim>& hole_tile,
  Triangulation<dim>& triangulation)
{
  const double angle = 0.5*dealii::numbers::PI;
  dealii::GridTools::rotate(angle,triangulation); // bottom
  dealii::GridTools::rotate(-angle,hole_tile);  // top

  Tensor<1,dim> shift_vector;
  shift_vector[0] = tile_width;
  shift_vector[1] = 0;
  // move bottom left corner back to origin
  dealii::GridTools::shift(shift_vector,triangulation);

  shift_vector[0] = 0;
  shift_vector[1] = tile_width;
  // move above bottom tile
  dealii::GridTools::shift(shift_vector,hole_tile);

  // merge:
  dealii::GridGenerator::merge_triangulations(hole_tile,triangulation,
                                              triangulation);
}

template<int dim>
void
addMixerEnds(double tile_width, double full_width,
  Triangulation<dim>& triangulation,
  unsigned int sphere_refinement)
{
  const double center_shift = 0.5*(full_width-tile_width);

  Tensor<1,dim> shift_vector;
  shift_vector[0] = center_shift;
  shift_vector[1] = 0;
  dealii::GridTools::shift(shift_vector,triangulation);

  std::vector<unsigned int> repetitions{2 + sphere_refinement, 2 + sphere_refinement};

  Point<dim> bottomLeft(0,0);
  Point<dim> topRight(center_shift,tile_width);

  Triangulation<dim> rectangle;
  dealii::GridGenerator::subdivided_hyper_rectangle(
                                          rectangle,
                                          repetitions,
                                          bottomLeft,
                                          topRight);

  dealii::GridGenerator::merge_triangulations(rectangle,triangulation,
                                            triangulation);

  // add right side
  shift_vector[0] = center_shift + tile_width;
  dealii::GridTools::shift(shift_vector,rectangle);
  dealii::GridGenerator::merge_triangulations(rectangle,triangulation,
                                            triangulation);

}


/*


void MyGridGenerator::boxCheese(Triangulation<2,2> &tria, const Geometry& geo)
{
  Triangulation<2>    square_tile;
  Triangulation<2>    hole_tile;

  const double tile_width = geo.getTileWidth(); // this should eventually be a member of geometry
  const double hole_radius = geo.circleAt(0)->getRadius();

  const Point<2> bottomLeft(0, 0);
  const Point<2> topRight(tile_width,tile_width);

  std::vector< unsigned int > repetitions;
  repetitions.push_back(2);
  repetitions.push_back(2); // 2x2 cells...

  GridGenerator::subdivided_hyper_rectangle(square_tile,
                                            repetitions,
                                            bottomLeft,
                                            topRight);

  GridGenerator::hyper_cube_with_cylindrical_hole(hole_tile,
                                                  hole_radius,
                                                  0.5*tile_width);

  // shift so corner at origin:
  Tensor<1,2> shift_vector;
  shift_vector[0] = 0.5*tile_width;
  shift_vector[1] = 0.5*tile_width;
  dealii::GridTools::shift(shift_vector,hole_tile);

  // construct base tile with hole on bottom left and top right:
  constructBaseTile(tile_width,hole_tile,square_tile,tria);


  const double full_width = geo.getXMax()-geo.getXMin();
  const double full_height = geo.getYMax()-geo.getYMin();
  const double base_width = 2*tile_width;

  baseTileToFullMesh(full_width,full_height,base_width,tria);

  // get these from geometry:
  const double left = 0; // geo.getXMin(); // should be 0
  const double bottom = 0; //get.getYMin();
  const double top = full_height;
  const double right = full_width;
  setEdgeBoundaries(left,right,top,bottom,tria);

}


void MyGridGenerator::constructBaseTile(double tile_width, Triangulation<2,2>& hole_tile,
   Triangulation<2,2>& square_tile, Triangulation<2,2>& tria)
{
  Tensor<1,2> shift_vector;
  shift_vector[0] = tile_width;
  shift_vector[1] = 0.0;

  // right of hole:
  dealii::GridTools::shift(shift_vector,square_tile);
  GridGenerator::merge_triangulations(hole_tile,square_tile,tria);

  // top left corner:
  shift_vector[0] = -tile_width;
  shift_vector[1] = tile_width;
  dealii::GridTools::shift(shift_vector,square_tile);
  GridGenerator::merge_triangulations(square_tile,tria,tria);

  // top right corner:
  shift_vector[0] = tile_width;
  shift_vector[1] = tile_width;
  dealii::GridTools::shift(shift_vector,hole_tile);
  GridGenerator::merge_triangulations(hole_tile,tria,tria);
}


void MyGridGenerator::baseTileToFullMesh(double full_width, double full_height,
    double base_width, Triangulation<2,2>& tria)
{
  const unsigned int n_x_stackings = full_width/base_width - 1;
  const unsigned int n_y_stackings = full_height/base_width - 1;

  // std::cout << n_x_stackings << " x stackings" << std::endl;
  // std::cout << n_y_stackings << " y stackings" << std::endl;

  Triangulation<2,2> base_tile;
  base_tile.copy_triangulation(tria);

  Tensor<1,2> shift_vector;

  // stack along width:
  shift_vector[0] = base_width;
  shift_vector[1] = 0;

  for(unsigned int i = 0; i < n_x_stackings; i++)
  {
    dealii::GridTools::shift(shift_vector,base_tile);
    GridGenerator::merge_triangulations(base_tile,tria,tria);
  }

  base_tile.clear();
  base_tile.copy_triangulation(tria);

  // stack along height:
  shift_vector[0] = 0;
  shift_vector[1] = base_width;

  for(unsigned int i = 0; i < n_y_stackings; i++)
  {
    dealii::GridTools::shift(shift_vector,base_tile);
    GridGenerator::merge_triangulations(base_tile,tria,tria);
  }

}

*/

// print out locations of boundary faces and their boundary id:
template<int dim>
void
querryBoundaryFaces(const Triangulation<dim>& triangulation)
{
  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
     cell != triangulation.end();
     ++cell)
  for (unsigned int f=0; f < dealii::GeometryInfo<dim>::faces_per_cell; ++f)
    if (cell->face(f)->at_boundary())
    {
      // const unsigned int b_id = cell->face(f)->boundary_indicator();
      std::cout << "boundary cell center: " << cell->face(f)->center() << std::endl;
      // std::cout << "\t boundary id: " << b_id << std::endl;
    }


  // std::vector<types::boundary_id> boundary_ids = triangulation.get_boundary_indicators();

  // for(unsigned int i = 0; i < boundary_ids.size(); i++)
  //   std::cout << "\t boundary id: " << boundary_ids[i] << std::endl;
  // if (cell->face(f)->center()[0] == -1)
  //   cell->face(f)->set_boundary_indicator (42);
}


/*

void MyGridGenerator::generateMergedSquare(const Point<2> bottomLeft,
  const Point<2> topRight, Triangulation<2,2>& triangulation)
{
  Triangulation<2,2> temp_mesh;

  std::vector< unsigned int > repetitions;
  repetitions.push_back(2);
  repetitions.push_back(2); // 2x2 cells...

  GridGenerator::subdivided_hyper_rectangle(temp_mesh,
                                            repetitions,
                                            bottomLeft,
                                            topRight);

  triangulation.copy_triangulation(temp_mesh);

  Tensor<1,2> shift_vector;
  const double width = topRight(0)-bottomLeft(0);
  shift_vector[0] = width;
  shift_vector[1] = 0;

  dealii::GridTools::shift(shift_vector,temp_mesh);

  GridGenerator::merge_triangulations(temp_mesh,triangulation,triangulation);

  // add a hole at end:
  // Triangulation<2,2> hole_mesh;

  // const double hole_radius = 0.25*width;
  // GridGenerator::hyper_cube_with_cylindrical_hole(hole_mesh,
  //                                                 hole_radius,
  //                                                 0.5*width);

  // shift_vector[0] = 0.5*width + 2*width;
  // shift_vector[1] = 0.5*width;
  // dealii::GridTools::shift(shift_vector,hole_mesh);

  // GridGenerator::merge_triangulations(hole_mesh,triangulation,triangulation);

  const double left = bottomLeft(0);
  const double bottom = bottomLeft(1);
  const double top = topRight(1);
  const double right = bottomLeft(0) + 2*width;
  setEdgeBoundaries(left,right,top,bottom,triangulation);
}



void MyGridGenerator::generateMergedSquarePlusHole(const Point<2> bottomLeft,
  const Point<2> topRight, Triangulation<2,2>& triangulation)
{
  Triangulation<2,2> temp_mesh;

  std::vector< unsigned int > repetitions;
  repetitions.push_back(2);
  repetitions.push_back(2); // 2x2 cells...

  GridGenerator::subdivided_hyper_rectangle(temp_mesh,
                                            repetitions,
                                            bottomLeft,
                                            topRight);

  triangulation.copy_triangulation(temp_mesh);

  Tensor<1,2> shift_vector;
  const double width = topRight(0)-bottomLeft(0);
  shift_vector[0] = width;
  shift_vector[1] = 0;

  dealii::GridTools::shift(shift_vector,temp_mesh);

  GridGenerator::merge_triangulations(temp_mesh,triangulation,triangulation);

  // add a hole at end:
  Triangulation<2,2> hole_mesh;

  const double hole_radius = 0.25*width;
  GridGenerator::hyper_cube_with_cylindrical_hole(hole_mesh,
                                                  hole_radius,
                                                  0.5*width);

  shift_vector[0] = 0.5*width + 2*width;
  shift_vector[1] = 0.5*width;
  dealii::GridTools::shift(shift_vector,hole_mesh);

  GridGenerator::merge_triangulations(hole_mesh,triangulation,triangulation);

  const double left = bottomLeft(0);
  const double bottom = bottomLeft(1);
  const double top = topRight(1);
  const double right = bottomLeft(0) + 3*width;
  setEdgeBoundaries(left,right,top,bottom,triangulation);
}


void MyGridGenerator::generateSquareHole(Triangulation<2,2>& triangulation)
{

  std::vector<unsigned int> holes(2,1);
  GridGenerator::cheese(triangulation,holes);

  printBoundaryVertices(triangulation);

  const double left = 0;
  const double bottom = 0;
  const double top = 3;
  const double right = 3;
  setEdgeBoundaries(left,right,top,bottom,triangulation);
}

void MyGridGenerator::printBoundaryVertices(const Triangulation<2,2>& tria)
{
  std::map<unsigned int, Point<2> > boundaryVertices(dealii::GridTools::get_all_vertices_at_boundary(tria));

  std::cout << "\n Boundary vertex locations: \n";
  for(std::map<unsigned int, Point<2> >::const_iterator it = boundaryVertices.begin();
    it != boundaryVertices.end(); ++it)
    std::cout << it->first << ": " << it->second << std::endl;
}

struct Grid6Func
{
  double trans(const double y) const
  {
    return std::tanh(2*y)/tanh(2);
  }
  Point<2> operator() (const Point<2> &in) const
  {
    return Point<2> (in(0),
                     trans(in(1)));
  }
};

void MyGridGenerator::generateWallAdaptedMesh(Triangulation<2,2>& triangulation)
{
  std::vector< unsigned int > repetitions(2);
  repetitions[0] = repetitions[1] = 40;
  GridGenerator::subdivided_hyper_rectangle (triangulation,
                                             repetitions,
                                             Point<2>(0.0,0.0),
                                             Point<2>(1.0,1.0));
  dealii::GridTools::transform(Grid6Func(), triangulation);


  const double left = 0;
  const double bottom = 0;
  const double top = 1;
  const double right = 1;
  setEdgeBoundaries(left,right,top,bottom,triangulation);
}

*/


/*
  later modify swissCheese for hexagonal sphere packing... need half circles at ends
  for this to work nicely
*/

// void MyGridGenerator::swissCheese(Triangulation<2,2> &tria, const Geometry& geo)
// {

//   Triangulation<2>    square_tile;
//   Triangulation<2>    copy_tria;

//   copy_tria.copy_triangulation(tria);

//   std::string copy_grid_out_file = "copy_tile_grid.eps";
//   std::ofstream out_copy(copy_grid_out_file);
//   GridOut copy_grid_out;
//   copy_grid_out.write_eps(copy_tria,out_copy);

//   Point<2> bottomLeft(2.75, 0);
//   Point<2> topRight(2*2.75, 2.75);

//   GridGenerator::hyper_rectangle(square_tile,
//                 bottomLeft,
//                 topRight);

//   // square_tile.refine_global(1);

//   std::string square_grid_out_file = "square_tile_grid.eps";
//   std::ofstream out_square(square_grid_out_file);
//   GridOut square_grid_out;
//   square_grid_out.write_eps(square_tile,out_square);

//   // merge tria with temp:
//   GridGenerator::merge_triangulations(copy_tria,square_tile,tria); // 1 + 2 = 3

//   // copy triangulation, rotate, shift, and merge again:
//   copy_tria.clear();
//   copy_tria.copy_triangulation(tria);
//   const double angle = numbers::PI;
//   dealii::GridTools::rotate(angle,copy_tria);

//   Tensor<1,2> shift_vector;
//   shift_vector[0] = 2*2.75;
//   shift_vector[1] = 2*2.75;
//   dealii::GridTools::shift(shift_vector,copy_tria);

//   // remerge:
//   GridGenerator::merge_triangulations(copy_tria,tria,tria);

//   // shift and tile:
//   copy_tria.clear();
//   copy_tria.copy_triangulation(tria);

//   // top corner:
//   shift_vector[0] = 2.75*2;
//   shift_vector[1] = 2.75*2;
//   dealii::GridTools::shift(shift_vector,copy_tria);
//   GridGenerator::merge_triangulations(copy_tria,tria,tria);


//   // top:

//   // first shift back and rotate:
//   shift_vector[0] = -2.75*2;
//   shift_vector[1] = -2.75*2;
//   dealii::GridTools::shift(shift_vector,copy_tria);

//   dealii::GridTools::rotate(numbers::PI/2.0,copy_tria);
//   shift_vector[0] = 2.75*2;
//   shift_vector[1] = 2.75*2;

//   dealii::GridTools::shift(shift_vector,copy_tria);
//   GridGenerator::merge_triangulations(copy_tria,tria,tria);

//   // bottom:
//   shift_vector[0] = 2.75*2;;
//   shift_vector[1] = -2.75*2;

//   dealii::GridTools::shift(shift_vector,copy_tria);
//   GridGenerator::merge_triangulations(copy_tria,tria,tria);


// // remove hanging nodes to refine coarse squares:
//   dealii::GridTools::remove_hanging_nodes(tria);

//   std::cout << "Minimal cell diameter: " <<  dealii::GridTools::minimal_cell_diameter(tria)
//     << std::endl;

//   std::cout << "Maximal cell diameter: " << dealii::GridTools::maximal_cell_diameter(tria)
//     << std::endl;

//   // // corner:
//   // shift_vector[1] = 2*2.75;
//   // shift_vector[1] = 0;

//   // dealii::GridTools::shift(shift_vector,copy_tria);
//   // GridGenerator::merge_triangulations(copy_tria,tria,tria);

//   // // bottom:
//   // shift_vector[0] = 0;
//   // shift_vector[1] = -2*2.75;

//   // dealii::GridTools::shift(shift_vector,copy_tria);
//   // GridGenerator::merge_triangulations(copy_tria,tria,tria);

//   // GridGenerator::create_union_triangulation(copy_tria,square_tile,tria);

//   std::string result_grid_out_file = "result_tile_grid.eps";
//   std::ofstream out_result(result_grid_out_file);
//   GridOut result_grid_out;
//   result_grid_out.write_eps(tria,out_result);

//   std::cout << "\n****MERGED TRIANGULATIONS****\n" << std::endl;

//   // Point<2> bottomLeft(geo.getXMin(),geo.getYMin());
//   // Point<2> topRight(geo.getXMax(),geo.getYMax());

//   unsigned int nc = geo.numCircles();
//   for(unsigned int ic = 0; ic < nc; ic++)
//   {
//     const Point<2> center( geo.circleAt(ic)->getCenterX(), geo.circleAt(ic)->getCenterY());
//     const double radius = geo.circleAt(ic)->getRadius();

//     std::cout << "Circle: " << center << "; " << radius << std::endl;
//   }
// } // swiss cheese from geometry


/*
  boxCheese constructs a square packed mesh of holes given by geometry
  boxCheese assumes circles in geometry are square packed and creates a
  base square tile and hole tile given by a circle element in geometry

  might want a function that calculates the bounding box and tile width first...
*/

template<int dim>
void
setEdgeBoundaries(double left, double right,
  double top, double bottom, Triangulation<dim>& triangulation)
{
	if(dim != 2)
		throw std::runtime_error("Function not implemented for dim != 2");

  const unsigned int id_left = 0;
  const unsigned int id_right = 1;
  const unsigned int id_top = 2;
  const unsigned int id_bottom = 3;
  const unsigned int id_other = 5;

  unsigned int left_count = 0;
  unsigned int right_count = 0;
  unsigned int top_count = 0;
  unsigned int bottom_count = 0;

  for (typename Triangulation<2>::active_cell_iterator
       cell = triangulation.begin_active();
     cell != triangulation.end();
     ++cell)
  for (unsigned int f=0; f<dealii::GeometryInfo<2>::faces_per_cell; ++f)
    if (cell->face(f)->at_boundary())
    {
      // std::cout << "yo" << std::endl;
        if (cell->face(f)->center()[0] == left){
          cell->face(f)->set_boundary_id(id_left);
          ++left_count;
        }
        else if (cell->face(f)->center()[0] == right){
          cell->face(f)->set_boundary_id(id_right);
          ++right_count;
        }
        else if (cell->face(f)->center()[1] == top){
          cell->face(f)->set_boundary_id(id_top);
          ++top_count;
        }
        else if (cell->face(f)->center()[1] == bottom){
          cell->face(f)->set_boundary_id(id_bottom);
          ++bottom_count;
        }
        else
        {
          cell->face(f)->set_boundary_id(id_other);
        }
    } // if at boundary

    std::cout << "Left boundaries: " << left_count << std::endl
    << "Right boundaries: " << right_count << std::endl
    << "Top boundaries: " << top_count << std::endl
    << "Bottom boundaries: " << bottom_count << std::endl;
} // setEdgeBoundaries()



template<int dim>
void
setEdgeBoundaries(const Geometry<dim>& geo,
  Triangulation<dim>& triangulation)
{
  if(dim != 2)
    throw std::runtime_error("Function not implemented for dim != 2");

  const double sphere_tolerance = 0.3; // may need to make this better...
  // can maybe make coarser, since coarse mesh shouldn't have a problem

  const double edge_tolerance = 1e-8;

  const unsigned int id_left = 0;
  const unsigned int id_right = 1;
  const unsigned int id_top = 2;
  const unsigned int id_bottom = 3;
  const unsigned int id_other = 5;

  for (typename Triangulation<2>::active_cell_iterator
       cell = triangulation.begin_active();
     cell != triangulation.end();
     ++cell)
  for (unsigned int f=0; f<dealii::GeometryInfo<2>::faces_per_cell; ++f)
    if (cell->face(f)->at_boundary())
    {
      // std::cout << "yo" << std::endl;
        if ( std::fabs( cell->face(f)->center()[0]
          - geo.getBottomLeftPoint()[0] ) < edge_tolerance )
          cell->face(f)->set_boundary_id(id_left);
        else if ( std::fabs( cell->face(f)->center()[0]
          - geo.getTopRightPoint()[0] ) < edge_tolerance )
          cell->face(f)->set_boundary_id(id_right);
        else if ( std::fabs( cell->face(f)->center()[1]
          - geo.getTopRightPoint()[1] ) < edge_tolerance )
          cell->face(f)->set_boundary_id(id_top);
        else if ( std::fabs( cell->face(f)->center()[1]
          - geo.getBottomLeftPoint()[1] ) < edge_tolerance )
          cell->face(f)->set_boundary_id(id_bottom);
        else
          cell->face(f)->set_boundary_id(id_other);
    } // if at boundary

  // set sphere boundaryies:
  unsigned int sphere_bid = 10; // start at 10 ... make global constant...

  for(unsigned int sphere = 0; sphere < geo.getNumberSpheres(); sphere++)
  {
    const Point<2> center = geo.getSphereAt(sphere).getCenter();
    const double radius = geo.getSphereAt(sphere).getRadius();

std::cout << "sphere: " << center << " ; " << radius << " ....." << std::endl;

    for(typename Triangulation<2>::active_cell_iterator
        cell = triangulation.begin_active();
        cell != triangulation.end();
        ++cell)
    {
      for(unsigned int f=0; f < dealii::GeometryInfo<2>::faces_per_cell; ++f)
      {
        const double distance_from_center = center.distance(cell->face(f)->center());

        if( std::fabs(distance_from_center - radius) < sphere_tolerance )
        {
          std::cout << "in sphere: " << sphere_bid << std::endl;
            cell->face(f)->set_boundary_id(sphere_bid);
        }
      } // for faces
    } // for cells

    ++sphere_bid;
  } // for spheres


}



}} // namespace

#endif
