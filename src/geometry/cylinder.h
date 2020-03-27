#ifndef MICROBESIMULATOR_CYLINDER_H
#define MICROBESIMULATION_CYLINDER_H

/** @file @todo can implement cylinder by inheriting from 2D sphere
 * reflection will then be handled by first projecting 3d coordinate
 * to 2d in appropriate direction and using a circle to reflect
 */

/** Cylinder can be a circle == sphere<2> if 2d or a pipe if 3d
* an enclosing cylinder reflects points on the interior
* whereas an interior cylinder acts as an obstacle that reflects points
* colliding from the outside. 
* A sphere should have the same functionality
* hence this should just be inherited from sphere. We can then also use a sphere
* in 3d to simulation inside a ball for example. 
* 
*
* FOR EXTERIOR BOUNDS:
* In terms of the geometry, bottom left and top right points correspond to the smallest
* box that contains the sphere/cyllinder. functions such as isInDomain return false if in
* space between box and exterior sphere/cylinder. getGridPoints still loops over the 
* outside box, though this will be probably slower than necessary. Since this is mainly
* used for debugging, it shouldn't be a huge issue for now. 
*
* boundary conditions for exterior sphere/cylinder is reflect, and 
* boundary conditions set by geometry shouldn't really matter
* they will matter for 3d cylinder in which case cylinder doesnt reflect along 
* the x axis anyways, and the chemical grid boundary should be handed appropriately
*/

#include <deal.II/base/point.h>
using dealii::Point;

#include <deal.II/base/tensor.h>
using dealii::Tensor;

namespace MicrobeSimulator{

// template<int dim> // automatically 3D, use sphere for 2D
	// or instead of inheriting, have cylinder contain a sphere object?
class Cylinder : public Sphere<2>{
public:
	Cylinder();


private:
	double length;
};

} // CLOSE NAMESPACE
#endif