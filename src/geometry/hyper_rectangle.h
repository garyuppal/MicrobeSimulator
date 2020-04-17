#ifndef MICROBESIMULATOR_HYPER_RECTANGLE_H
#define MICROBESIMULATOR_HYPER_RECTANGLE_H

#include <deal.II/base/point.h>
using dealii::Point;

#include <deal.II/base/tensor.h>
using dealii::Tensor;

namespace MicrobeSimulator{

// if we implement ``lines'' or hyper-plane in general, rectangle can also be constructed from these...

/** \brief HyperRectangle class for interior rectangular obstacles */
template<int dim>
class HyperRectangle{
public:
    HyperRectangle();
    HyperRectangle(const Point<dim>& lower,
      const Point<dim>& upper);

    HyperRectangle(const Point<dim>& lower,
  				const Point<dim>& upper, 
  				double theta);

    // accessors:
    Point<dim> getBottomLeft() const;
    Point<dim> getTopRight() const;
    double getAngleXY() const;

    // mutators:
    void setBottomLeft(const Point<dim>& bl);
    void setTopRight(const Point<dim>& tr);
    void setAngleXY(double theta);

    double distance_from_border(const Point<dim>& p, double buffer=0) const; // add angle***

    Tensor<1, dim> getNormalVector(const Point<dim>& p) const; // add angle***

    void reflectPoint(const Point<dim>& old_point,
                      Point<dim>& new_point,
                      const double buffer = 0.) const; // add angle*** could just switch coordinates
private:
    Point<dim> bottom_left;
    Point<dim> top_right;

    double 		angle_xy; /**< Angle in X-Y plane measured from x axis */

    Point<dim> rotate_point() const; // for using angle..., about origin??? or about bottom left point?

}; // class HyperRectangle{}


// IMPLEMENTATION
// -------------------------------------------------------------------

/** \brief Default constructor */
template<int dim>
HyperRectangle<dim>::HyperRectangle()
	:
	angle_xy(0)
{}

/** \brief Constuctor given endpoints */
template<int dim>
HyperRectangle<dim>::HyperRectangle(const Point<dim>& lower,
									const Point<dim>& upper)
	:
	bottom_left(lower),
	top_right(upper),
	angle_xy(0)
{}

/** \brief Constuctor given endpoints and angle */
template<int dim>
HyperRectangle<dim>::HyperRectangle(const Point<dim>& lower,
									const Point<dim>& upper,
									double theta)
	:
	bottom_left(lower),
	top_right(upper),
	angle_xy(theta)
{}

// accessors:

/** \brief Return bottom left corner of rectangle */
template<int dim>
Point<dim> 
HyperRectangle<dim>::getBottomLeft() const
{
	return bottom_left;
}

/** \brief Return top right corner of rectangle */
template<int dim>
Point<dim> 
HyperRectangle<dim>::getTopRight() const
{
	return top_right;
}

/** \brief Return x-y plane angle */
template<int dim>
double 
HyperRectangle<dim>::getAngleXY() const
{
	return angle_xy;	
}


// mutators: (do we use these? need also to check validity... )

/** \brief Set bottom left point of rectangle */
template<int dim>
void 
HyperRectangle<dim>::setBottomLeft(const Point<dim>& bl)
{
	bottom_left = bl;
}

/** \brief Set top right point of rectangle */
template<int dim>
void 
HyperRectangle<dim>::setTopRight(const Point<dim>& tr)
{
	top_right = tr;
}

/** \brief Set XY plane angle */
template<int dim>
void
HyperRectangle<dim>::setAngleXY(double theta)
{
	angle_xy = theta;
}



/** \brief Return distance to outer boundary of rectangle */
template<int dim>
double 
HyperRectangle<dim>::distance_from_border(const Point<dim>& p, double buffer) const
{
	// find minimal distance outside rectangle
	double distance = 0.; 

	Point<dim> center = 0.5*(top_right + bottom_left);
	Point<dim> half_width = 0.5*(top_right + (-1.)*bottom_left);

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
	{
		const double dx = std::max( 
			std::fabs(p[dim_itr] - center[dim_itr]) - (half_width[dim_itr] + buffer), 0.);

		distance += (dx*dx);
	}

	/** @todo double check ... */
	return std::sqrt(distance);
}

/** \brief Reflect point off of rectangle
*/
template<int dim>
void 
HyperRectangle<dim>::reflectPoint(const Point<dim>& old_point,
                              Point<dim>& new_point,
                              const double buffer) const
{	
    for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
    {
    	// buffer should be consistent with check
	    const double edge_tolerance = 1e-8 + buffer; 

	    // dont add buffer to left (x = dim 0 and from below):
	    // double fb_edge_tolerance = edge_tolerance;
	    // if(dim_itr == 0) 
	    // 	fb_edge_tolerance = 1e-8;

		const bool from_below = (old_point[dim_itr] < (bottom_left[dim_itr] - edge_tolerance) ) && 
		                  (new_point[dim_itr] > (bottom_left[dim_itr] - edge_tolerance) ); 

		const bool from_above = (new_point[dim_itr] < (top_right[dim_itr] + edge_tolerance) ) &&
		                  (old_point[dim_itr] > (top_right[dim_itr] + edge_tolerance) );

		if(from_below) // including from left
		{ 
			// delta should already be positive, but in case not, use fabs:
			const double delta = std::fabs(new_point[dim_itr] - bottom_left[dim_itr] + edge_tolerance) ;
			new_point[dim_itr] = new_point[dim_itr] - 2.*delta;
		}
		else if(from_above)
		{
			const double delta = std::fabs(top_right[dim_itr] + edge_tolerance - new_point[dim_itr]);
			new_point[dim_itr] = new_point[dim_itr] + 2.*delta;
		}
		// else do nothing 
	} // for each dimension
} // reflect point incident on rectangle

/** \brief Get unit vector normal to rectangle in direction of supplied point */
template<int dim>
Tensor<1, dim> 
HyperRectangle<dim>::getNormalVector(const Point<dim>& p) const
{
	Tensor<1, dim> normal;

	// find direction:
	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		if( p[dim_itr] >= top_right[dim_itr] )
			normal[dim_itr] = 1.;
		else if( p[dim_itr] <= bottom_left[dim_itr])
			normal[dim_itr] = -1.;

	normal /= normal.norm();
	return normal; 
}


} // close namespace

#endif 
