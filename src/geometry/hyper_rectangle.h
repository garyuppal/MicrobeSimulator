#ifndef MICROBESIMULATOR_HYPER_RECTANGLE_H
#define MICROBESIMULATOR_HYPER_RECTANGLE_H

#include <deal.II/base/point.h>
using dealii::Point;

#include <deal.II/base/tensor.h>
using dealii::Tensor;

namespace MicrobeSimulator{

  template<int dim>
  class HyperRectangle{
  private:
    Point<dim> bottom_left;
    Point<dim> top_right;
  public:
    HyperRectangle();
    HyperRectangle(const Point<dim>& lower,
      const Point<dim>& upper);
    HyperRectangle(const HyperRectangle& rect);

    // accessors:
    Point<dim> getBottomLeft() const;
    Point<dim> getTopRight() const;

    // mutators:
    void setBottomLeft(const Point<dim>& bl);
    void setTopRight(const Point<dim>& tr);

    double distance_from_border(const Point<dim>& p) const;

    Tensor<1, dim> getNormalVector(const Point<dim>& p) const;

    void reflectPoint(const Point<dim>& old_point,
                      Point<dim>& new_point,
                      const double buffer = 0.) const;

  }; // class Sphere{}

// IMPLEMENTATION
// -------------------------------------------------------------------
  template<int dim>
  HyperRectangle<dim>::HyperRectangle()
  {}

  template<int dim>
  HyperRectangle<dim>::HyperRectangle(const Point<dim>& lower,
      const Point<dim>& upper)
    :
    bottom_left(lower),
    top_right(upper)
  {}

  template<int dim>
  HyperRectangle<dim>::HyperRectangle(const HyperRectangle& rect)
  {
    bottom_left = rect.bottom_left;
    top_right = rect.top_right;
  }

  // accessors:
  template<int dim>
  Point<dim> 
  HyperRectangle<dim>::getBottomLeft() const
  {
    return bottom_left;
  }

  template<int dim>
  Point<dim> 
  HyperRectangle<dim>::getTopRight() const
  {
    return top_right;
  }

  // mutators:
  template<int dim>
  void 
  HyperRectangle<dim>::setBottomLeft(const Point<dim>& bl)
  {
    bottom_left = bl;
  }
  
  template<int dim>
  void 
  HyperRectangle<dim>::setTopRight(const Point<dim>& tr)
  {
    top_right = tr;
  }

  template<int dim>
  double 
  HyperRectangle<dim>::distance_from_border(const Point<dim>& p) const
  {
    /// find minimal distance outside rectangle
    double distance = 0.; 

    Point<dim> center = 0.5*(top_right + bottom_left);
    Point<dim> half_width = 0.5*(top_right + (-1.)*bottom_left);

    for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
    {
      const double dx = std::max( std::fabs(p[dim_itr] - center[dim_itr]) -
          half_width[dim_itr], 0.);

      distance += (dx*dx);
    }

    return distance;

//     dx = max(abs(px - x) - width / 2, 0);
// dy = max(abs(py - y) - height / 2, 0);
// return dx * dx + dy * dy;
  }


  template<int dim>
  void 
  HyperRectangle<dim>::reflectPoint(const Point<dim>& old_point,
                                  Point<dim>& new_point,
                                  const double buffer) const
  {
    const double edge_tolerance = 1e-8 + buffer; 
    
    /// @ todo implement buffer for rectangle reflect point
    for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
    {
      /// check if point came from lower or above end along this dimension
      const bool from_below = (old_point[dim_itr] < (bottom_left[dim_itr] - edge_tolerance) ) && 
                              (new_point[dim_itr] > (bottom_left[dim_itr] - edge_tolerance) ); 

      const bool from_above = (new_point[dim_itr] < (top_right[dim_itr] + edge_tolerance) ) &&
                              (old_point[dim_itr] > (top_right[dim_itr] + edge_tolerance) );

      if(from_below) // including from left
      { 
        // delta should already be positive, but in case not, use fabs:
        double use_buffer = buffer;
        if(dim_itr == 0)
          use_buffer = 0; // dont add buffer to left (x = dim 0)
        const double delta = std::fabs(new_point[dim_itr] - bottom_left[dim_itr]) + 0.5*use_buffer;
        new_point[dim_itr] = new_point[dim_itr] - 2.*delta;
      }
      else if(from_above)
      {
        const double delta = std::fabs(top_right[dim_itr] - new_point[dim_itr]) + 0.5*buffer;
        new_point[dim_itr] = new_point[dim_itr] + 2.*delta;
      }
      // else do nothing
    } // for each dimension

  } // reflect point incident on rectangle


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
