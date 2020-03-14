#ifndef SPHERE_H
#define SPHERE_H

#include <deal.II/base/point.h>
using dealii::Point;

namespace MicrobeSimulator{
  template<int dim>
  class Sphere{
  private:
    Point<dim> center;
    double radius;
  public:
    Sphere();
    Sphere(Point<dim> c, double r);
    Sphere(const Sphere &s);

    // accessors:
    Point<dim> getCenter() const;
    double getRadius() const;

    // mutators:
    void setCenter(Point<dim> c);
    void setRadius(double r);

  }; // class Sphere{}

// IMPLEMENTATION
// -------------------------------------------------------------------
    template<int dim>
  Sphere<dim>::Sphere()
    : center(), radius(0.0)
    {}

  template<int dim>
  Sphere<dim>::Sphere(Point<dim> c, double r)
    : center(c), radius(r)
    {}

  template<int dim>
  Sphere<dim>::Sphere(const Sphere &s)
  {
    center = s.center;
    radius = s.radius;
  }

  // accessors:
  template<int dim>
  Point<dim> Sphere<dim>::getCenter() const
  {
    return center;
  }

  template<int dim>
  double Sphere<dim>::getRadius() const
  {
    return radius;
  }

  // mutators:
  template<int dim>
  void Sphere<dim>::setCenter(Point<dim> c)
  {
    center = c;
  }

  template<int dim>
  void Sphere<dim>::setRadius(double r)
  {
    radius = r;
  }


}

#endif 
