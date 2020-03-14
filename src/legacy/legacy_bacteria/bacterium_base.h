#ifndef MICROBESIMULATOR_BACTERIUM_BASE_H
#define MICROBESIMULATOR_BACTERIUM_BASE_H

#include <deal.II/base/point.h>

#include "../geometry/geometry.h"
#include "../advection/advection_handler.h"
#include "../fitness/fitness_base.h"

#include <iostream>
#include <string>
#include <array>

namespace MicrobeSimulator{ namespace MultiBacteria{

/** Bacterium base
* forms the foundation for normal, intermitently cheating, and altruisitic suicide bacteria...
* could also not include numchem and have specialists inherit from this type...
* for now we let number of chemicals stay a template parameter
*
* This class on its own random walk in a given geometry, ... 
*/
template<int dim>
class BacteriumBase{
public:
	BacteriumBase();
	BacteriumBase(double pg, double w);
	BacteriumBase(const Point<dim>& p, double pg, double w);
	virtual ~BacteriumBase() {}

	void randomStep(double time_step, double diffusion_constant,
		const Geometry<dim>& geometry, 
		const Velocity::AdvectionHandler<dim>& velocity); // done

	void setLocation(const Point<dim>& p); // done
	void setGoodSecretion(double pg); // done
	void setWasteSecretion(double w); // done

	// accessors:
	virtual double getFitness(const FitnessBase<dim, 2>& fitness_function) const; // done

	Point<dim> getLocation() const; // done
	double getGoodSecretion() const; // done
	double getWasteSecretion() const; // done

	virtual void print(std::ostream& out) const; 

protected:
	Point<dim> location;
	double public_good_secretion;
	double waste_secretion;
};

// IMPLEMENTATION
// ---------------------------------------------------------------------
template<int dim>
BacteriumBase<dim>::BacteriumBase()
	:
	location(),
	public_good_secretion(0),
	waste_secretion(0)
{}

template<int dim>
BacteriumBase<dim>::BacteriumBase(double pg, double w)
	:
	public_good_secretion(pg),
	waste_secretion(w)
{}

template<int dim>
BacteriumBase<dim>::BacteriumBase(const Point<dim>& p, double pg, double w)
	:
	location(p),
	public_good_secretion(pg),
	waste_secretion(w)
{}

template<int dim>
void
BacteriumBase<dim>::randomStep(double time_step, double diffusion_constant,
		const Geometry<dim>& geometry, const Velocity::AdvectionHandler<dim>& velocity)
{
	Point<dim> old_location = location;

	const double theta = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));
	const double phi = dealii::numbers::PI*((double) rand() / (RAND_MAX));

	Point<dim> randomPoint = (dim == 2) ? Point<dim>(std::cos(theta),std::sin(theta))
	                        : Point<dim>(std::cos(phi), std::sin(phi)*std::cos(theta),
	                            std::sin(phi)*std::sin(theta));

	location += std::sqrt(2*dim*time_step*diffusion_constant)*randomPoint
	    + time_step*velocity.value(location);

	geometry.checkBoundaries(old_location, location); 
}

// virtual:
template<int dim>
double 
BacteriumBase<dim>::getFitness(const FitnessBase<dim, 2>& fitness_function) const
{
	return fitness_function.value(location, 
		std::array<double, 2>({public_good_secretion, waste_secretion}) );
}

// MUTATORS:
template<int dim>
void 
BacteriumBase<dim>::setLocation(const Point<dim>& p)
{
	location = p;
}

template<int dim>
void 
BacteriumBase<dim>::setGoodSecretion(double pg)
{
	public_good_secretion = pg;
}

template<int dim>
void 
BacteriumBase<dim>::setWasteSecretion(double w)
{
	waste_secretion = w;
}

// ACCESSORS:
template<int dim>
Point<dim> 
BacteriumBase<dim>::getLocation() const
{
	return location;
}

template<int dim>
double 
BacteriumBase<dim>::getGoodSecretion() const
{
	return public_good_secretion;
}

template<int dim>
double 
BacteriumBase<dim>::getWasteSecretion() const
{
	return waste_secretion;
}

template<int dim>
void
BacteriumBase<dim>::print(std::ostream& out) const
{
    out << location << " "
    	<< public_good_secretion << " "
    	<< waste_secretion << std::endl;
}

}} // close namespaces
#endif