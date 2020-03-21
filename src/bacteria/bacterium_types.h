#ifndef MICROBE_SIMULATOR_BACTERIUM_TYPES_H
#define MICROBE_SIMULATOR_BACTERIUM_TYPES_H

#include <deal.II/base/point.h>
using dealii::Point;

#include "../geometry/geometry.h"
#include "../advection/advection_handler.h"
// #include "../fitness/fitness_base.h"
#include "./bacteria_fitness.h"

#include "../utility/utility.h"

#include <vector>
#include <memory>


namespace MicrobeSimulator{ namespace BacteriaNew{

template<int dim>
class BacteriumBase{
public:
	BacteriumBase();
	BacteriumBase(const Point<dim>& p);
	BacteriumBase(const std::vector<double>& rates);
	BacteriumBase(const Point<dim>& p, 
					const std::vector<double>& rates);

	virtual ~BacteriumBase() {}

	// VIRTUAL METHODS:
// virtual since can be overridden by chemotaxis for example
	virtual void randomStep(double time_step, 
						double diffusion_constant,
						const Geometry<dim>& geo,
						const Velocity::AdvectionHandler<dim>& velocity,
						double buffer=0); 

	virtual double getFitness(const FitnessBase<dim>& fitness_function) const;
		// need new fitness base class -- only dim as templated
	virtual std::unique_ptr<BacteriumBase<dim> > clone() const;
	
	// MUTATORS:
	void setLocation(const Point<dim>& p);
	void setSecretionRates(const std::vector<double> rates);
	void setSecretionRate(unsigned int index, double value);

	// ACCESSORS:
	Point<dim> getLocation() const;
	std::vector<double> getSecretionRates() const;
	double getSecretionRate(unsigned int i) const;
	unsigned int getNumberChemicals() const;

	virtual void print(std::ostream& out) const;
protected:
	Point<dim> location;
	std::vector<double> secretion_rates; // use dynamic array?
};

// IMPL
// ----------------------------------------------------------------
template<int dim>
BacteriumBase<dim>::BacteriumBase()
{}

template<int dim>
BacteriumBase<dim>::BacteriumBase(const Point<dim>& p)
	:
	location(p)
{}

template<int dim>
BacteriumBase<dim>::BacteriumBase(const std::vector<double>& rates)
	:
	location(),
	secretion_rates(rates)
{}

template<int dim>
BacteriumBase<dim>::BacteriumBase(const Point<dim>& p, 
					const std::vector<double>& rates)
	:
	location(p),
	secretion_rates(rates)
{}

/** \brief random walk step */
/** buffer set to 0 as default. Non-zero buffer corresponds to an extra
*	reflection length against solid boundaries 
*/
template<int dim>
void 
BacteriumBase<dim>::randomStep(double time_step, 
						double diffusion_constant,
						const Geometry<dim>& geometry,
						const Velocity::AdvectionHandler<dim>& velocity,
						double buffer)
{
	Point<dim> old_location(location);

	const double theta = 2*dealii::numbers::PI*Utility::getRand();
	const double phi = dealii::numbers::PI*Utility::getRand();

	Point<dim> randomPoint = (dim == 2) ? Point<dim>(std::cos(theta),std::sin(theta))
	                        : Point<dim>(std::cos(phi), std::sin(phi)*std::cos(theta),
	                            std::sin(phi)*std::sin(theta));

	location += std::sqrt(2*dim*time_step*diffusion_constant)*randomPoint
	    + time_step*velocity.value(location);

	geometry.checkBoundaries(old_location, location, buffer); 
}

template<int dim>
double 
BacteriumBase<dim>::getFitness(const FitnessBase<dim>& fitness_function) const
{
	return fitness_function.value(location, secretion_rates);
}

template<int dim>	
std::unique_ptr<BacteriumBase<dim> > 
BacteriumBase<dim>::clone() const
{
	return std::unique_ptr<BacteriumBase<dim> >(
			new BacteriumBase(location,secretion_rates) );
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
BacteriumBase<dim>::setSecretionRates(const std::vector<double> rates)
{
	assert(secretion_rates.size() == rates.size());
	secretion_rates = rates;
}

template<int dim>
void 
BacteriumBase<dim>::setSecretionRate(unsigned int index, double value)
{
	assert(index < secretion_rates.size());
	secretion_rates[index] = value;
}

// ACCESSORS:
template<int dim>
Point<dim> 
BacteriumBase<dim>::getLocation() const
{
	return location;
}

template<int dim>
std::vector<double> 
BacteriumBase<dim>::getSecretionRates() const
{
	return secretion_rates;
}

template<int dim>
double 
BacteriumBase<dim>::getSecretionRate(unsigned int i) const
{
	assert(i < secretion_rates.size());
	return secretion_rates[i];
}

template<int dim>
unsigned int
BacteriumBase<dim>::getNumberChemicals() const
{
	return secretion_rates.size();
}

template<int dim>
void 
BacteriumBase<dim>::print(std::ostream& out) const
{
	out << location << " ";
	for(unsigned int i = 0; i < secretion_rates.size()-1; ++i)
		out << secretion_rates[i] << " ";
	out << secretion_rates[secretion_rates.size()-1];
}

// ------------------------------------------------------------------------------
// IC-AS BACTERIA
// ------------------------------------------------------------------------------

template<int dim>
class IC_Bacterium : public BacteriumBase<dim>{
	IC_Bacterium();
	IC_Bacterium(const Point<dim>& p);
	IC_Bacterium(const std::vector<double>& rates);
	IC_Bacterium(const Point<dim>& p,
				const std::vector<double>& rates);
	IC_Bacterium(const Point<dim>& p,
				const std::vector<double>& rates,
				double t);
	// 	BacteriumBase(const Point<dim>& p);
	// BacteriumBase(const std::vector<double>& rates);
	// BacteriumBase(const Point<dim>& p, 
	// 				const std::vector<double>& rates);


	// VIRTUAL METHODS:
	// virtual since can be overridden by chemotaxis for example
	// virtual void randomStep(double time_step, 
	// 					double diffusion_constant,
	// 					const Geometry<dim>& geo,
	// 					const Velocity::AdvectionHandler<dim>& velocity);  // same as base

	// virtual double getFitness(const FitnessBase<dim>& fitness_function) const; 
//TODO:		// modify???
	std::unique_ptr<BacteriumBase<dim> > clone() const override;

protected:
	double time;
};

// IMPL
// ---------------------------------------------------------------------
template<int dim>
IC_Bacterium<dim>::IC_Bacterium()
	:
	BacteriumBase<dim>(),
	time(0)
{}

template<int dim>
IC_Bacterium<dim>::IC_Bacterium(const Point<dim>& p)
	:
	BacteriumBase<dim>(p),
	time(0)
{}

template<int dim>
IC_Bacterium<dim>::IC_Bacterium(const std::vector<double>& rates)
	:
	BacteriumBase<dim>(rates),
	time(0)
{}

template<int dim>
IC_Bacterium<dim>::IC_Bacterium(const Point<dim>& p,
			const std::vector<double>& rates)
	:
	BacteriumBase<dim>(p,rates),
	time(0)
{}

template<int dim>
IC_Bacterium<dim>::IC_Bacterium(const Point<dim>& p,
			const std::vector<double>& rates,
			double t)
	:
	BacteriumBase<dim>(p,rates),
	time(t)
{}

template<int dim>
std::unique_ptr<BacteriumBase<dim> > 
IC_Bacterium<dim>::clone() const
{
	return std::unique_ptr<BacteriumBase<dim> >(
		new IC_Bacterium(this->location,this->secretion_rates,time));
}


}} // CLOSE NAMESPACES
#endif