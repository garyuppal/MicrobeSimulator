#pragma once

#include <deal.II/base/point.h>
using dealii::Point;

#include "../utility/parameter_handler.h"
#include "./fitness.h"

namespace MicrobeSimulator{ namespace Aging{

// template<int dim>
// std::vector<Point<dim> >
// getCellLocations(unsigned int n, Point<dim> lower, Point<dim>, upper);

// template<>
std::vector<Point<1> >
getCellLocations(unsigned int n, Point<1> lower, Point<1> upper)
{
	const double length = upper[0] - lower[0];

	const double dx = length/((double)n);

	std::vector<Point<1> > locs;
	locs.reserve(n);
	for(double x = 0; x <= length; x += dx)
		locs.emplace_back(Point<1>(x));

	return locs;
}


template<int dim>
class Cells{
public:
	Cells();

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm);

	void update(double dt, const Aging::Fitness<dim>& fit);

	std::vector<Point<dim> > getLocations() const;
	
	// secretion/consumption:

	// double getConsumptionRate() const;
	std::vector<double> getSecretionRates() const;
	double getSecretionRate(unsigned int i) const;

	std::vector<double> getConsumptionRates() const;
	double getConsumptionRate(unsigned int i) const;

	// checking/output
	unsigned int getSize() const;
	bool isAlive() const;

	void printInfo(std::ostream& out) const;
	void print(std::ostream& out) const;
private:
	std::vector<Point<dim> > locations;

	std::vector<double> secretion_rates;
	std::vector<double> consumption_rates;
};

// IMPL
// -------------------------------------------------------------------
template<int dim>
Cells<dim>::Cells()
{}

template<int dim>
void 
Cells<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Cells");
		prm.declare_entry("Number cells", "1", Patterns::Unsigned());
		prm.declare_entry("Save cells", "True",Patterns::Bool());
		prm.declare_entry("Save counts", "False",Patterns::Bool());
		
		prm.declare_entry("Secretion rates",
				"{0}",
				Patterns::List(Patterns::Double()));
		prm.declare_entry("Consumption rates",
				"{0}",
				Patterns::List(Patterns::Double()));
	prm.leave_subsection();
}

template<int dim>
void
Cells<dim>::init(const ParameterHandler& prm)
{
	const std::string section = "Cells";

	secretion_rates = prm.get_double_vector(section,"Secretion rates");
	consumption_rates = prm.get_double_vector(section,"Consumption rates");

	const unsigned int number_cells = prm.get_unsigned(section, "Number cells");

	// uniformly spread out cells:
	std::vector<Point<dim> > bounds = prm.get_oneD_point_list("Chemicals", "Bounds");
	Point<dim> lower = bounds[0];
	Point<dim> upper = bounds[1];
	locations = getCellLocations(number_cells, lower, upper);
}

template<int dim>
void 
Cells<dim>::update(double dt, const Aging::Fitness<dim>& fit)
{
	for(auto it = locations.begin(); it != locations.end(); )
	{
		const double death_rate = dt*std::fabs(fit.value(*it)); // probability of death
		double prob = Utility::getRand();

		if(prob < death_rate )
			it = locations.erase(it);
		else
			++it;
	}
}

template<int dim>
std::vector<Point<dim> > 
Cells<dim>::getLocations() const
{
	return locations;
}

template<int dim>
std::vector<double> 
Cells<dim>::getSecretionRates() const
{
	return secretion_rates;
}

template<int dim>
double 
Cells<dim>::getSecretionRate(unsigned int i) const
{
	// assert(i < secretion_rates.size());
	return secretion_rates[i];
}

template<int dim>
std::vector<double> 
Cells<dim>::getConsumptionRates() const
{
	return consumption_rates;
}

template<int dim>
double 
Cells<dim>::getConsumptionRate(unsigned int i) const
{
	// assert(i < consumption_rates.size());
	return consumption_rates[i];
}

template<int dim>
unsigned int 
Cells<dim>::getSize() const
{
	return locations.size();
}

template<int dim>
bool 
Cells<dim>::isAlive() const
{
	return !locations.empty();
}


template<int dim>
void
Cells<dim>::printInfo(std::ostream& out) const
{
	out << "\n\n" << Utility::medium_line << std::endl 
		<< "\t\t CELLS INFO:" << std::endl
		<< Utility::medium_line << std::endl;

	out << "Secretion rates: ";
	for(unsigned int i = 0; i < secretion_rates.size(); ++i)
		out << secretion_rates[i] << " ";
	out << std::endl;

	out << "Consumption rates: ";
	for(unsigned int i = 0; i < consumption_rates.size(); ++i)
		out << consumption_rates[i] << " ";
	out << std::endl;

	out << std::endl << Utility::medium_line << std::endl
		<< std::endl << std::endl;
}

template<int dim>
void
Cells<dim>::print(std::ostream& out) const
{
	const unsigned int n = locations.size();
	for(unsigned int i = 0; i < n; ++i)
		out << locations[i] << std::endl;
}

}} // CLOSE NAMESPACE
/* chemicals.h */