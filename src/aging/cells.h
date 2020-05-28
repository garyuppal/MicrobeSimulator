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
	double getConsumptionRate() const;

	bool isAlive() const;

	void printInfo(std::ostream& out) const;
	void print(std::ostream& out) const;
private:
	std::vector<Point<dim> > locations;

	double consumption_rate;
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
		prm.declare_entry("Consumption rate","0",Patterns::Double());
	prm.leave_subsection();
}

template<int dim>
void
Cells<dim>::init(const ParameterHandler& prm)
{
	const std::string section = "Cells";
	consumption_rate = prm.get_double(section, "Consumption rate");

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
		const double death_rate = dt*(fit.value(*it)); // probability of death
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
double
Cells<dim>::getConsumptionRate() const
{
	return consumption_rate;
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
		<< Utility::medium_line << std::endl
		<< "Consumption rate: " << consumption_rate << std::endl
		<< std::endl << Utility::medium_line << std::endl
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