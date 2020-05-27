#pragma once

#include <deal.II/base/point.h>
using dealii::Point;

#include "../utility/parameter_handler.h"

namespace MicrobeSimulator{ namespace Aging{

template<int dim>
class Cells{
public:
	Cells();

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm);

	bool isAlive() const;

	void printInfo(std::ostream& out) const;
	void print(std::ostream& out) const;
private:
	std::vector<Point<dim> > locations;
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
	prm.leave_subsection();
}

template<int dim>
void
Cells<dim>::init(const ParameterHandler& prm)
{

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