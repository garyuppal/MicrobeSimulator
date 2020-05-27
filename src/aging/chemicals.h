#pragma once

#include "../utility/parameter_handler.h"

namespace MicrobeSimulator{ namespace Aging{

template<int dim>
class Chemicals{
public:
	Chemicals();

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm);

	void printInfo(std::ostream& out) const;

private:

};

// IMPL
// -------------------------------------------------------------------
template<int dim>
Chemicals<dim>::Chemicals()
{}

template<int dim>
void 
Chemicals<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Chemicals");
		prm.declare_entry("Number chemicals", "1", Patterns::Unsigned());
		prm.declare_entry("Diffusion","{1}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Decay rate","{0}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Save chemicals", "False",Patterns::Bool());
	prm.leave_subsection();
}

template<int dim>
void
Chemicals<dim>::init(const ParameterHandler& prm)
{}

template<int dim>
void
Chemicals<dim>::printInfo(std::ostream& out) const
{
	
}

}} // CLOSE NAMESPACE
/* chemicals.h */