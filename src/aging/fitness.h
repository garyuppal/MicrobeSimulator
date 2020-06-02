#pragma once

#include <deal.II/base/point.h>
using dealii::Point;

#include "../utility/utility.h"
#include "../utility/parameter_handler.h"
#include "./chemicals.h"

namespace MicrobeSimulator{ namespace Aging{

// ------------------------------------------------------------------
// FITNESS BASE:
// ------------------------------------------------------------------ 
template<int dim>
class FitnessBase{
public:
	FitnessBase(const ParameterHandler& prm, const Aging::Chemicals<dim>& ch);
	virtual ~FitnessBase() {}

	virtual double value(const Point<dim>& location) const=0;

	virtual void printInfo(std::ostream& out) const=0;

protected:
	const Aging::Chemicals<dim> * const chemicals; 
};

template<int dim>
FitnessBase<dim>::FitnessBase(const ParameterHandler& /*prm*/,
								 const Aging::Chemicals<dim>& ch)
	:
	chemicals(&ch)
{}

// ------------------------------------------------------------------
// Cooperative factors FITNESS:
// ------------------------------------------------------------------ 
template<int dim>
class CF_Fitness : public FitnessBase<dim>{
public:
	CF_Fitness(const ParameterHandler& prm, const Aging::Chemicals<dim>& ch);

	static void declare_parameters(ParameterHandler& prm);

	double value(const Point<dim>& location) const override;
	void printInfo(std::ostream& out) const override;
private:
	double death_rate;
	double saturation_constant;
	double hill_constant;
};

template<int dim>
CF_Fitness<dim>::CF_Fitness(const ParameterHandler& prm, const Aging::Chemicals<dim>& ch)
	:
	FitnessBase<dim>(prm, ch)
{
	const std::string section = "Fitness.CF";
	death_rate = prm.get_double(section, "Death rate");
	saturation_constant = prm.get_double(section, "Saturation constant");
	hill_constant = prm.get_double(section, "Hill constant");
}

template<int dim>
void 
CF_Fitness<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Fitness");
		prm.enter_subsection("CF");
			prm.declare_entry("Death rate", "0", Patterns::Double());
			prm.declare_entry("Saturation constant", "1", Patterns::Double());
			prm.declare_entry("Hill constant","1",Patterns::Double());
		prm.leave_subsection();
	prm.leave_subsection();
}

template<int dim>
double 
CF_Fitness<dim>::value(const Point<dim>& location) const
{
	const double f = (*(this->chemicals))[0].value(location);

	// const double rval = death_rate/(1. + std::pow( f/saturation_constant ,hill_constant) );
	// std::cout << "chemical: " << f << " return: " << rval << std::endl;

	return death_rate/(1. + std::pow( f/saturation_constant ,hill_constant) );
}

template<int dim>
void 
CF_Fitness<dim>::printInfo(std::ostream& out) const
{
	out << "\n\n" << Utility::medium_line << std::endl
	    << "\t\tCF FITNESS FUNCTION" << std::endl
	    << Utility::medium_line << std::endl
	    << "\t Death rate: " << death_rate << std::endl
	    << "\t Saturation constant: " << saturation_constant << std::endl
	    << "\t Hill constant: " << hill_constant << std::endl
		<< Utility::medium_line
		<< std::endl << std::endl << std::endl;
}

// ------------------------------------------------------------------
// FITNESS:
// ------------------------------------------------------------------ 
template<int dim>
class Fitness{
public:
	Fitness();
	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm, const Aging::Chemicals<dim>& chem);

	double value(const Point<dim>& p) const;

	void printInfo(std::ostream& out) const;
private:
	std::shared_ptr<FitnessBase<dim> > 		fit;
};

// IMPL
// ------------------------------------------------------------------

template<int dim>
Fitness<dim>::Fitness()
{}

template<int dim>
void 
Fitness<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Fitness");
		prm.declare_entry("Fitness type", "CF", Patterns::Selection("CF"));
	prm.leave_subsection();

	CF_Fitness<dim>::declare_parameters(prm);
}

template<int dim>
void 
Fitness<dim>::init(const ParameterHandler& prm, const Aging::Chemicals<dim>& chem)
{
	const std::string section = "Fitness";
	const std::string fitness_type = prm.get_string(section, "Fitness type");

	if ( boost::iequals(fitness_type, "CF") )
		fit = std::make_shared<CF_Fitness<dim> >(prm, chem);
	else
		throw std::runtime_error("Invalid fitness type: <" + fitness_type + ">");
}

template<int dim>
double 
Fitness<dim>::value(const Point<dim>& p) const
{
	return fit->value(p);
}

template<int dim>
void 
Fitness<dim>::printInfo(std::ostream& out) const
{
	fit->printInfo(out);
}

}} // CLOSE NAMESPACE
/* fitness.h */