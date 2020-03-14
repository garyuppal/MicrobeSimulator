#ifndef MICROBE_SIMULATOR_NEW_BACTERIA_FITNESS_TYPES_H
#define MICROBE_SIMULATOR_NEW_BACTERIA_FITNESS_TYPES_H

#include <deal.II/base/point.h>
using dealii::Point;

// #include "../chemicals/chemical_interface.h"
#include "../chemicals/chemical_handler.h"

#include <vector>
#include <memory>

namespace MicrobeSimulator{ namespace BacteriaNew{

namespace Fitness{
	void declare_parameters(ParameterHandler& prm)
	{
		prm.enter_subsection("Fitness");
			prm.declare_entry("Fitness type",
								"OR",
								Patterns::Selection("OR|AND|SUM"));
			prm.declare_entry("Chemical fitness","{0,0}",Patterns::List(Patterns::Double()));
			// prm.declare_entry("Chemical harm","{0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Secretion cost","{0,0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Chemical saturation","{1,1}",Patterns::List(Patterns::Double()));
		prm.leave_subsection();
	}
} // fitness namespace
// ----------------------------------------------------------------------
// FITNESS BASE: (Interface)
// ----------------------------------------------------------------------
template<int dim>
class FitnessBase{
public:
	FitnessBase();
	FitnessBase(const Chemicals::ChemicalHandler<dim>& ch);
	virtual ~FitnessBase() {}

	void attach_chemicals(const Chemicals::ChemicalHandler<dim>& ch);

	virtual double value(const Point<dim>& location,
					const std::vector<double>& rates) const =0;

	virtual void setup_fitness_constants(const std::vector<double>& cfit,
									double sc,
									const std::vector<double> csat) = 0;

	virtual void printInfo(std::ostream& out) const;
protected:
	Chemicals::ChemicalHandler<dim> const * chemicals; 
	// std::vector<Chemicals::ChemicalInterface<dim>* > chemicals; // maybe use a chemical handler instead

};

// IMPL
// -----------------------------------------------------------
template<int dim>
FitnessBase<dim>::FitnessBase()
{}

template<int dim>
FitnessBase<dim>::FitnessBase(const Chemicals::ChemicalHandler<dim>& ch)
	:
	chemicals(&ch)
{}

template<int dim>
void 
FitnessBase<dim>::attach_chemicals(const Chemicals::ChemicalHandler<dim>& ch)
{
	chemicals = &ch;
}

template<int dim>
void 
FitnessBase<dim>::printInfo(std::ostream& out) const
{
	out << "Fitness base" << std::endl;
}

// ----------------------------------------------------------------------
// OR FITNESS:
// ----------------------------------------------------------------------
template<int dim>
class OR_Fitness : public FitnessBase<dim>{
public:
	OR_Fitness();
	OR_Fitness(const Chemicals::ChemicalHandler<dim>& ch);
	OR_Fitness(const Chemicals::ChemicalHandler<dim>& ch,
		double bc, double hc, double bs, double hs, double sc);

	double value(const Point<dim>& location,
				const std::vector<double>& rates) const override;

	void setup_fitness_constants(const std::vector<double>& cfit,
								double sc,
								const std::vector<double> csat) override;

	void printInfo(std::ostream& out) const override;
private:
	double benefit_constant;
	double harm_constant;

	double benefit_saturation;
	double harm_saturation;

	double secretion_cost;
};

// IMPL
// -----------------------------------------------------------
template<int dim>
OR_Fitness<dim>::OR_Fitness()
	: 
	benefit_constant(0), harm_constant(0), benefit_saturation(1), harm_saturation(1),
	secretion_cost(0)
{}

template<int dim>
OR_Fitness<dim>::OR_Fitness(
	const Chemicals::ChemicalHandler<dim>& ch)
	:
	FitnessBase<dim>(ch),
	benefit_constant(0), harm_constant(0), benefit_saturation(1), harm_saturation(1),
	secretion_cost(0)	
{}

template<int dim>
OR_Fitness<dim>::OR_Fitness(
	const Chemicals::ChemicalHandler<dim>& ch,
	double bc, double hc, double bs, double hs, double sc)
	:
	FitnessBase<dim>(ch),
	benefit_constant(bc),
	harm_constant(hc),
	benefit_saturation(bs),
	harm_saturation(hs),
	secretion_cost(sc)
{}

template<int dim>
double 
OR_Fitness<dim>::value(const Point<dim>& location,
			const std::vector<double>& rates) const
{
	assert(rates.size() == this->chemicals->getNumberChemicals());
	const unsigned int numchem = this->chemicals->getNumberChemicals();

	double total_goods = 0;
	double total_secretion_rate = 0;

	for(unsigned int i = 0; i < numchem-1; ++i)
	{
		total_goods += (*this->chemicals)[i].value(location);
		total_secretion_rate += rates[i];
	}
	const double waste = (*this->chemicals)[numchem-1].value(location);

	return benefit_constant * total_goods / ( total_goods + benefit_saturation )
		- harm_constant * waste / ( waste + harm_saturation )
		- secretion_cost * total_secretion_rate; 
}

template<int dim>
void 
OR_Fitness<dim>::setup_fitness_constants(const std::vector<double>& cfit,
							double sc,
							const std::vector<double> csat)
{
	benefit_constant = cfit[0];
	harm_constant = -cfit[1];
	secretion_cost = sc;
	benefit_saturation = csat[0];
	harm_saturation = csat[1];
}

template<int dim>
void 
OR_Fitness<dim>::printInfo(std::ostream& out) const
{
	const unsigned int numchem = this->chemicals->getNumberChemicals();

	out << "\n\n-----------------------------------------------------" << std::endl
	    << "\t\tOR FITNESS FUNCTION (for " << numchem << " chemicals)" << std::endl
	    << "-----------------------------------------------------" << std::endl
	    << "\t public good benefit: " << benefit_constant << std::endl
	    << "\t waste harm: " << harm_constant << std::endl
	    << "\t public good saturation: " << benefit_saturation << std::endl
	    << "\t waste saturation: " << harm_saturation << std::endl
	    << "\t secretion cost: " << secretion_cost << std::endl;
}

}} // CLOSE NAMESPACE
#endif