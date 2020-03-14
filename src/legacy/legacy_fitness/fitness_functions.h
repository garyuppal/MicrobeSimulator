#ifndef MICROBESIMULATOR_FITNESS_FUNCTIONS_H
#define MICROBESIMULATOR_FITNESS_FUNCTIONS_H

#include "./fitness_base.h"
#include "../chemicals/chemical_interface.h"

namespace MicrobeSimulator{ namespace FitnessFunctions{

/** Fitness function for two chemicals,
* 1st being public good and 2nd the waste chemical
*/

template<int dim>
class TwoChemFitness : public FitnessBase<dim, 2>{
public:
	TwoChemFitness();

	virtual double value(const Point<dim>& location, 
		const std::array<double, 2>& secretion_rates) const;

	void attach_chemicals(const Chemicals::ChemicalInterface<dim>& c1, 
						const Chemicals::ChemicalInterface<dim>& c2);

	void set_fitness_constants(double benefit, double harm, double good_sat,
		double waste_sat, double cost);

	void printInfo(std::ostream& out) const;

private:
	const Chemicals::ChemicalInterface<dim>*  good_chemical;
	const Chemicals::ChemicalInterface<dim>*  waste_chemical;

	double 		benefit_constant;
	double 		harm_constant;

	double		benefit_saturation;
	double 		harm_saturation;

	double 		secretion_cost;
};

// IMPLEMENTATION:
//-----------------------------------------------------------------------------------------------
template<int dim>
TwoChemFitness<dim>::TwoChemFitness() 
{}

template<int dim>
double 
TwoChemFitness<dim>::value(const Point<dim>& location, 
		const std::array<double, 2>& secretion_rates) const
{
	const double goods = good_chemical->value(location);
	const double waste = waste_chemical->value(location);

	return benefit_constant * goods / ( goods + benefit_saturation )
		- harm_constant * waste / ( waste + harm_saturation )
		- secretion_cost * secretion_rates[0]; 
}

template<int dim>
void 
TwoChemFitness<dim>::attach_chemicals(const Chemicals::ChemicalInterface<dim>& c1, 
							const Chemicals::ChemicalInterface<dim>& c2)
{
	good_chemical = &c1;
	waste_chemical = &c2;
}

template<int dim>
void 
TwoChemFitness<dim>::set_fitness_constants(double benefit, double harm, double good_sat,
	double waste_sat, double cost)
{
	benefit_constant = benefit;
	harm_constant = harm;
	benefit_saturation = good_sat;
	harm_saturation = waste_sat;
	secretion_cost = cost;
}

template<int dim>
void 
TwoChemFitness<dim>::printInfo(std::ostream& out) const
{
  out << "\n\n-----------------------------------------------------" << std::endl
    << "\t\t FITNESS FUNCTION FOR 2 CHEMICALS" << std::endl
    << "-----------------------------------------------------" << std::endl
    << "\t public good benefit: " << benefit_constant << std::endl
    << "\t waste harm: " << harm_constant << std::endl
    << "\t public good saturation: " << benefit_saturation << std::endl
    << "\t waste saturation: " << harm_saturation << std::endl
    << "\t secretion cost: " << secretion_cost << std::endl;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------









/** A general OR type fitness function
* benefit is proportional to sum of (numchem-1) chemical values
* waste chemical is taken as last chemical in array
* secretion cost is given as same for first (numchem-1) rates
*/


template<int dim, int numchem>
class OR_Fitness : public FitnessBase<dim, numchem>
{
public:
	OR_Fitness() {}
	~OR_Fitness() {}

	virtual double value(const Point<dim>& location,
						const std::array<double, numchem>& secretion_rates) const;

	void attach_chemicals(std::array<Chemicals::ChemicalInterface<dim>*, numchem> chem_ptrs);

	void set_fitness_constants(
		double bene, double harm, double bsat, double hsat, double cost);

	void printInfo(std::ostream& out) const;

private:
	std::array<Chemicals::ChemicalInterface<dim>*, numchem>		chemicals;
	
	double 		benefit_constant;
	double 		harm_constant;

	double		benefit_saturation;
	double 		harm_saturation;

	double 		secretion_cost;
};

// IMPLEMENTATION:
//---------------------------------------------------------------------------

template<int dim, int numchem>
void 
OR_Fitness<dim,numchem>::attach_chemicals(std::array<Chemicals::ChemicalInterface<dim>*, numchem> chem_ptrs)
{
	chemicals = chem_ptrs;
}


template<int dim, int numchem>
void 
OR_Fitness<dim, numchem>::set_fitness_constants(
	double bene, double harm, double bsat, double hsat, double cost)
{
	benefit_constant = bene;
	harm_constant = harm;

	benefit_saturation = bsat;
	harm_saturation = hsat;

	secretion_cost = cost;
}


template<int dim, int numchem>
void 
OR_Fitness<dim, numchem>::printInfo(std::ostream& out) const
{
  out << "\n\n-----------------------------------------------------" << std::endl
    << "\t\tOR FITNESS FUNCTION (for " << numchem << " chemicals)" << std::endl
    << "-----------------------------------------------------" << std::endl
    << "\t public good benefit: " << benefit_constant << std::endl
    << "\t waste harm: " << harm_constant << std::endl
    << "\t public good saturation: " << benefit_saturation << std::endl
    << "\t waste saturation: " << harm_saturation << std::endl
    << "\t secretion cost: " << secretion_cost << std::endl;
}


template<int dim, int numchem>
double 
OR_Fitness<dim, numchem>::value(const Point<dim>& location, 
        const std::array<double, numchem>& secretion_rates) const
{
	double total_goods = 0.;
	double total_secretion_rate = 0.;

	for(unsigned int i = 0; i < numchem-1; ++i)
	{
		total_goods += chemicals[i]->value(location);
		total_secretion_rate += secretion_rates[i];
	} // for all public goods
	const double waste = chemicals[numchem-1]->value(location);

	// const double return_value =
	return	benefit_constant * total_goods / ( total_goods + benefit_saturation )
		- harm_constant * waste / ( waste + harm_saturation )
		- secretion_cost * total_secretion_rate; 

	// return return_value;
}



//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------







/** A general AND type fitness function
* benefit is proportional to product of (numchem-1) chemical values
* waste chemical is taken as last chemical in array
* secretion cost is given as same for first (numchem-1) rates
*/

template<int dim, int numchem>
class AND_Fitness : public FitnessBase<dim, numchem>
{
public:
	AND_Fitness() {}
	~AND_Fitness() {}

	virtual double value(const Point<dim>& location,
						const std::array<double, numchem>& secretion_rates) const;

	void attach_chemicals(std::array<Chemicals::ChemicalInterface<dim>*, numchem> chem_ptrs);

	void set_fitness_constants(
		double bene, double harm, double bsat, double hsat, double cost);

	void printInfo(std::ostream& out) const;

private:
	std::array<Chemicals::ChemicalInterface<dim>*, numchem>		chemicals;

	double 		benefit_constant;
	double 		harm_constant;

	double		benefit_saturation;
	double 		harm_saturation;

	double 		secretion_cost;
};

// IMPLEMENTATION:
//---------------------------------------------------------------------------

template<int dim, int numchem>
void 
AND_Fitness<dim,numchem>::attach_chemicals(std::array<Chemicals::ChemicalInterface<dim>*, numchem> chem_ptrs)
{
	chemicals = chem_ptrs;
}


template<int dim, int numchem>
void 
AND_Fitness<dim, numchem>::set_fitness_constants(
	double bene, double harm, double bsat, double hsat, double cost)
{
	benefit_constant = bene;
	harm_constant = harm;

	benefit_saturation = bsat;
	harm_saturation = hsat;

	secretion_cost = cost;
}


template<int dim, int numchem>
void 
AND_Fitness<dim, numchem>::printInfo(std::ostream& out) const
{
  out << "\n\n-----------------------------------------------------" << std::endl
    << "\t\tAND FITNESS FUNCTION (for " << numchem << " chemicals)" << std::endl
    << "-----------------------------------------------------" << std::endl
    << "\t public good benefit: " << benefit_constant << std::endl
    << "\t waste harm: " << harm_constant << std::endl
    << "\t public good saturation: " << benefit_saturation << std::endl
    << "\t waste saturation: " << harm_saturation << std::endl
    << "\t secretion cost: " << secretion_cost << std::endl;
}


template<int dim, int numchem>
double 
AND_Fitness<dim, numchem>::value(const Point<dim>& location, 
        const std::array<double, numchem>& secretion_rates) const
{
	double total_goods = 1.;
	double total_secretion_rate = 0.;

	for(unsigned int i = 0; i < numchem-1; ++i)
	{
		total_goods *= chemicals[i]->value(location);
		total_secretion_rate += secretion_rates[i];
	} // for all public goods
	const double waste = chemicals[numchem-1]->value(location);

	const double return_value =
		benefit_constant * total_goods / ( total_goods + benefit_saturation )
		- harm_constant * waste / ( waste + harm_saturation )
		- secretion_cost * total_secretion_rate; 

  return return_value;
}


//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------









/** Fintess function for aging simulations
* rather than benefit, this gives a probability of cell death
* the probability is decreased based on concentration of cooperative factors
* there is technically only one cooperative factor, so this should be a numchem = 1
* class, but since the aging cells are derived from a base class with numchem = 2
* we set numchem =2 for now. Should restucture this later ...
*/

template<int dim>
class AgingDeathProb : public FitnessBase<dim, 2>{
public:
	AgingDeathProb();

	double value(const Point<dim>& location, 
		const std::array<double, 2>& secretion_rates) const override;

	void attach_chemicals(Chemicals::ChemicalInterface<dim>& cfs);

	void set_fitness_constants(double a, double k, double phi);

	void printInfo(std::ostream& out) const;

private:
	const Chemicals::ChemicalInterface<dim>*  cooperative_factors;

	double 		alpha;
	double 		hill_constant_k;
	double 		saturation_constant;

};

// IMPLEMENTATION:
// --------------------------------------------------------------------------
template<int dim>
AgingDeathProb<dim>::AgingDeathProb()
	:
	alpha(0),
	hill_constant_k(1),
	saturation_constant(1)
{}

template<int dim>
double 
AgingDeathProb<dim>::value(const Point<dim>& location, 
		const std::array<double, 2>& /* secretion_rates */) const
{
	return alpha*( 
		std::pow(saturation_constant,hill_constant_k) / 
			( std::pow(saturation_constant,hill_constant_k) +
				std::pow(cooperative_factors->value(location), hill_constant_k) ) 
		);
}

template<int dim>
void 
AgingDeathProb<dim>::attach_chemicals(Chemicals::ChemicalInterface<dim>& cfs)
{
	cooperative_factors = &cfs;
}

template<int dim>
void 
AgingDeathProb<dim>::set_fitness_constants(double a, double k, double phi)
{
	alpha = a;
	hill_constant_k = k;
	saturation_constant = phi;
}

template<int dim>
void 
AgingDeathProb<dim>::printInfo(std::ostream& out) const
{
  out << "\n\n-----------------------------------------------------" << std::endl
    << "\t\t AGING DEATH FUNCTION" << std::endl
    << "-----------------------------------------------------" << std::endl
    << "\t alpha: " << alpha << std::endl
    << "\t hill constant: " << hill_constant_k << std::endl
    << "\t saturation constant: " << saturation_constant << std::endl;
}

}} // close namespace
#endif