#ifndef MICROBE_SIMULATOR_BACTERIA_FITNESS_TYPES_H
#define MICROBE_SIMULATOR_BACTERIA_FITNESS_TYPES_H

#include <deal.II/base/point.h>
using dealii::Point;

#include "../chemicals/chemical_handler.h"
#include "../refactored_chemicals/chemical_handler.h"
#include "../utility/utility.h"

#include <vector>
#include <memory>

namespace MicrobeSimulator{ namespace Bacteria{

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

	out << "\n\n" << Utility::medium_line << std::endl
	    << "\t\tOR FITNESS FUNCTION (for " << numchem << " chemicals)" << std::endl
	    << Utility::medium_line << std::endl
	    << "\t public good benefit: " << benefit_constant << std::endl
	    << "\t waste harm: " << harm_constant << std::endl
	    << "\t public good saturation: " << benefit_saturation << std::endl
	    << "\t waste saturation: " << harm_saturation << std::endl
	    << "\t secretion cost: " << secretion_cost << std::endl;
	out << Utility::medium_line
		<< std::endl << std::endl << std::endl;
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

namespace TestNewFitness{

// ---------------------------------------------------------------------------------
// FITNESS BASE: 
// ---------------------------------------------------------------------------------

/** \brief Fitness function base */
template<int dim>
class FitnessBase{
public:
	FitnessBase(const RefactoredChemicals::ChemicalHandler<dim>& ch);
	virtual ~FitnessBase() {}

	virtual double value(const Point<dim>& location,
					const std::vector<double>& rates) const=0;

	virtual void printInfo(std::ostream& out) const=0;

protected:
	RefactoredChemicals::ChemicalHandler<dim> const * chemicals; 
};

// IMPL
//------------------------------------------------------------

/** \brief Fitness base constructor */
template<int dim>
FitnessBase<dim>::FitnessBase(const RefactoredChemicals::ChemicalHandler<dim>& ch)
	:
	chemicals(&ch)
{}



// ---------------------------------------------------------------------------------
// OR FITNESS:
// ---------------------------------------------------------------------------------

/** \brief Fitness base constructor */
template<int dim>
class OR_Fitness : public FitnessBase<dim>{
public:
	OR_Fitness(const RefactoredChemicals::ChemicalHandler<dim>& ch, 
			const ParameterHandler& prm);

	static void declare_parameters(ParameterHandler& prm);

	double value(const Point<dim>& location,
			const std::vector<double>& rates) const override;

	void printInfo(std::ostream& out) const override;
private:
	double benefit;
	double harm;

	double benefit_saturation;
	double harm_saturation;

	double secretion_cost;
};

// IMPL
// ------------------------------------------------

/** \brief OR type fitness constructor */
template<int dim>
OR_Fitness<dim>::OR_Fitness(const RefactoredChemicals::ChemicalHandler<dim>& ch, 
		const ParameterHandler& prm)
	:
	FitnessBase<dim>(ch)
{
	std::string subsection = "Fitness.OR";

	benefit = prm.get_double(subsection, "Benefit");
	harm = prm.get_double(subsection, "Harm");
	benefit_saturation = prm.get_double(subsection, "Benefit saturation");
	harm_saturation = prm.get_double(subsection, "Harm saturation");
	secretion_cost = prm.get_double(subsection, "Secretion cost");	
}

/** \brief Declare OR fitness type parameters */
template<int dim>
void 
OR_Fitness<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Fitness");
		prm.enter_subsection("OR");
			prm.declare_entry("Benefit","0",Patterns::Double());
			prm.declare_entry("Harm","0",Patterns::Double());
			prm.declare_entry("Benefit saturation","0",Patterns::Double());
			prm.declare_entry("Harm saturation","0",Patterns::Double());
			prm.declare_entry("Secretion cost","0",Patterns::Double());
		prm.leave_subsection();
	prm.leave_subsection();
}
	
/** \brief Return OR type fitness value */
template<int dim>
double 
OR_Fitness<dim>::value(const Point<dim>& location,
		const std::vector<double>& rates) const
{
	const unsigned int numchem = this->chemicals->getNumberChemicals();

	assert(rates.size() == numchem);

	double total_goods = 0;
	double total_secretion_rate = 0;

	for(unsigned int i = 0; i < numchem-1; ++i)
	{
		total_goods += (*this->chemicals)[i].value(location);
		total_secretion_rate += rates[i];
	}
	const double waste = (*this->chemicals)[numchem-1].value(location);

	return benefit*total_goods / ( total_goods + benefit_saturation )
		- harm*waste / ( waste + harm_saturation )
		- secretion_cost * total_secretion_rate; 
}

/** \brief Display OR type fitness info */
template<int dim>
void 
OR_Fitness<dim>::printInfo(std::ostream& out) const
{
	const unsigned int numchem = this->chemicals->getNumberChemicals();

	out << "\n\n" << Utility::medium_line << std::endl
	    << "\t\tOR FITNESS FUNCTION (for " << numchem << " chemicals)" << std::endl
	    << Utility::medium_line << std::endl
	    << "\t public good benefit: " << benefit << std::endl
	    << "\t waste harm: " << harm << std::endl
	    << "\t public good saturation: " << benefit_saturation << std::endl
	    << "\t waste saturation: " << harm_saturation << std::endl
	    << "\t secretion cost: " << secretion_cost << std::endl;
	out << Utility::medium_line
		<< std::endl << std::endl << std::endl;
}

// ---------------------------------------------------------------------------------
// AND FITNESS:
// ---------------------------------------------------------------------------------

/** \brief AND type fitness class */
template<int dim>
class AND_Fitness : public FitnessBase<dim>{
public:
	AND_Fitness(const RefactoredChemicals::ChemicalHandler<dim>& ch, 
			const ParameterHandler& prm);
	
	static void declare_parameters(ParameterHandler& prm);
	
	double value(const Point<dim>& location,
			const std::vector<double>& rates) const override;

	void printInfo(std::ostream& out) const override;
private:
	double benefit;
	double harm;

	double benefit_saturation;
	double harm_saturation;

	double secretion_cost;
};

// IMPL
// ----------------------------------------------------------

/** \brief OR type fitness constructor */
template<int dim>
AND_Fitness<dim>::AND_Fitness(const RefactoredChemicals::ChemicalHandler<dim>& ch, 
		const ParameterHandler& prm)
	:
	FitnessBase<dim>(ch)
{
	std::string subsection = "Fitness.OR";

	benefit = prm.get_double(subsection, "Benefit");
	harm = prm.get_double(subsection, "Harm");
	benefit_saturation = prm.get_double(subsection, "Benefit saturation");
	harm_saturation = prm.get_double(subsection, "Harm saturation");
	secretion_cost = prm.get_double(subsection, "Secretion cost");	
}

/** \brief Declare AND fitness type parameters */
template<int dim>
void 
AND_Fitness<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Fitness");
		prm.enter_subsection("AND");
			prm.declare_entry("Benefit","0",Patterns::Double());
			prm.declare_entry("Harm","0",Patterns::Double());
			prm.declare_entry("Benefit saturation","0",Patterns::Double());
			prm.declare_entry("Harm saturation","0",Patterns::Double());
			prm.declare_entry("Secretion cost","0",Patterns::Double());
		prm.leave_subsection();
	prm.leave_subsection();
}

template<int dim>
double 
AND_Fitness<dim>::value(const Point<dim>& location,
		const std::vector<double>& rates) const
{
	const unsigned int numchem = this->chemicals->getNumberChemicals();

	assert(rates.size() == numchem);

	double total_goods = 1;
	double total_secretion_rate = 0;

	for(unsigned int i = 0; i < numchem-1; ++i)
	{
		total_goods *= (*this->chemicals)[i].value(location); // product of all goods
		total_secretion_rate += rates[i];
	}
	const double waste = (*this->chemicals)[numchem-1].value(location);

	return benefit*total_goods / ( total_goods + benefit_saturation )
		- harm*waste / ( waste + harm_saturation )
		- secretion_cost * total_secretion_rate; 
}

/** \brief Display AND type fitness info */
template<int dim>
void 
AND_Fitness<dim>::printInfo(std::ostream& out) const
{
	const unsigned int numchem = this->chemicals->getNumberChemicals();

	out << "\n\n" << Utility::medium_line << std::endl
	    << "\t\tAND FITNESS FUNCTION (for " << numchem << " chemicals)" << std::endl
	    << Utility::medium_line << std::endl
	    << "\t public good benefit: " << benefit << std::endl
	    << "\t waste harm: " << harm << std::endl
	    << "\t public good saturation: " << benefit_saturation << std::endl
	    << "\t waste saturation: " << harm_saturation << std::endl
	    << "\t secretion cost: " << secretion_cost << std::endl;
	out << Utility::medium_line
		<< std::endl << std::endl << std::endl;
}


// ---------------------------------------------------------------------------------
// FITNESS HANDLER:
// ---------------------------------------------------------------------------------

/** \brief Fitness function class */
template<int dim>
class Fitness_Function{
public:
	Fitness_Function();

	void init(const ParameterHandler& prm,
		const RefactoredChemicals::ChemicalHandler<dim>& ch);

	static void declare_parameters(ParameterHandler& prm);

	double value(const Point<dim>& location,
					const std::vector<double>& rates) const;

	void printInfo(std::ostream& out) const;

private:
	std::shared_ptr<FitnessBase<dim> > 		fitness;
}; 

// IMPL
// ----------------------------------------------------------------

/** \brief Fitness function constructor */
template<int dim>
Fitness_Function<dim>::Fitness_Function()
{}

template<int dim>
void
Fitness_Function<dim>::init(const ParameterHandler& prm,
		const RefactoredChemicals::ChemicalHandler<dim>& ch)
{
	const std::string section = "Fitness";
	std::string fitness_type = prm.get_string(section, "Fitness type");

	if( boost::iequals(fitness_type, "OR") )
		fitness = std::make_shared<OR_Fitness<dim> >(ch, prm);
	else if( boost::iequals(fitness_type, "AND") )
		fitness = std::make_shared<AND_Fitness<dim> >(ch, prm);
	else
		throw std::runtime_error("Invalid fitness type: <" + fitness_type + ">");
}

/** \brief Declare fitness parameters */
template<int dim>
void 
Fitness_Function<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Fitness");
		prm.declare_entry("Fitness type", "OR", Patterns::Selection("OR|AND"));
	prm.leave_subsection();

	OR_Fitness<dim>::declare_parameters(prm);
	AND_Fitness<dim>::declare_parameters(prm);
}

/** \brief Return fitness function value at given location for given secretion rates */
template<int dim>
double 
Fitness_Function<dim>::value(const Point<dim>& location,
				const std::vector<double>& rates) const
{
	return fitness->value(location, rates);
}

/** \brief Display fitness function info */
template<int dim>
void 
Fitness_Function<dim>::printInfo(std::ostream& out) const
{
	fitness->printInfo(out);
}


} // test new fitness








}} // CLOSE NAMESPACE
#endif