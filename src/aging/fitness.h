#pragma once

#include <deal.II/base/point.h>
using dealii::Point;

#include "../utility/utility.h"
#include "../utility/parameter_handler.h"
#include "./chemicals.h"

namespace MicrobeSimulator{ namespace Aging{

template<class T>
void
print_list(std::ostream& out, const std::vector<T>& list)
{
	if(list.empty())
	{
		out << std::endl;
	}
	else if( list.size() == 1)
	{
		out << list[0] << std::endl;
	}
	else
	{
		const unsigned int n = list.size()-1;

		for(unsigned int i = 0; i < n; ++i)
	    	out << list[i] << ", ";
	    out << list[n] << std::endl;
    }
}
    
 
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





// --------------------------------------------------------------------------------
// MULTI CHEM AGING FITNESS:
// --------------------------------------------------------------------------------
template<int dim>
class Multi_Aging : public FitnessBase<dim>{
public:
	Multi_Aging(const ParameterHandler& prm, const Aging::Chemicals<dim>& ch);
	
	static void declare_parameters(ParameterHandler& prm);
	
	double value(const Point<dim>& location) const override;

	void printInfo(std::ostream& out) const override;
private:
	unsigned int n_pos, n_neg;
	std::vector<double> pos_coefs;
	std::vector<double> neg_coefs;
	std::vector<double> pos_sats;
	std::vector<double> neg_sats;
	std::vector<double> pos_hills;
	std::vector<double> neg_hills;
	double constant;

	void check_dimensions() const;
};

// IMPL
// ----------------------------------------------------------

/** \brief Multi aging fitness type constructor */
template<int dim>
Multi_Aging<dim>::Multi_Aging(const ParameterHandler& prm,
	 const Aging::Chemicals<dim>& ch)
	:
	FitnessBase<dim>(prm, ch)
{
	std::string section = "Fitness.Multi Aging";


	const bool mixing = prm.get_bool(section,"Mixing");
	if(mixing)
	{
		n_pos = 1;
		n_neg = 1;

		const double alpha = prm.get_double(section,"Alpha");
		const double gamma = prm.get_double(section,"Gamma");

		std::vector<double> hill_constants = prm.get_double_vector(section, "Hill constants");

		pos_coefs.emplace_back(alpha*gamma);
		neg_coefs.emplace_back(alpha*(1. - gamma));
		pos_sats.emplace_back(prm.get_double(section,"Saturation constant"));
		neg_sats.emplace_back(prm.get_double(section,"Saturation constant"));
		// pos_hills.emplace_back(prm.get_double(section,"Hill constant"));
		// neg_hills.emplace_back(prm.get_double(section,"Hill constant"));
		pos_hills.emplace_back(hill_constants[0]);
		neg_hills.emplace_back(hill_constants[1]);
		constant = prm.get_double(section, "Constant");

	}
	else
	{
		n_pos = prm.get_unsigned(section,"Number CFs");
		n_neg = prm.get_unsigned(section, "Number toxins");

		pos_coefs = prm.get_double_vector(section,"CF coefs");
		neg_coefs = prm.get_double_vector(section, "Toxin coefs");
		pos_sats = prm.get_double_vector(section, "Benefit saturation");
		neg_sats = prm.get_double_vector(section, "Harm saturation");
		pos_hills = prm.get_double_vector(section, "Benefit hill");
		neg_hills = prm.get_double_vector(section, "Harm hill");
		constant = prm.get_double(section, "Constant");
	}

	check_dimensions();
}

template<int dim>
void
Multi_Aging<dim>::check_dimensions() const
{
	assert( pos_coefs.size() == n_pos );
	assert( pos_sats.size() == n_pos );
	assert( pos_hills.size() == n_pos );

	assert( neg_coefs.size() == n_neg );
	assert( neg_sats.size() == n_neg );
	assert( neg_hills.size() == n_neg );
}

/** \brief Declare AND fitness type parameters */
template<int dim>
void 
Multi_Aging<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Fitness");
		prm.enter_subsection("Multi Aging");
			prm.declare_entry("Number CFs", "0", Patterns::Unsigned());
			prm.declare_entry("Number toxins", "0", Patterns::Unsigned());
			prm.declare_entry("CF coefs","{0}", Patterns::List(Patterns::Double()));
			prm.declare_entry("Toxin coefs","{0}", Patterns::List(Patterns::Double()));
			prm.declare_entry("Benefit saturation","{0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Harm saturation","{0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Benefit hill","{0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Harm hill","{0}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Constant","0",Patterns::Double());

			prm.declare_entry("Mixing","False",Patterns::Bool());
			prm.declare_entry("Alpha","0",Patterns::Double());
			prm.declare_entry("Saturation constant","1",Patterns::Double());
			prm.declare_entry("Hill constants","{1,1}",Patterns::List(Patterns::Double()));
			prm.declare_entry("Gamma","0",Patterns::Double());

		prm.leave_subsection();
	prm.leave_subsection();
}

template<int dim>
double 
Multi_Aging<dim>::value(const Point<dim>& location) const
{
	double fit = 0.;

	// chemicals ordered so benefits first, negs later
	for(unsigned int i = 0; i < n_pos; ++i)
	{
		const double chemk = std::pow( (*this->chemicals)[i].value(location), pos_hills[i]);

		fit += pos_coefs[i]*(
			 	std::pow(pos_sats[i],pos_hills[i])/
				( chemk + std::pow(pos_sats[i],pos_hills[i]) ) 
			);
	}

	for(unsigned int i = 0; i < n_neg; ++i)
	{
		const double chemk = std::pow( (*this->chemicals)[i + n_pos].value(location), neg_hills[i]);

		fit += neg_coefs[i]*(
			 	 chemk /
				(  chemk + std::pow(neg_sats[i],neg_hills[i]) ) 
			);
	}

	return -fit + constant;
}

template<int dim>
void 
Multi_Aging<dim>::printInfo(std::ostream& out) const
{
	out << "\n\n" << Utility::medium_line << std::endl
	    << "\t\tMULTI AGING FITNESS (for " << n_pos << " CFS and " 
	    		<< n_neg << " toxins)" << std::endl
	    << Utility::medium_line << std::endl;
	    
    out << "\t public good benefits: ";
    print_list(out, pos_coefs);

    out << "\t public good saturation constants: ";
    print_list(out, pos_sats);

    out << "\t public good hill constants: ";
    print_list(out, pos_hills);

    out << "\t toxin coefficients: ";
    print_list(out, neg_coefs);

    out << "\t toxin saturation constants: ";
    print_list(out, neg_sats);

    out << "\t toxin hill constants: ";
    print_list(out, neg_hills);

	out << "\t constant: " << constant << std::endl;
	out << Utility::medium_line
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
	Multi_Aging<dim>::declare_parameters(prm);
}

template<int dim>
void 
Fitness<dim>::init(const ParameterHandler& prm, const Aging::Chemicals<dim>& chem)
{
	const std::string section = "Fitness";
	const std::string fitness_type = prm.get_string(section, "Fitness type");

	if ( boost::iequals(fitness_type, "CF") )
		fit = std::make_shared<CF_Fitness<dim> >(prm, chem);
	else if( boost::iequals(fitness_type, "Multi Aging") )
		fit = std::make_shared<Multi_Aging<dim> >(prm, chem);
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