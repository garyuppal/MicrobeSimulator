#ifndef MICROBE_SIMULATOR_BACTERIA_HANDLER_2_H
#define MICROBE_SIMULATOR_BACTERIA_HANDLER_2_H

#include "./bacterium_types.h"
#include "./bacteria_fitness.h"

#include "../utility/parameter_handler.h"
#include "../simulators/simulation_tools.h" // move this method to this file...
#include "../geometry/geometry.h"

#include <vector>
#include <memory>

namespace MicrobeSimulator{ namespace BacteriaNew{

/** \brief Bacteria handler class
*/
template<int dim>
class BacteriaHandler{
public:
	BacteriaHandler();

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm, const Geometry<dim>& geo);

	void reintro(const ParameterHandler& prm, const Geometry<dim>& geo);

	// temporary initialization for testing:
	void init_reg(double db, unsigned int n_bact, double pg_rate, double waste_rate);

	// MODIFIERS:
	void randomWalk(double dt, const Geometry<dim>& geometry,
			const Velocity::AdvectionHandler<dim>& velocity, double buffer=0);
	void randomWalk(double dt, const Geometry<dim>& geometry,
			const Velocity::AdvectionHandler<dim>& velocity,
			std::vector<double>& pg_rates, double buffer=0);

	void reproduce(double dt, const FitnessBase<dim>& fitness_function); // need to change this
	// void mutate(double dt); // general mutation method ... 

	// simple version for now:
	void mutate_simple(double dt, double mutation_rate, double ds);

	// ACCESSORS:
	double getDiffusionConstant() const;
	unsigned int getTotalNumber() const; 

	std::vector<Point<dim> > 		getAllLocations() const;
	std::vector<double> 			getAllRates(unsigned int index) const;
	std::vector<std::vector<double> > getAllRates() const;

	bool isAlive() const;

	// OUTPUT:
	void print(std::ostream& out) const;
	void printInfo(std::ostream& out) const;

private:
	std::vector<std::unique_ptr<BacteriumBase<dim> > > bacteria;

	double diffusion_constant;

	void remove_fallen_bacteria(double right_edge);
	void remove_and_capture_fallen_bacteria(
		double right_edge, std::vector<double>& pg_rates);

	void add_bacteria(const ParameterHandler& prm, const Geometry<dim>& geo);
};

// IMPL
// ----------------------------------------------------------------
template<int dim>
BacteriaHandler<dim>::BacteriaHandler()
{}

template<int dim>
void
BacteriaHandler<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Bacteria");
		prm.declare_entry("Number bacteria","100",Patterns::Unsigned());
		prm.declare_entry("Number groups","1",Patterns::Unsigned());
		prm.declare_entry("Diffusion","0.1",Patterns::Double());
		prm.declare_entry("Edge buffer","0",Patterns::Double());
		prm.declare_entry("Secretion rate",
							"{100,100}",
							Patterns::List(Patterns::Double()));
		prm.declare_entry("Mutation rate","0",Patterns::Double());
		prm.declare_entry("Mutation strength", "0", Patterns::Double());
		prm.declare_entry("Initial number cheaters","0",Patterns::Unsigned());
		prm.declare_entry("Deterministic number mutate","0",Patterns::Unsigned());
		prm.declare_entry("Deterministic mutate time","0",Patterns::Double());
		prm.declare_entry("Initial locations",
							"{{}}",
							Patterns::List(Patterns::List(Patterns::Double())));
		prm.declare_entry("Reintroducing","False",Patterns::Bool());
		prm.declare_entry("Reintroduction period","0",Patterns::Double());
		// left length, for reintoduction
		prm.declare_entry("Left subdomain length","-1",Patterns::Double());
	prm.leave_subsection();
}

template<int dim>
void 
BacteriaHandler<dim>::init(const ParameterHandler& prm, const Geometry<dim>& geo)
{
	const unsigned int n_bact = prm.get_unsigned("Bacteria", "Number bacteria");
	bacteria.clear();
	bacteria.reserve(n_bact);

	add_bacteria(prm,geo);
} 

template<int dim>
void
BacteriaHandler<dim>::reintro(const ParameterHandler& prm, const Geometry<dim>& geo)
{
	add_bacteria(prm,geo);
}

template<int dim>
void
BacteriaHandler<dim>::add_bacteria(const ParameterHandler& prm, const Geometry<dim>& geo)
{
	const std::string section = "Bacteria";

	diffusion_constant = prm.get_double(section,"Diffusion");
	const unsigned int n_bact = prm.get_unsigned(section, "Number bacteria");

	std::vector<double> rates = prm.get_double_vector(section, "Secretion rate");
	unsigned int number_groups = prm.get_unsigned(section, "Number groups");

	std::vector<Point<2> > initial_locations = prm.get_point_list(section, "Initial locations");

	if(initial_locations.empty()) // if not given from parameters
	{
		if(number_groups == 0)
			number_groups = n_bact;
		const double left_length = prm.get_double(section, "Left subdomain length");
		initial_locations = SimulationTools::get_bacteria_locations(geo, number_groups, left_length);
	}

	for(unsigned int i = 0; i < n_bact; ++i)
	{
		unsigned int group_index = i % number_groups;	
		Point<dim> location = initial_locations[group_index];

		bacteria.emplace_back(new BacteriumBase<dim>(location, rates));
	}	
}

template<int dim>
void
BacteriaHandler<dim>::init_reg(double db, 
	unsigned int n_bact, double pg_rate, double waste_rate)
{
	diffusion_constant = db;
	bacteria.clear();
	bacteria.reserve(n_bact);

	Point<dim> location = (dim == 2) ? Point<dim>(0,0)
									: Point<dim>(0,0,0);

	std::vector<double> rates = {pg_rate,waste_rate};

	for(unsigned int i = 0; i < n_bact; ++i)
		bacteria.emplace_back(new BacteriumBase<dim>(location,rates) );
}

template<int dim>
void 
BacteriaHandler<dim>::randomWalk(double dt, const Geometry<dim>& geometry,
			const Velocity::AdvectionHandler<dim>& velocity,
			double buffer)
{
	// std::cout << "doing random walk" << std::endl;
	for(unsigned int i = 0; i < bacteria.size(); ++i)
		bacteria[i]->randomStep(dt, diffusion_constant, geometry, velocity, buffer);

	// if boundary is open, remove fallen bacteria:
	if(geometry.getBoundaryConditions()[0] == BoundaryCondition::OPEN)
		remove_fallen_bacteria(geometry.getTopRightPoint()[0]);	
}

template<int dim>
void 
BacteriaHandler<dim>::randomWalk(double dt, const Geometry<dim>& geometry,
		const Velocity::AdvectionHandler<dim>& velocity,
		std::vector<double>& pg_rates,
		double buffer)
{
	// std::cout << "doing random walk" << std::endl;
	for(unsigned int i = 0; i < bacteria.size(); ++i)
		bacteria[i]->randomStep(dt, diffusion_constant, geometry, velocity, buffer);

	// if boundary is open, remove fallen bacteria:
	if(geometry.getBoundaryConditions()[0] == BoundaryCondition::OPEN)
		remove_and_capture_fallen_bacteria(geometry.getTopRightPoint()[0],
											pg_rates);	
}

/** \brief Remove bacteria that move past right boundary*/
template<int dim>
void 
BacteriaHandler<dim>::remove_fallen_bacteria(double right_edge)
{
	const double tolerance = 1e-4;
	for(auto it = bacteria.begin(); it != bacteria.end(); )
	{
		if( (*it)->getLocation()[0] > (right_edge - tolerance) ) 
			it = bacteria.erase(it);
		else
			++it;
	} // for all bacteria
} // remove_fallen_bacteria()

/** \brief Remove and record public good secretion rates of fallen bacteria */
template<int dim>
void 
BacteriaHandler<dim>::remove_and_capture_fallen_bacteria(
	double right_edge, std::vector<double>& pg_rates)
{
	const double tolerance = 1e-4;
	for(auto it = bacteria.begin(); it != bacteria.end(); )
	{
		if( (*it)->getLocation()[0] > (right_edge - tolerance) ) 
		{
			pg_rates.emplace_back((*it)->getSecretionRate(0)); 
				/** @todo generalize to multiple public goods */
			it = bacteria.erase(it);
		}
		else
		{
			++it;
		}
	} // for all bacteria
} // remove_and_capture_fallen_bacteria()

// simple implementaiton to get working, speed up later:
// *** instead of recloning, and for future improvement,
// can we move the pointer instead? (there should be  a move function for unique_ptr)
template<int dim>
void 
BacteriaHandler<dim>::reproduce(
	double dt, const FitnessBase<dim>& fitness_function)
{
	std::vector<std::unique_ptr<BacteriumBase<dim> > > offspring;

	for(auto it = bacteria.begin(); it != bacteria.end(); )
	{
		const double fit = (*it)->getFitness(fitness_function)*dt; 
		if(fit < 0)
		{
			double prob = Utility::getRand();
			if(prob < (-fit) )
			{
				it = bacteria.erase(it);
			}
			else
				++it;
		} // if possibly dying
		else
		{
			double prob = Utility::getRand();
			if(prob < fit)
			{
				offspring.emplace_back( (*it)->clone() );
			}
			++it;
		} // else possibly reproduce
	} // for all bacteria

	// add offspring to end: 
	for(unsigned int i = 0; i < offspring.size(); ++i)
		bacteria.emplace_back( offspring[i]->clone() );
	// bacteria.insert(bacteria.end(), offspring.begin(), offspring.end());
} // reproduce()

template<int dim>
void 
BacteriaHandler<dim>::mutate_simple(double dt, double mutation_rate, double ds)
{
	// loop through bacteria:
	for(auto it = bacteria.begin(); it != bacteria.end(); ++it)
	{
		double prob = dt*Utility::getRand();
		if(prob < mutation_rate)
		{
			double sec = (*it)->getSecretionRate(0) + ds*(2.0*Utility::getRand()-1.0);
			if(sec < 0)
				sec = 0;
			(*it)->setSecretionRate(0, sec);
		} // if mutating
	} // for all bacteria
} // mutate_simple()

// ACCESSORS:
template<int dim>
double 
BacteriaHandler<dim>::getDiffusionConstant() const
{
	return diffusion_constant;
}

template<int dim>
unsigned int 
BacteriaHandler<dim>::getTotalNumber() const
{
	return bacteria.size();
}

template<int dim>
std::vector<Point<dim> > 		
BacteriaHandler<dim>::getAllLocations() const
{
	std::vector<Point<dim> > all_loc;
	all_loc.reserve(bacteria.size());

	for(unsigned int i = 0; i < bacteria.size(); ++i)
		all_loc.emplace_back(bacteria[i]->getLocation());

	return all_loc;
}

template<int dim>
std::vector<double> 			
BacteriaHandler<dim>::getAllRates(unsigned int index) const
{
	std::vector<double> all_rates;

	all_rates.reserve(bacteria.size());

	for(unsigned int i = 0; i < bacteria.size(); ++i)
		all_rates.emplace_back(bacteria[i]->getSecretionRate(index));

	return all_rates;
}

template<int dim>
std::vector<std::vector<double> > 
BacteriaHandler<dim>::getAllRates() const
{
	std::vector<std::vector<double> > all_rates;

	if(bacteria.empty())
	{
		all_rates.clear();
		return all_rates;
	}

	const unsigned int num_chem = bacteria[0]->getSecretionRates().size();
	// std::cout << "\n\n\n\n!!!!number chemicals = " << num_chem << "\n\n\n" << std::endl;

	all_rates.reserve(num_chem);
	for(unsigned int c = 0; c < num_chem; ++c)
		all_rates.emplace_back(getAllRates(c));

	// all_rates.reserve(num_chem);

	// for(unsigned int c = 0; c < num_chem; ++c)
	// {
	// 	all_rates[c].reserve(bacteria.size());
	// 	for(unsigned int i = 0; i < bacteria.size(); ++i)
	// 		all_rates[c].emplace_back(bacteria[i]->getSecretionRate(c)); 
	// }

	return all_rates;
}


// template<int dim>
// std::vector<std::vector<double> > 
// BacteriaHandler<dim>::getAllRates() const
// {

// }

template<int dim>
bool 
BacteriaHandler<dim>::isAlive() const
{
	return !(bacteria.empty());
}


template<int dim>
void
BacteriaHandler<dim>::print(std::ostream& out) const
{
	for(unsigned int i = 0; i < bacteria.size(); ++i)
	{
		bacteria[i]->print(out);
		out << std::endl;
	}
}

template<int dim>
void 
BacteriaHandler<dim>::printInfo(std::ostream& out) const
{
	out << "\n\n" << Utility::medium_line << std::endl 
		<< "\t\t BACTERIA INFO:" << std::endl
		<< Utility::medium_line << std::endl
		<< "\t Diffusion constant: " << diffusion_constant << std::endl
		<< "\t Number bacteria: " << bacteria.size() << std::endl
		<< "\t First Bacterium: " << bacteria[0]->getLocation() << std::endl
		<< "\t First secretion rates:";

	std::vector<double> srates = bacteria[0]->getSecretionRates();
	for(unsigned int i = 0; i < srates.size(); ++i)
		out << srates[i] << " ";
	out << std::endl << Utility::medium_line << std::endl
		<< std::endl << std::endl;

}

}} // CLOSE NAMESPACE
#endif