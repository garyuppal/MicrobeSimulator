#ifndef MICROBE_SIMULATOR_BACTERIA_HANDLER_H
#define MICROBE_SIMULATOR_BACTERIA_HANDLER_H

#include "./bacterium_types.h"
#include "./bacteria_fitness.h"

#include "../utility/parameter_handler.h"
#include "../geometry/geometry.h"

#include <vector>
#include <memory>

namespace MicrobeSimulator{ 
	/** \brief Bacteria classes and methods
	*/
	namespace Bacteria{

/** \brief Method to get valid locations in geometry domain
*  for seeding bacteria groups
*/
template<int dim>
std::vector<Point<dim> >
get_bacteria_locations(const Geometry<dim>& geometry, unsigned int number_groups,
	double buffer, double left_start_width, double left_start_buffer)
{
	std::cout << "...Finding " << number_groups
		<< " group positions" << std::endl;

	std::vector<Point<dim> > group_locations;
	group_locations.reserve(number_groups);

	Point<dim> temp_point;
	for(unsigned int i = 0; i < number_groups; ++i)
	{
		bool found = false;

		while(!found)
		{
		  for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		  {
		  	double width = geometry.getWidth(dim_itr) - 2.0*buffer;
		  	assert(width>0); // error if buffer too big
		  	// though this should be caught by geometry builder already

		  	// possibly initialize in subdomain:
		  	if( (dim_itr == 0) && (left_start_width > 0) )
		  		width = left_start_width; 

		  	if(left_start_buffer < 0)
		  		left_start_buffer = 0;

		    temp_point[dim_itr] = (width)*((double)rand() / RAND_MAX) 
		      + geometry.getBottomLeftPoint()[dim_itr] + buffer + left_start_buffer;
		  }

		  if( geometry.isInDomain(temp_point, buffer) )
		  {
		    group_locations.emplace_back(temp_point);
		    found = true;
		  }
		} // while not found
	} // for group locations
	std::cout << "...Group positions found." << std::endl;

	return group_locations;
}


/** \brief Bacteria handler class */
template<int dim>
class BacteriaHandler{
public:
	BacteriaHandler();

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm, const Geometry<dim>& geo); 
	void reintro(const ParameterHandler& prm, const Geometry<dim>& geo);

	// MODIFIERS:
	void move(double dt, const Geometry<dim>& geometry, 
		const Velocity::AdvectionHandler<dim>& velocity); 

	void force_mutate(int n_mutate, double ds);
	void mutate(double dt);

	// legacy: (remove)
	void reproduce(double dt, const FitnessBase<dim>& fitness_function); 
	
	// using new chemical handler and fitness handler:
	void reproduce(double dt, const TestNewFitness::Fitness_Function<dim>& fitness_function); 

	// ACCESSORS:
	double getDiffusionConstant() const;
	unsigned int getTotalNumber() const; 

	std::vector<Point<dim> > 		getAllLocations() const;
	std::vector<double> 			getAllRates(unsigned int index) const;
	std::vector<std::vector<double> > getAllRates() const;

	std::vector<double> get_pg_rates();

	bool isAlive() const;

	// OUTPUT:
	void print(std::ostream& out) const;
	void printInfo(std::ostream& out) const;

private:
	std::vector<std::unique_ptr<BacteriumBase<dim> > > bacteria;

	double diffusion_constant;
	double mutation_rate;
	double mutation_strength;
	bool binary_mutation;
	double original_rate; /** @todo need to generalize for multiple chemicals */
	double edge_buffer; 
	double right_open_buffer;

	// recording fallen bacteria:
	std::vector<double> pg_rates;

	void remove_and_capture_fallen_bacteria(
		double right_edge, std::vector<double>& pg_rates);

	void add_bacteria(const ParameterHandler& prm, const Geometry<dim>& geo);
};

// IMPL
// ----------------------------------------------------------------

/** \brief Constructor for BacteriaHandler */
template<int dim>
BacteriaHandler<dim>::BacteriaHandler()
	:
	diffusion_constant(0),
	mutation_rate(0),
	mutation_strength(0),
	binary_mutation(false),
	original_rate(0),
	edge_buffer(0),
	right_open_buffer(0)
{}

/** \brief Declare parameters needed to constuct bacteria */
template<int dim>
void
BacteriaHandler<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Bacteria");
		prm.declare_entry("Number bacteria","100",Patterns::Unsigned());
		prm.declare_entry("Number groups","1",Patterns::Unsigned());
		prm.declare_entry("Initial growth time","0",Patterns::Double());
		prm.declare_entry("Diffusion","0.1",Patterns::Double());
		prm.declare_entry("Edge buffer","0",Patterns::Double());
		prm.declare_entry("Right open buffer","0",Patterns::Double());
		prm.declare_entry("Secretion rate",
							"{100,100}",
							Patterns::List(Patterns::Double()));
		prm.declare_entry("Mutation rate","0",Patterns::Double());
		prm.declare_entry("Mutation strength", "0", Patterns::Double());
		prm.declare_entry("Binary mutation","False",Patterns::Bool());
		prm.declare_entry("Initial number cheaters","0",Patterns::Unsigned());
		prm.declare_entry("Deterministic number mutate","0",Patterns::Unsigned());
		prm.declare_entry("Deterministic mutate time","0",Patterns::Double());
		prm.declare_entry("Deterministic mutation strength","10000",Patterns::Double());
		prm.declare_entry("Initial locations",
							"{{}}",
							Patterns::List(Patterns::List(Patterns::Double())));
		prm.declare_entry("Reintroducing","False",Patterns::Bool());
		prm.declare_entry("Reintroduction period","0",Patterns::Double());
		// left length, for reintoduction
		prm.declare_entry("Left start width","-1",Patterns::Double());
		prm.declare_entry("Left start buffer","0",Patterns::Double());
	prm.leave_subsection();
}

/** \brief Initialization of bacteria from system parameters and bounding geometry */
/** Note geometry should be constructed before initializing bacteria */
template<int dim>
void 
BacteriaHandler<dim>::init(const ParameterHandler& prm, const Geometry<dim>& geo)
{
	const std::string section = "Bacteria";

	// assign constants:
	mutation_rate = prm.get_double(section, "Mutation rate");
	mutation_strength = prm.get_double(section, "Mutation strength");
	binary_mutation = prm.get_bool(section, "Binary mutation");
	/** @todo generalize this to multiple possible chemicals */
	original_rate = prm.get_double_vector(section, "Secretion rate")[0];
	edge_buffer = prm.get_double(section,"Edge buffer");
	right_open_buffer = prm.get_double(section, "Right open buffer");

	// add bacteria:
	const unsigned int n_bact = prm.get_unsigned(section, "Number bacteria");
	bacteria.clear();
	bacteria.reserve(n_bact);

	add_bacteria(prm,geo);
} 

/** \brief Method to introduce new bacteria groups during a simulation */
template<int dim>
void
BacteriaHandler<dim>::reintro(const ParameterHandler& prm, const Geometry<dim>& geo)
{
	add_bacteria(prm,geo);
}

/** \brief Method to support adding bacteria to handler class. Called by init() and reintro() methods. */
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

	double edge_buffer = prm.get_double("Bacteria","Edge buffer");

	if(initial_locations.empty()) // if not given from parameters
	{
		if(number_groups == 0)
			number_groups = n_bact;
		
		const double left_start_width = prm.get_double(section, "Left start width");
		const double left_start_buffer = prm.get_double(section, "Left start buffer");

		initial_locations = get_bacteria_locations(geo, number_groups, 
			edge_buffer, left_start_width, left_start_buffer);
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
BacteriaHandler<dim>::move(double dt, const Geometry<dim>& geometry, 
	const Velocity::AdvectionHandler<dim>& velocity)
{
	for(unsigned int i = 0; i < bacteria.size(); ++i)
		bacteria[i]->randomStep(dt, diffusion_constant, 
			geometry, velocity, edge_buffer);

	// if boundary is open, remove fallen bacteria:
	if(geometry.getBoundaryConditions()[0] == BoundaryCondition::OPEN)
		remove_and_capture_fallen_bacteria(geometry.getTopRightPoint()[0] - right_open_buffer,
											pg_rates);	

}

template<int dim>
void
BacteriaHandler<dim>::force_mutate(int n_mutate, double ds)
{
	unsigned int mutated = 0;
	unsigned int n_bact = bacteria.size();
	n_mutate = (n_mutate < n_bact)? n_mutate : n_bact;

	auto it = bacteria.begin();
	do{
		double current_sec = (*it)->getSecretionRate(0);
		if( current_sec > 0)
		{
			double set_sec = current_sec -  ds;
			set_sec = (set_sec < 0)? 0. : set_sec; 
			(*it)->setSecretionRate(0, set_sec);
			++mutated;
		}
		++it;
	}while( (mutated < n_mutate) && (it !=bacteria.end() ) );
}

template<int dim>
void 
BacteriaHandler<dim>::mutate(double dt)
{
	// loop through bacteria:
	for(auto it = bacteria.begin(); it != bacteria.end(); ++it)
	{
		double prob = Utility::getRand(); 
		if(prob < (dt*mutation_rate) )
		{
			double sec = original_rate;
			if(binary_mutation)
			{
				double current_sec = (*it)->getSecretionRate(0);	
				if(current_sec > 0)
					sec = 0;
			}
			else
			{
				sec = (*it)->getSecretionRate(0) 
					+ mutation_strength*(2.0*Utility::getRand()-1.0);
				if(sec < 0)
					sec = 0;
			}
			(*it)->setSecretionRate(0, sec);
		} // if mutating
	} // for all bacteria
}

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

// LEGACY.... REMOVE
/** \brief Reproduce bacteria according to supplied fitness function and time step */
/** @todo Can probably speed up */
template<int dim>
void 
BacteriaHandler<dim>::reproduce(
	double dt, const FitnessBase<dim>& fitness_function)
{
// simple implementaiton to get working, speed up later:
// instead of recloning, and for future improvement, can we move the pointer instead? (there should be  a move function for unique_ptr)
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
} // reproduce()

/** \brief Reproduce bacteria according to supplied fitness function and time step */
// USING NEW CHEM AND FITNESS:
/** @todo Can probably speed up */
template<int dim>
void 
BacteriaHandler<dim>::reproduce(
	double dt, const TestNewFitness::Fitness_Function<dim>& fitness_function)
{
// simple implementaiton to get working, speed up later:
// instead of recloning, and for future improvement, can we move the pointer instead? (there should be  a move function for unique_ptr)
	std::vector<std::unique_ptr<BacteriumBase<dim> > > offspring;

	for(auto it = bacteria.begin(); it != bacteria.end(); )
	{
		const double fit = fitness_function.value( (*it)->getLocation(),
												(*it)->getSecretionRates()
												)*dt; 
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
} // reproduce()

// ACCESSORS:
// -------------------------------------------------------

/** \brief Returns diffusion constant for bacteria */
template<int dim>
double 
BacteriaHandler<dim>::getDiffusionConstant() const
{
	return diffusion_constant;
}

/** \brief Returns total number of current bacteria */
template<int dim>
unsigned int 
BacteriaHandler<dim>::getTotalNumber() const
{
	return bacteria.size();
}

/** \brief Returns all locations of bacteria */
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

/** \brief Return all secretion rates for ith chemical */
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

/** \brief Returns vector of vector of doubles for all secretion rates
* of all chemicals */
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

	all_rates.reserve(num_chem);
	for(unsigned int c = 0; c < num_chem; ++c)
		all_rates.emplace_back(getAllRates(c));

	return all_rates;
}

/** \brief Return public good secretion rates of captured bacteria */
template<int dim>
std::vector<double>
BacteriaHandler<dim>::get_pg_rates()
{
	return pg_rates;
}

/** \brief Returns true is bacteria not empty */
template<int dim>
bool 
BacteriaHandler<dim>::isAlive() const
{
	return !(bacteria.empty());
}

/** \brief Prints out bacteria info for each bacteria to given ostream object */
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

/** \brief Prints out general info for bacteria to provided ostream object */
template<int dim>
void 
BacteriaHandler<dim>::printInfo(std::ostream& out) const
{
	std::string mstren;
	if(binary_mutation)
		mstren = "Binary";
	else
		mstren = std::to_string(mutation_strength);

	out << "\n\n" << Utility::medium_line << std::endl 
		<< "\t\t BACTERIA INFO:" << std::endl
		<< Utility::medium_line << std::endl
		<< "\t Diffusion constant: " << diffusion_constant << std::endl
		<< "\t Mutation rate: " << mutation_rate << std::endl
		<< "\t Mutation strength: " << mstren << std::endl
		<< "\t Edge buffer: " << edge_buffer << std::endl
		<< "\t Right open buffer: " << right_open_buffer << std::endl
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