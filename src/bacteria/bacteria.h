#pragma once

#include "./bacteria_fitness.h"
#include "../geometry/geometry.h"
#include "../advection/advection_handler.h"

#include "../utility/parameter_handler.h" 
#include "../utility/utility.h"

#include <vector>
#include <memory>

namespace MicrobeSimulator{ 
	/** \brief Bacteria classes and methods
	*/
	namespace BacteriaTools{
		// namespace MicrobeSimulator::Bacteria::TestNewFitness = tf; 
/** \brief Method to get valid locations in geometry domain
*  for seeding bacteria groups
*/
template<int dim>
std::vector<Point<dim> >
get_bacteria_locations(const Geometry<dim>& geometry, unsigned int number_groups,
	double buffer, double left_start_width, double left_start_buffer, double y_start_buffer)
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

		  	if( (dim_itr == 1) && (y_start_buffer > 0) )
		  		width -= 2*y_start_buffer;

		  	if(left_start_buffer < 0)
		  		left_start_buffer = 0;

		    temp_point[dim_itr] = (width)*((double)rand() / RAND_MAX) 
		      + geometry.getBottomLeftPoint()[dim_itr] + buffer;
		  } // set temp point
	  	temp_point[0] = temp_point[0] + left_start_buffer; // buffer for left side (x axis)
	  	temp_point[1] = temp_point[1] + y_start_buffer;

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



// ------------------------------------------------------------------------------
// BACTERIUM BASE CLASS
// ------------------------------------------------------------------------------
/** \brief Bacterium base class */
template<int dim>
class BacteriumBase{
public:
	BacteriumBase(const ParameterHandler& prm, const Point<dim>& p);
	// BacteriumBase(const Point<dim>& p);
	// BacteriumBase(const std::vector<double>& rates);
	// BacteriumBase(const Point<dim>& p, 
	// 				const std::vector<double>& rates);

	virtual ~BacteriumBase() {}

	static void declare_parameters(ParameterHandler& prm);

	virtual std::unique_ptr<BacteriumBase<dim> > clone() const;

	// possibly motify for run and tumble motion:
	virtual void move(double dt, 
						double diff, 
						const Geometry<dim>& geo,
						const Velocity::AdvectionHandler<dim>& velocity,
						double buffer=0); 

	// modify for intermittent cheaters:
	virtual double getFitness(const MicrobeSimulator::Bacteria::TestNewFitness::Fitness_Function<dim>& fitness_function) const;
	virtual std::vector<double> getSecretionRates() const;
	virtual double getSecretionRate(unsigned int i) const;

	// mainly used by intermittent cheaters:
	virtual void update_state(double dt, 
		const RefactoredChemicals::ChemicalHandler<dim>& chemicals); 
	
	// MUTATORS:
	void setLocation(const Point<dim>& p);
	void setSecretionRates(const std::vector<double> rates);
	void setSecretionRate(unsigned int index, double value);

	// ACCESSORS:
	Point<dim> getLocation() const;
	unsigned int getNumberChemicals() const;

	virtual void print(std::ostream& out) const;

protected:
	Point<dim> location;
	std::vector<double> secretion_rates;
};

// IMPL
// ------------------------------------------------------------------------------

template<int dim>
BacteriumBase<dim>::BacteriumBase(const ParameterHandler& prm, const Point<dim>& p)
{
	const std::string section = "Bacteria.Base";

	secretion_rates = prm.get_double_vector(section,"Secretion rate");
	location = p;
}

template<int dim>
void 
BacteriumBase<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Bacteria");
		prm.enter_subsection("Base");
			prm.declare_entry("Secretion rate",
					"{100,100}",
					Patterns::List(Patterns::Double()));
			prm.declare_entry("Number bacteria","0",Patterns::Unsigned());
		prm.leave_subsection();
	prm.leave_subsection();
}

/** \brief random walk step */
/** buffer set to 0 as default. Non-zero buffer corresponds to an extra
*	reflection length against solid boundaries 
*/
template<int dim>
void 
BacteriumBase<dim>::move(double time_step, 
						double diffusion_constant,
						const Geometry<dim>& geometry,
						const Velocity::AdvectionHandler<dim>& velocity,
						double buffer)
{
	Point<dim> old_location(location);

	const double theta = 2*dealii::numbers::PI*Utility::getRand();
	const double phi = dealii::numbers::PI*Utility::getRand();

	Point<dim> randomPoint = (dim == 2) ? Point<dim>(std::cos(theta),std::sin(theta))
	                        : Point<dim>(std::cos(phi), std::sin(phi)*std::cos(theta),
	                            std::sin(phi)*std::sin(theta));

	location += std::sqrt(2*dim*time_step*diffusion_constant)*randomPoint
	    + time_step*velocity.value(location);

	geometry.checkBoundaries(old_location, location, buffer); 
}

template<int dim>
double 
BacteriumBase<dim>::getFitness(const Bacteria::TestNewFitness::Fitness_Function<dim>& fitness_function) const
{
	return fitness_function.value(location, secretion_rates);
}

template<int dim>	
std::unique_ptr<BacteriumBase<dim> > 
BacteriumBase<dim>::clone() const
{
	return std::unique_ptr<BacteriumBase<dim> >(
			new BacteriumBase(*this) );
}

template<int dim>
void 
BacteriumBase<dim>::update_state(double /* dt */, 
	const RefactoredChemicals::ChemicalHandler<dim>& /* chemicals*/)
{
	return;
}
	
// MUTATORS:
template<int dim>
void 
BacteriumBase<dim>::setLocation(const Point<dim>& p)
{
	location = p;
}

template<int dim>
void 
BacteriumBase<dim>::setSecretionRates(const std::vector<double> rates)
{
	assert(secretion_rates.size() == rates.size());
	secretion_rates = rates;
}

template<int dim>
void 
BacteriumBase<dim>::setSecretionRate(unsigned int index, double value)
{
	assert(index < secretion_rates.size());
	secretion_rates[index] = value;
}

// ACCESSORS:
template<int dim>
Point<dim> 
BacteriumBase<dim>::getLocation() const
{
	return location;
}

template<int dim>
std::vector<double> 
BacteriumBase<dim>::getSecretionRates() const
{
	return secretion_rates;
}

template<int dim>
double 
BacteriumBase<dim>::getSecretionRate(unsigned int i) const
{
	assert(i < secretion_rates.size());
	return secretion_rates[i];
}

template<int dim>
unsigned int
BacteriumBase<dim>::getNumberChemicals() const
{
	return secretion_rates.size();
}

template<int dim>
void 
BacteriumBase<dim>::print(std::ostream& out) const
{
	out << location << " ";
	for(unsigned int i = 0; i < secretion_rates.size()-1; ++i)
		out << secretion_rates[i] << " ";
	out << secretion_rates[secretion_rates.size()-1];
}





// ------------------------------------------------------------------------------
// PERIODICALLY INTERMITTENTLY CHEATING BACTERIA
// ------------------------------------------------------------------------------
/** \brief Intermittently cheating bacteria */
template<int dim>
class PICBacterium : public BacteriumBase<dim>{
public:
	PICBacterium(const ParameterHandler& prm, const Point<dim>& p);

	static void declare_parameters(ParameterHandler& prm);

	std::unique_ptr<BacteriumBase<dim> > clone() const override;

	// modify for intermittent cheaters:
	double getFitness(const Bacteria::TestNewFitness::Fitness_Function<dim>& fitness_function) const override;
	std::vector<double> getSecretionRates() const override;
	double getSecretionRate(unsigned int i) const override;

	// mainly used by intermittent cheaters:
	void update_state(double dt, 
		const RefactoredChemicals::ChemicalHandler<dim>& chemicals) override; 
	
	void print(std::ostream& out) const override;

private:
	// inherited:
	// Point<dim> location;
	// std::vector<double> secretion_rates;	

	double time; // internal clock
	double delay;
	double on_period;
	double off_period;

	double get_current_public_good() const;
};

// IMPL
// ------------------------------------------------------------------------------
template<int dim>
PICBacterium<dim>::PICBacterium(const ParameterHandler& prm, const Point<dim>& p)
	:
	BacteriumBase<dim>(prm,p),
	time(0)
{
	const std::string section = "Bacteria.PIC";

	delay = prm.get_double(section,"Delay");
	on_period = prm.get_double(section,"On period");
	off_period = prm.get_double(section,"Off period");

	this->secretion_rates = prm.get_double_vector(section,"Secretion rate");
	this->location = p;
}

template<int dim>
void 
PICBacterium<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Bacteria");
		prm.enter_subsection("PIC");
			prm.declare_entry("Number bacteria","0",Patterns::Unsigned());
			prm.declare_entry("Secretion rate",		
							"{100,100}",
							Patterns::List(Patterns::Double())); 
			prm.declare_entry("Delay","0",Patterns::Double());
			prm.declare_entry("On period","0",Patterns::Double());
			prm.declare_entry("Off period","0",Patterns::Double());
		prm.leave_subsection();
	prm.leave_subsection();
}

template<int dim>
std::unique_ptr<BacteriumBase<dim> > 
PICBacterium<dim>::clone() const 
{
	return std::unique_ptr<BacteriumBase<dim> >(
		new PICBacterium<dim>(*this) );
}

// modify for intermittent cheaters:
template<int dim>
double 
PICBacterium<dim>::getFitness(const Bacteria::TestNewFitness::Fitness_Function<dim>& fitness_function) const
{
	return fitness_function.value(this->location, this->getSecretionRates() );
}

template<int dim>
std::vector<double> 
PICBacterium<dim>::getSecretionRates() const
{
	std::vector<double> current_rates = this->secretion_rates;
	current_rates[0] = this->getSecretionRate(0);
	return current_rates;
}

template<int dim>
double 
PICBacterium<dim>::getSecretionRate(unsigned int i) const
{
	assert(i < this->secretion_rates.size());
	if(i == 0)
		return get_current_public_good();

	return this->secretion_rates[i];
}

template<int dim>
double
PICBacterium<dim>::get_current_public_good() const
{
	const double period = on_period + off_period;
	const double tcomp = std::fmod(time, period);

	if( (time >= delay) && (tcomp < on_period) )
		return this->secretion_rates[0];
	else
		return 0;
}

template<int dim>
void 
PICBacterium<dim>::update_state(double dt, 
	const RefactoredChemicals::ChemicalHandler<dim>& /* chemicals */)
{
	time += dt;
}

template<int dim>
void 
PICBacterium<dim>::print(std::ostream& out) const
{
	out << this->location << " ";
	out << this->getSecretionRate(0) << " ";
	for(unsigned int i = 1; i < this->secretion_rates.size()-1; ++i)
		out << this->secretion_rates[i] << " ";
	out << this->secretion_rates[this->secretion_rates.size()-1];

	out << " " << time << " " << delay;
}








// ***************
// also add:
// - stochastic switching
// - feedback switching (history dependent with integration kernel, c.f chemotaxis)



// ------------------------------------------------------------------------------
// BACTERIA CLASS
// ------------------------------------------------------------------------------
/** \brief Bacteria class */
template<int dim>
class Bacteria{
public:
	Bacteria();

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm, const Geometry<dim>& geo); 
	void reintro(const ParameterHandler& prm, const Geometry<dim>& geo);

	
	void move(double dt,
				const Geometry<dim>& geo,
				const Velocity::AdvectionHandler<dim>& velocity); 

	// check using right fitness, remove/move legacy
	// void reproduce(double dt, const TestNewFitness::Fitness_Function<dim>& fitness_function); 

	void force_mutate(int n_mutate, double ds);
	void mutate(double dt);
	void update_state(double dt, 
				const RefactoredChemicals::ChemicalHandler<dim>& chemicals);

	void reproduce(
		double dt, 
		const MicrobeSimulator::Bacteria::TestNewFitness::Fitness_Function<dim>& fitness_function);


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
		double right_edge, std::vector<double>& pg_rates, int direction);

	void add_bacteria(const ParameterHandler& prm, const Geometry<dim>& geo); 
		// mixing all types for now
};

// IMPL
// ------------------------------------------------------------------------------
template<int dim>
Bacteria<dim>::Bacteria() 
	:
	diffusion_constant(0),
	mutation_rate(0),
	mutation_strength(0),
	binary_mutation(false),
	original_rate(0),
	edge_buffer(0),
	right_open_buffer(0)
{}

template<int dim>
void 
Bacteria<dim>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Bacteria");
		// prm.declare_entry("Number bacteria","100",Patterns::Unsigned());
		prm.declare_entry("Number groups","1",Patterns::Unsigned());
		prm.declare_entry("Initial growth time","-1",Patterns::Double()); // add delays instead...
		prm.declare_entry("Diffusion","0.1",Patterns::Double());
		prm.declare_entry("Edge buffer","0",Patterns::Double());
		prm.declare_entry("Right open buffer","0",Patterns::Double());
		// prm.declare_entry("Secretion rate",
		// 					"{100,100}",
		// 					Patterns::List(Patterns::Double())); // make local to type
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
		prm.declare_entry("Y start buffer","0",Patterns::Double());
	prm.leave_subsection();

	// declare parameters for each type:
	BacteriumBase<dim>::declare_parameters(prm);
	PICBacterium<dim>::declare_parameters(prm);
}

template<int dim>
void 
Bacteria<dim>::init(const ParameterHandler& prm, const Geometry<dim>& geo)
{
	const std::string section = "Bacteria";

	// assign constants:
	diffusion_constant = prm.get_double(section,"Diffusion");

	mutation_rate = prm.get_double(section, "Mutation rate");
	mutation_strength = prm.get_double(section, "Mutation strength");
	binary_mutation = prm.get_bool(section, "Binary mutation");
	/** @todo generalize this to multiple possible chemicals */
	original_rate = prm.get_double_vector(section + ".Base", "Secretion rate")[0];
	edge_buffer = prm.get_double(section,"Edge buffer");
	right_open_buffer = prm.get_double(section, "Right open buffer");

	// add bacteria:
	const unsigned int n_bact_base = 
		prm.get_unsigned(section + ".Base", "Number bacteria");
	const unsigned int n_bact_PIC = 
		prm.get_unsigned(section + ".PIC", "Number bacteria");
	bacteria.clear();
	bacteria.reserve(n_bact_base + n_bact_PIC);

	add_bacteria(prm,geo);
}

template<int dim>
void 
Bacteria<dim>::reintro(const ParameterHandler& prm, const Geometry<dim>& geo)
{
	add_bacteria(prm,geo);
}

template<int dim>
void 
Bacteria<dim>::add_bacteria(const ParameterHandler& prm, const Geometry<dim>& geo)
{
	const std::string section = "Bacteria";

	unsigned int number_groups = prm.get_unsigned(section, "Number groups");
	std::vector<Point<2> > initial_locations = 
		prm.get_point_list(section, "Initial locations");

	if(initial_locations.empty()) // if not given from parameters
	{
		if(number_groups == 0){
			const unsigned int n_bact_base = 
				prm.get_unsigned(section + ".Base", "Number bacteria");
			const unsigned int n_bact_PIC = 
				prm.get_unsigned(section + ".PIC", "Number bacteria");
			number_groups = n_bact_base + n_bact_PIC;
		}
		
		const double left_start_width = prm.get_double(section, "Left start width");
		const double left_start_buffer = prm.get_double(section, "Left start buffer");
		const double y_start_buffer = prm.get_double(section, "Y start buffer");

		initial_locations = get_bacteria_locations(geo, number_groups, 
			edge_buffer, left_start_width, left_start_buffer, y_start_buffer);
	}

	const unsigned int n_bact_base = 
		prm.get_unsigned(section + ".Base", "Number bacteria");
	const unsigned int n_bact_PIC = 
		prm.get_unsigned(section + ".PIC", "Number bacteria");
	// const unsigned int n_bact = n_bact_base + n_bact_PIC;

	// add Base bacteria
	for(unsigned int i = 0; i < n_bact_base; ++i)
	{
		unsigned int group_index = i % number_groups;	
		Point<dim> location = initial_locations[group_index];

		bacteria.emplace_back(new BacteriumBase<dim>(prm, location));
	}	

	// add PIC bacteria
	for(unsigned int i = 0; i < n_bact_PIC; ++i)
	{
		unsigned int group_index = i % number_groups;	
		Point<dim> location = initial_locations[group_index];

		bacteria.emplace_back(new PICBacterium<dim>(prm, location));
	}	
}

template<int dim>
void 
Bacteria<dim>::move(double dt,
			const Geometry<dim>& geometry,
			const Velocity::AdvectionHandler<dim>& velocity)
{
	for(unsigned int i = 0; i < bacteria.size(); ++i)
		bacteria[i]->move(dt, diffusion_constant, 
			geometry, velocity, edge_buffer);

	const int direction = velocity.get_direction();
	// if boundary is open, remove fallen bacteria:
	if(geometry.getBoundaryConditions()[0] == BoundaryCondition::OPEN)
		if(direction > 0) 
			remove_and_capture_fallen_bacteria(geometry.getTopRightPoint()[0] - right_open_buffer,
												pg_rates, direction);
		else
			remove_and_capture_fallen_bacteria(geometry.getBottomLeftPoint()[0] + right_open_buffer,
												pg_rates, direction);
}

// check using right fitness, remove/move legacy
// void reproduce(double dt, const TestNewFitness::Fitness_Function<dim>& fitness_function); 

template<int dim>
void 
Bacteria<dim>::force_mutate(int n_mutate, double ds)
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
Bacteria<dim>::mutate(double dt)
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

template<int dim>
void 
Bacteria<dim>::update_state(double dt, 
			const RefactoredChemicals::ChemicalHandler<dim>& chemicals)
{
	for(unsigned int i = 0; i < bacteria.size(); ++i)
		bacteria[i]->update_state(dt, chemicals);
}


/** \brief Remove and record public good secretion rates of fallen bacteria */
template<int dim>
void 
Bacteria<dim>::remove_and_capture_fallen_bacteria(
	double right_edge, std::vector<double>& pg_rates, int direction)
{
	const double tolerance = 1e-4;

	if(direction > 0)
	{
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
	}
	else
	{
		for(auto it = bacteria.begin(); it != bacteria.end(); )
		{
			if( (*it)->getLocation()[0] < (right_edge + tolerance) ) 
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
	}

} // remove_and_capture_fallen_bacteria()

/** \brief Reproduce bacteria according to supplied fitness function and time step */
// USING NEW CHEM AND FITNESS:
/** @todo Can probably speed up */
template<int dim>
void 
Bacteria<dim>::reproduce(
	double dt, const MicrobeSimulator::Bacteria::TestNewFitness::Fitness_Function<dim>& fitness_function)
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
Bacteria<dim>::getDiffusionConstant() const
{
	return diffusion_constant;
}

/** \brief Returns total number of current bacteria */
template<int dim>
unsigned int 
Bacteria<dim>::getTotalNumber() const
{
	return bacteria.size();
}



/** \brief Returns all locations of bacteria */
template<int dim>
std::vector<Point<dim> > 		
Bacteria<dim>::getAllLocations() const
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
Bacteria<dim>::getAllRates(unsigned int index) const
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
Bacteria<dim>::getAllRates() const
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
Bacteria<dim>::get_pg_rates()
{
	return pg_rates;
}

/** \brief Returns true is bacteria not empty */
template<int dim>
bool 
Bacteria<dim>::isAlive() const
{
	return !(bacteria.empty());
}

/** \brief Prints out bacteria info for each bacteria to given ostream object */
template<int dim>
void
Bacteria<dim>::print(std::ostream& out) const
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
Bacteria<dim>::printInfo(std::ostream& out) const
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
		<< "\t Right open buffer: " << right_open_buffer << std::endl;
	// bacteria[0].printInfo(out);

	out << std::endl << Utility::medium_line << std::endl
		<< std::endl << std::endl;
}


}} // CLOSE NAMESPACES
/* bacteria.h */