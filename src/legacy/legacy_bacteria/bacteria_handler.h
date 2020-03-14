#ifndef MICROBESIMULATOR_BACTERIA_HANDLER_H
#define MICROBESIMULATOR_BACTERIA_HANDLER_H

#include "./ic_bacterium.h"
#include "./as_bacterium.h"
#include "./reg_bacterium.h"
#include "../utility/bool_functions.h"

#include <vector>
// #include <array>

namespace MicrobeSimulator{ namespace MultiBacteria{

template<int dim>
class BacteriaHandler
{
public:
	BacteriaHandler();

	void init(double db, unsigned int n_reg, unsigned int n_ic, unsigned int n_as,
		double switching_rate, double pg_rate, double waste_rate); // done
	// intialized at origin

	void init(double db, unsigned int n_reg, unsigned int n_ic, unsigned int n_as,
		double on_period, double off_period, 
		double pg_rate, double waste_rate, 
		const std::vector<Point<dim> >& group_locations); // done
	// partition each type equally into each group

	void init(double db, unsigned int n_reg, unsigned int n_ic, unsigned int n_as,
		double switching_rate, double pg_rate, double waste_rate, 
		const std::vector<Point<dim> >& group_locations); // done

	// accessors:
	double getDiffusionConstant() const; // done
	unsigned int getNumberRegular() const; // done
	unsigned int getNumberIC() const; // done
	unsigned int getNumberAS() const; // done -- for now
	unsigned int getTotalNumber() const; // done - for now

	std::vector<Point<dim> >	getAllLocations() const; // done
	std::vector<double> 		getAllGoodRates() const; // done
	std::vector<double> 		getAllWasteRates() const; // done

	bool isAlive() const; // done

	void print(std::ostream& out) const; // done
	void printInfo(std::ostream& out) const; // done

	// modifiers:
	void randomWalk(double dt, const Geometry<dim>& geometry,
		const Velocity::AdvectionHandler<dim>& velocity); // done

	void updateTimeAndSecretion(double dt); // done -- may use to update AS suicide too

	void reproduce(double dt, const FitnessBase<dim, 2>& fitness_function); //***
	void reproduce_simple(double dt, const FitnessBase<dim, 2>& fitness_function); // done

	/**
	* Not the most elegant solution, but could use the second ``rate'' for the chemical
	* ratio threshold in a separate fitness function for AS_bacteria...
	*/

	// mutation:
	void mutate(double dt, double sp_mu, double reg_mu, double ds); 
	void mutate_binary(double dt, double sp_mu, double reg_mu, 
		double original_secretion_rate); 
private:
	double diffusion_constant;

	BoolFunctions::SquareWave 	switch_function;

	std::vector<RegBacterium<dim> > reg_bacteria;
	std::vector<ICBacterium<dim> > ic_bacteria;
	// std::vector<ASBacterium<dim> > as_bacterium; // worry about later...
	void mutate_base(double dt, double sp_mu, 
		double reg_mu, double delta_secretion, 
		bool mutate_binary, double original_secretion_rate);
};

// IMPLEMENTATION:
// ----------------------------------------------------------------------------
template<int dim>
BacteriaHandler<dim>::BacteriaHandler()
	:
	diffusion_constant(0)
	// switching_rate(0) 
{}

template<int dim>
void 
BacteriaHandler<dim>::init(double db, 
	unsigned int n_reg, unsigned int n_ic, unsigned int n_as,
	double switching_rate, double pg_rate, double waste_rate)
{
	reg_bacteria.clear();
	ic_bacteria.clear();

	reg_bacteria.reserve(n_reg);
	ic_bacteria.reserve(n_ic);

	diffusion_constant = db;

	for(unsigned int i = 0; i < n_reg; ++i)
		reg_bacteria.emplace_back( RegBacterium<dim>(pg_rate, waste_rate) );

	for(unsigned int i = 0; i < n_ic; ++i)
		ic_bacteria.emplace_back( ICBacterium<dim>(pg_rate, waste_rate) );
} // intialized at origin


template<int dim>
void 
BacteriaHandler<dim>::init(double db, 
	unsigned int n_reg, unsigned int n_ic, unsigned int n_as,
	double switching_rate, double pg_rate, double waste_rate, 
	const std::vector<Point<dim> >& group_locations)
{
	reg_bacteria.clear();
	ic_bacteria.clear();

	reg_bacteria.reserve(n_reg);
	ic_bacteria.reserve(n_ic);

	diffusion_constant = db;

	const unsigned int number_groups = group_locations.size();

	for(unsigned int i = 0; i < n_reg; ++i)
	{
		unsigned int group_index = i % number_groups;
		reg_bacteria.emplace_back( 
			RegBacterium<dim>(group_locations[group_index],pg_rate, waste_rate) );
	}

	for(unsigned int i = 0; i < n_ic; ++i)
	{
		unsigned int group_index = i % number_groups;
		ic_bacteria.emplace_back( 
			ICBacterium<dim>(group_locations[group_index], pg_rate, waste_rate) );
	}

} // partition each type equally into each group


template<int dim>
void 
BacteriaHandler<dim>::init(double db, unsigned int n_reg, unsigned int n_ic, unsigned int n_as,
	double on_period, double off_period, 
	double pg_rate, double waste_rate, 
	const std::vector<Point<dim> >& group_locations)
{
	reg_bacteria.clear();
	ic_bacteria.clear();

	reg_bacteria.reserve(n_reg);
	ic_bacteria.reserve(n_ic);

	diffusion_constant = db;

	const unsigned int number_groups = group_locations.size();

	for(unsigned int i = 0; i < n_reg; ++i)
	{
		unsigned int group_index = i % number_groups;
		reg_bacteria.emplace_back( 
			RegBacterium<dim>(group_locations[group_index],pg_rate, waste_rate) );
	}

	for(unsigned int i = 0; i < n_ic; ++i)
	{
		unsigned int group_index = i % number_groups;
		ic_bacteria.emplace_back( 
			ICBacterium<dim>(group_locations[group_index], pg_rate, waste_rate) );
	}

	switch_function.setOnPeriod(on_period);
	switch_function.setOffPeriod(off_period);
}


// accessors:
template<int dim>
double 
BacteriaHandler<dim>::getDiffusionConstant() const
{
	return diffusion_constant;
}

template<int dim>
unsigned int 
BacteriaHandler<dim>::getNumberRegular() const
{
	return reg_bacteria.size();
}

template<int dim>
unsigned int 
BacteriaHandler<dim>::getNumberIC() const
{
	return ic_bacteria.size();
}

template<int dim>
unsigned int 
BacteriaHandler<dim>::getNumberAS() const
{
	return 0;
}

template<int dim>
unsigned int 
BacteriaHandler<dim>::getTotalNumber() const
{
	return reg_bacteria.size() + ic_bacteria.size() + 0;
}

template<int dim>
std::vector<Point<dim> >	
BacteriaHandler<dim>::getAllLocations() const
{
	const unsigned int total_size = getTotalNumber();
	std::vector<Point<dim> > all_locations;
	all_locations.reserve(total_size);

	// regular bacteria first:
	for(unsigned int i = 0; i < reg_bacteria.size(); ++i)
		all_locations.emplace_back(reg_bacteria[i].getLocation());

	// ic bacteria next:
	for(unsigned int i = 0; i < ic_bacteria.size(); ++i)
		all_locations.emplace_back(ic_bacteria[i].getLocation());

	// finally as bacteria:
	// ...

	return all_locations;
}

template<int dim>
std::vector<double> 		
BacteriaHandler<dim>::getAllGoodRates() const
{
	const unsigned int total_size = getTotalNumber();
	std::vector<double> all_good_rates;
	all_good_rates.reserve(total_size);

	// regular bacteria first:
	for(unsigned int i = 0; i < reg_bacteria.size(); ++i)
		all_good_rates.emplace_back(reg_bacteria[i].getGoodSecretion());

	// ic bacteria next:
	for(unsigned int i = 0; i < ic_bacteria.size(); ++i)
		all_good_rates.emplace_back(ic_bacteria[i].getGoodSecretion());

	// finally as bacteria:
	// ...

	return all_good_rates;
}

template<int dim>
std::vector<double> 		
BacteriaHandler<dim>::getAllWasteRates() const
{
	const unsigned int total_size = getTotalNumber();
	std::vector<double> all_waste_rates;
	all_waste_rates.reserve(total_size);

	// regular bacteria first:
	for(unsigned int i = 0; i < reg_bacteria.size(); ++i)
		all_waste_rates.emplace_back(reg_bacteria[i].getWasteSecretion());

	// ic bacteria next:
	for(unsigned int i = 0; i < ic_bacteria.size(); ++i)
		all_waste_rates.emplace_back(ic_bacteria[i].getWasteSecretion());

	// finally as bacteria:
	// ...

	return all_waste_rates;
}

template<int dim>
bool 
BacteriaHandler<dim>::isAlive() const
{
	return !(reg_bacteria.empty() && ic_bacteria.empty());
}

template<int dim>
void 
BacteriaHandler<dim>::print(std::ostream& out) const
{
	out << reg_bacteria.size() << std::endl;
	for(unsigned int i = 0; i < reg_bacteria.size(); ++i)
		reg_bacteria[i].print(out);

	out << ic_bacteria.size() << std::endl;
	for(unsigned int i = 0; i < ic_bacteria.size(); ++i)
		ic_bacteria[i].print(out);
}

template<int dim>
void 
BacteriaHandler<dim>::printInfo(std::ostream& out) const
{
	out << "Mutli-type bacteria" << std::endl
		<< "\t Diffusion constant: " << diffusion_constant << std::endl
		<< "\t Public good secretion: " << reg_bacteria[0].getGoodSecretion() << std::endl
		<< "\t Waste secretion: " << reg_bacteria[0].getWasteSecretion() << std::endl
		<< "\t Number regular bacteria: " << reg_bacteria.size() << std::endl
		<< "\t Number IC bacteria: " << ic_bacteria.size() << std::endl;

	// *** can print switch function info ***
	// if(!ic_bacteria.empty())
	// 	out << "\t Switching rate: " << ic_bacteria[0].getSwitchingRate() << std::endl;

	out << "\t Number AS bacteria: " << 0 << std::endl;
}

// modifiers:
template<int dim>
void 
BacteriaHandler<dim>::randomWalk(double dt, const Geometry<dim>& geometry,
	const Velocity::AdvectionHandler<dim>& velocity)
{
	// first regular bacteria:
	for(unsigned int i = 0; i < reg_bacteria.size(); ++i)
		reg_bacteria[i].randomStep(dt, diffusion_constant,
			geometry, velocity);

	// then IC bacteria:
	for(unsigned int i = 0; i < ic_bacteria.size(); ++i)
		ic_bacteria[i].randomStep(dt, diffusion_constant,
			geometry, velocity);

	// finally AS bacteria:
	// ...
}

template<int dim>
void 
BacteriaHandler<dim>::updateTimeAndSecretion(double dt) // may use to update AS suicide too
{
	for(unsigned int i = 0; i < ic_bacteria.size(); ++i)
		ic_bacteria[i].updateTimeAndSecretion(dt, switch_function);
}

template<int dim>
void 
BacteriaHandler<dim>::reproduce(double dt,
	const FitnessBase<dim, 2>& fitness_function)
{

}

template<int dim>
void 
BacteriaHandler<dim>::reproduce_simple(double dt,
	const FitnessBase<dim, 2>& fitness_function)
{
	std::vector<RegBacterium<dim> > reg_offspring;

	for(auto it_reg = reg_bacteria.begin();
		it_reg != reg_bacteria.end(); )
	{
		const double fit = it_reg->getFitness(fitness_function)*dt;
		if(fit < 0)
		{
			it_reg = reg_bacteria.erase(it_reg);
		}
		else
		{
			double prob = ((double) rand() / (RAND_MAX));
			if(prob<fit)
			{
				reg_offspring.emplace_back(*it_reg);
			}
			++it_reg;
		}
	}

	// add offspring to end:
	reg_bacteria.insert(reg_bacteria.end(),reg_offspring.begin(),reg_offspring.end());

	// same for ic bacteria:
	std::vector<ICBacterium<dim> > ic_offspring;

	for(auto it_ic = ic_bacteria.begin();
		it_ic != ic_bacteria.end(); )
	{
		const double fit = it_ic->getFitness(fitness_function)*dt;
		if(fit < 0)
		{
			it_ic = ic_bacteria.erase(it_ic);
		}
		else
		{
			double prob = ((double) rand() / (RAND_MAX));
			if(prob<fit)
			{
				ic_offspring.emplace_back(*it_ic);
			}
			++it_ic;
		}
	}

	// add offspring to end:
	ic_bacteria.insert(ic_bacteria.end(),ic_offspring.begin(),ic_offspring.end());
}

// would be quicker to have a type that we change ...
// instead of killing and reproducing ... (can use swap to end to speed up too...)
template<int dim>
void 
BacteriaHandler<dim>::mutate(double dt, double sp_mu, 
	double reg_mu, double delta_secretion)
{
	mutate_base(dt, sp_mu, reg_mu, delta_secretion, false, 0); 
}

template<int dim>
void 
BacteriaHandler<dim>::mutate_binary(double dt, double sp_mu, 
	double reg_mu, double original_secretion_rate)
{
	mutate_base(dt, sp_mu, reg_mu, 0, true, original_secretion_rate); 
}

template<int dim>
void 
BacteriaHandler<dim>::mutate_base(double dt, double sp_mu, 
	double reg_mu, double delta_secretion, 
	bool mutate_binary, double original_secretion_rate)
{
	std::vector<RegBacterium<dim> > 	reg_mutants;
	std::vector<ICBacterium<dim> > 		ic_mutants;
	// as_mutants ...

	// loop over regular bacteria:
	for(auto it_reg = reg_bacteria.begin();
			it_reg != reg_bacteria.end(); )
	{
		// mutate secretion rate:
		double prob = ((double) rand() / (RAND_MAX));
		if(prob < (reg_mu*dt) )
		{
			if(mutate_binary == false)
			{
				const double ds = delta_secretion*(
					2.*((double) rand() / (RAND_MAX)) - 1); 

				it_reg->setGoodSecretion( it_reg->getGoodSecretion() + ds );
				if(it_reg->getGoodSecretion() < 0)
					it_reg->setGoodSecretion(0); // absorbing lower bc
			}
			else
			{
				if(it_reg->getGoodSecretion() != 0)
					it_reg->setGoodSecretion(0);
				else
					it_reg->setGoodSecretion(original_secretion_rate);
			} // else binary mutation

		} // if mutating secretion rate

		double prob2 = ((double) rand() / (RAND_MAX));
		if(prob2 < (sp_mu*dt) )
		{
			ic_mutants.emplace_back( 
				ICBacterium<dim>(it_reg->getLocation(), it_reg->getGoodSecretion(),
					it_reg->getWasteSecretion() ) );

			it_reg = reg_bacteria.erase(it_reg);
		}
		else
		{
			++it_reg;
		}
	} // for each regular bacteria

	// loop over IC bacteria:
	for(auto it_ic = ic_bacteria.begin();
			it_ic != ic_bacteria.end(); )
	{
		double prob2 = ((double) rand() / (RAND_MAX));
		if(prob2 < (sp_mu*dt) )
		{
			reg_mutants.emplace_back( 
				RegBacterium<dim>(it_ic->getLocation(), it_ic->getGoodSecretion(),
					it_ic->getWasteSecretion() ) );

			it_ic = ic_bacteria.erase(it_ic);
		}
		else
		{
			++it_ic;
		}
	}

	// add mutants to end:
	reg_bacteria.insert(reg_bacteria.end(), reg_mutants.begin(), reg_mutants.end());	
	ic_bacteria.insert(ic_bacteria.end(), ic_mutants.begin(), ic_mutants.end());

} // mutate()


}} // close namespace
#endif