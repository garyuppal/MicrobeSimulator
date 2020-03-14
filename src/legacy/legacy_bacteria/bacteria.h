#ifndef MICROBESIMULATOR_BACTERIA_H
#define MICROBESIMULATOR_BACTERIA_H

#include <deal.II/base/point.h>
using dealii::Point;

#include "../geometry/geometry.h"
#include "../advection/advection_handler.h"
#include "../fitness/fitness_base.h"

#include <iostream>
#include <string>
#include <sstream>
#include <array>

namespace MicrobeSimulator{ 

template<int dim, int numchem>
class Bacteria{
public:
	Bacteria();

	void init(double db, unsigned int number_bacteria,
		const std::array<double, numchem>& rates); // done
	void init(double db, unsigned int number_bacteria,
		const std::array<double, numchem>& rates,
		const std::vector<Point<dim> >& groups); // done
	void init(double db, 
		const std::vector<std::pair<Point<dim>, unsigned int> >& coarse_bacteria,
		const std::array<double, numchem>& rates); // done

	// reinitialize (for 2 chem only)
	void reinit(double db, 
		const std::vector<double>& pg_rates, // of length number_bacteria
		double waste_rate,
		const std::vector<Point<dim> >& groups); 

	// reintroduce new bacteria (for 2 chem only)
	void reintro(unsigned int n, double pg_rate, double waste_rate, 
		const std::vector<Point<dim> >& p,
		const Geometry<dim>& geometry,
		double time_step);

	// accesors:
	const std::vector<Point<dim> >& getLocations() const; // done
	const std::vector<double>& 		getSecretionRates(unsigned int i) const; // done
	std::array<double, numchem> 	getBacteriaSecretionRates(unsigned int b) const;
	unsigned int getSize() const; // done
	bool isAlive() const; // done

	void print(std::ostream& out) const; // done
	void printInfo(std::ostream& out) const; // done

	// modifiers:
	void randomWalk(double dt, 
		const Geometry<dim>& geometry, const Velocity::AdvectionHandler<dim>& velocity); // done
	void randomWalk(double dt, 
		const Geometry<dim>& geometry, const Velocity::AdvectionHandler<dim>& velocity,
		std::vector<double>& pg_rates); // done

	void reproduce_simple(double dt, const FitnessBase<dim, numchem>& fitness_function);
	void reproduce_simple2(double dt, const FitnessBase<dim, numchem>& fitness_function);

	void reproduce(double dt, const FitnessBase<dim, numchem>& fitness_function);
	void mutate(unsigned int id, 
		double dt, double mutation_rate, double delta_secretion);
	void mutate_binary(unsigned int id,
		double dt, double mutation_rate, double original_secretion_rate);
	// void make_cheaters(unsigned int number_cheaters); 
			// what about for multiple secretion rates?

private:
	double diffusion_constant;

	std::vector<Point<dim> > locations;
	std::array<std::vector<double>, numchem> secretion_rates;

	// methods:
	void set_secretion_rates(const std::array<double, numchem>& rates); // done
	void remove_fallen_bacteria(double right_edge); // done -- double checked, triple check...
	void remove_and_store_fallen(double right_edge, std::vector<double>& pg_rates);

	// for reproducing:
	void copy_bacteria_to_this(unsigned int id_this, unsigned int id_copy); // done
	void copy_bacteria_to_end(unsigned int i); // done
	void clear_bacteria(); // done
	void reserve_bacteria(unsigned int number_bacteria);
	void remove_end(unsigned int last_living); // done
};


// IMPLEMENTATION:
//---------------------------------------------------------------------------------------------
template<int dim, int numchem>
Bacteria<dim, numchem>::Bacteria()
{}

// INTIALIZATION:
template<int dim, int numchem>
void 
Bacteria<dim, numchem>::init(double db, unsigned int number_bacteria,
	const std::array<double, numchem>& rates)
{
	diffusion_constant = db;

	locations.clear();
	locations.reserve(number_bacteria);

	for(unsigned int i = 0; i < number_bacteria; ++i)
		locations.emplace_back(Point<dim>()); // init at origin

	set_secretion_rates(rates);
}

template<int dim, int numchem>
void 
Bacteria<dim, numchem>::init(double db, unsigned int number_bacteria,
	const std::array<double, numchem>& rates,
	const std::vector<Point<dim> >& groups)
{
	diffusion_constant = db;

	locations.clear();
	locations.reserve(number_bacteria);

	const unsigned int number_groups = groups.size();

	for(unsigned int i=0; i < number_bacteria; ++i)
	{
		unsigned int group_index = i % number_groups;	
		locations.emplace_back( groups[group_index] );
	} // for number of bacteria

	set_secretion_rates(rates);
}

template<int dim, int numchem>
void 
Bacteria<dim, numchem>::init(double db, 
	const std::vector<std::pair<Point<dim>, unsigned int> >& coarse_bacteria,
	const std::array<double, numchem>& rates)
{
	diffusion_constant = db;

	unsigned int number_bacteria = 0;
	const unsigned int n_points = coarse_bacteria.size();
	for(unsigned int i = 0; i < n_points; ++i)
		number_bacteria += coarse_bacteria[i].second;

	locations.clear();

	if(number_bacteria == 0)
	{
		std::cout << "No bacteria intialized from continuous solution" << std::endl
			<< "Double check continuous parameters and solution" << std::endl;
		return;
	}

	std::cout << "... intializing " << number_bacteria << " total bacteria" << std::endl;

	locations.reserve(number_bacteria);

	for(unsigned int i = 0; i < n_points; ++i)
	{
		const unsigned int n_bact = coarse_bacteria[i].second;
		for(unsigned int b = 0; b < n_bact; ++b)
			locations.emplace_back( coarse_bacteria[i].first );
	}

	set_secretion_rates(rates);
}


template<int dim, int numchem>
void
Bacteria<dim, numchem>::set_secretion_rates(const std::array<double, numchem>& rates)
{
	const unsigned int number_bacteria = locations.size();

	for(unsigned int c = 0; c < numchem; ++c)
	{
		secretion_rates[c].clear();
		secretion_rates[c].reserve(number_bacteria);

		for(unsigned i = 0; i < number_bacteria; ++i)
			secretion_rates[c].emplace_back(rates[c]);
	}
}

// reinitialize (for 2 chem only)
template<int dim, int numchem>
void 
Bacteria<dim, numchem>::reinit(double db, 
	const std::vector<double>& pg_rates, // of length number_bacteria
	double waste_rate,
	const std::vector<Point<dim> >& groups)
{
	if(numchem != 2)
		throw std::runtime_error("Re-intialization setup for 2 chemicals only right now");

	clear_bacteria(); 
	const unsigned int number_bacteria = pg_rates.size();
	if(number_bacteria == 0)
		return;

	reserve_bacteria(number_bacteria);
	diffusion_constant = db; // should really already be set, but need not be

	const unsigned int number_groups = groups.size();

	for(unsigned int i = 0; i < number_bacteria; ++i)
	{
		unsigned int group_index = i % number_groups;
		locations.emplace_back( groups[group_index] );

		secretion_rates[0].emplace_back( pg_rates[i] ); // add public good rate
		secretion_rates[1].emplace_back( waste_rate ); // add waste rate
	}

}

// reintroduce new groups:
template<int dim, int numchem>
void 
Bacteria<dim, numchem>::reintro(unsigned int n, 
	double pg_rate, double waste_rate, const std::vector<Point<dim> >& p,
	const Geometry<dim>& geometry, double time_step)
{
	const double dt = 0.1*time_step;

	const unsigned int num_groups = p.size();
	std::cout << "reintroducing " << n << " bacteria in "
		<< num_groups << " groups..." << std::endl;

	// create new bacteria in n groups:
	for(unsigned int i = 0 ; i < n; ++i)
	{
		unsigned int group_index = i % num_groups;
		Point<dim> loc = p[group_index];

		// spread out locations:
		for(unsigned int j = 0; j < 20; ++j)
		{
			const double theta = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));
			const double phi = dealii::numbers::PI*((double) rand() / (RAND_MAX));

			Point<dim> randomPoint = (dim == 2) ? Point<dim>(std::cos(theta),std::sin(theta))
			                        : Point<dim>(std::cos(phi), std::sin(phi)*std::cos(theta),
			                            std::sin(phi)*std::sin(theta));

			loc += std::sqrt(2*dim*dt*diffusion_constant)*randomPoint;

			geometry.checkBoundaries(p[group_index], loc); 
		}

		locations.emplace_back(loc);
		secretion_rates[0].emplace_back(pg_rate);
		secretion_rates[1].emplace_back(waste_rate);
	}

	// spread out new ones only: ... might be easier to spread points first ...

}

// MODIFIERS:
//------------------------------------------------------------------------------------
template<int dim, int numchem>
void 
Bacteria<dim, numchem>::randomWalk(double dt, 
	const Geometry<dim>& geometry, const Velocity::AdvectionHandler<dim>& velocity)
{
	for(unsigned int i = 0; i < locations.size(); ++i)
	{
        Point<dim> old_location = locations[i];

		const double theta = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));
		const double phi = dealii::numbers::PI*((double) rand() / (RAND_MAX));

		Point<dim> randomPoint = (dim == 2) ? Point<dim>(std::cos(theta),std::sin(theta))
		                        : Point<dim>(std::cos(phi), std::sin(phi)*std::cos(theta),
		                            std::sin(phi)*std::sin(theta));

		locations[i] += std::sqrt(2*dim*dt*diffusion_constant)*randomPoint
			+ dt*velocity.value(locations[i]);

		geometry.checkBoundaries(old_location, locations[i]); 
	}

	// check right boundary if needed
	if(geometry.getBoundaryConditions()[0] == BoundaryCondition::OPEN)
		remove_fallen_bacteria(geometry.getTopRightPoint()[0]);	
}

template<int dim, int numchem>
void 
Bacteria<dim, numchem>::randomWalk(double dt, 
	const Geometry<dim>& geometry, 
	const Velocity::AdvectionHandler<dim>& velocity,
	std::vector<double>& pg_rates)
{
	for(unsigned int i = 0; i < locations.size(); ++i)
	{
        Point<dim> old_location = locations[i];

		const double theta = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));
		const double phi = dealii::numbers::PI*((double) rand() / (RAND_MAX));

		Point<dim> randomPoint = (dim == 2) ? Point<dim>(std::cos(theta),std::sin(theta))
		                        : Point<dim>(std::cos(phi), std::sin(phi)*std::cos(theta),
		                            std::sin(phi)*std::sin(theta));

		locations[i] += std::sqrt(2*dim*dt*diffusion_constant)*randomPoint
			+ dt*velocity.value(locations[i]);

		geometry.checkBoundaries(old_location, locations[i]); 
	}

	// check right boundary if needed
	if(geometry.getBoundaryConditions()[0] == BoundaryCondition::OPEN)
		remove_and_store_fallen(geometry.getTopRightPoint()[0], pg_rates);	
}

template<int dim, int numchem>
void 
Bacteria<dim,numchem>::reproduce_simple(double dt, 
	const FitnessBase<dim, numchem>& fitness_function)
{
	// for 2 chemicals
	auto it_loc = locations.begin();
	auto it_good = secretion_rates[0].begin();
	auto it_waste = secretion_rates[1].begin();

	std::vector<Point<dim> > offspring_locations;
	std::vector<double> offspring_good;
	std::vector<double> offspring_waste;

	for( ; it_loc != locations.end(); )
	{
		const std::array<double, 2> rates({*it_good, *it_waste});
		const double fit = dt*fitness_function.value( *it_loc, rates );

		 if(fit < 0)
		 {
		 	it_loc = locations.erase(it_loc);
		 	it_good = secretion_rates[0].erase(it_good);
		 	it_waste = secretion_rates[1].erase(it_waste);
		 }
		 else
		 {
		 	double prob = ((double) rand() / (RAND_MAX));
			if(prob<fit)
			{
				offspring_locations.emplace_back(*it_loc);
				offspring_good.emplace_back(*it_good);
				offspring_waste.emplace_back(*it_waste);
			}
			++it_loc;
			++it_good;
			++it_waste;
		 } // else possibly reproduce
	} // for all bacteria

	// add offspring to end:
	locations.insert(locations.end(), offspring_locations.begin(), offspring_locations.end());
	secretion_rates[0].insert(secretion_rates[0].end(), offspring_good.begin(), offspring_good.end());
	secretion_rates[1].insert(secretion_rates[1].end(), offspring_waste.begin(), offspring_waste.end());
}

template<int dim, int numchem>
void 
Bacteria<dim,numchem>::reproduce_simple2(double dt, 
	const FitnessBase<dim, numchem>& fitness_function)
{
	// for 2 chemicals
	auto it_loc = locations.begin();
	auto it_good_one = secretion_rates[0].begin();
	auto it_good_two = secretion_rates[1].begin();
	auto it_waste = secretion_rates[2].begin();

	std::vector<Point<dim> > offspring_locations;
	std::vector<double> offspring_good_one;
	std::vector<double> offspring_good_two;
	std::vector<double> offspring_waste;

	for( ; it_loc != locations.end(); )
	{
		const std::array<double, 3> rates({*it_good_one, *it_good_two, *it_waste});
		const double fit = dt*fitness_function.value( *it_loc, rates );

		 if(fit < 0)
		 {
		 	it_loc = locations.erase(it_loc);
		 	it_good_one = secretion_rates[0].erase(it_good_one);
		 	it_good_two = secretion_rates[1].erase(it_good_two);
		 	it_waste = secretion_rates[2].erase(it_waste);
		 }
		 else
		 {
		 	double prob = ((double) rand() / (RAND_MAX));
			if(prob<fit)
			{
				offspring_locations.emplace_back(*it_loc);
				offspring_good_one.emplace_back(*it_good_one);
				offspring_good_two.emplace_back(*it_good_two);
				offspring_waste.emplace_back(*it_waste);
			}
			++it_loc;
			++it_good_one;
			++it_good_two;
			++it_waste;
		 } // else possibly reproduce
	} // for all bacteria

	// add offspring to end:
	locations.insert(locations.end(), offspring_locations.begin(), offspring_locations.end());
	secretion_rates[0].insert(secretion_rates[0].end(), offspring_good_one.begin(), offspring_good_one.end());
	secretion_rates[1].insert(secretion_rates[1].end(), offspring_good_two.begin(), offspring_good_two.end());
	secretion_rates[2].insert(secretion_rates[2].end(), offspring_waste.begin(), offspring_waste.end());
}

template<int dim, int numchem>
void 
Bacteria<dim, numchem>::reproduce(double dt, 
	const FitnessBase<dim, numchem>& fitness_function)
{
	const unsigned int size = locations.size();
	unsigned int last_index = size - 1; 
	unsigned int last_living = last_index;

	// only need to loop over original n bacteria (not over those moved to end)
	for(unsigned int i = 0; i != last_index; ) 
	{
		const double fit = dt*fitness_function.value(locations[i], 
			getBacteriaSecretionRates(i));
		// decide what to do:
		if(fit < 0)
		{
			// swap to end -- not including those already swapped: 
				//(don't need to copy this to end, just delete)
			if(last_living > last_index) // if new offspring added
			{
				copy_bacteria_to_this(i, last_living);
				++i; // since an already reproduced bacteria is replacing this
				--last_living; // since moved dead to end 
			}
			else // should have last index == last_living
			{
				copy_bacteria_to_this(i, last_living); //bacteria[i] = bacteria[last_living];
				--last_index; // since its processed 
					//-- instead of decrementing i, want to reprocess this
				--last_living; // to keep living == last 
			}
		}
		else
		{
			double prob =  ((double) rand() / (RAND_MAX));
			if(prob<fit)
			{
				// reproduce: (either emplace copy to end or replace dead)
				// careful ... size changes ...

				// need to distinguised between pushed back offspring and those at end

				if(last_living < (locations.size()-1) ) //if dead already at end
				{
					copy_bacteria_to_this(last_living + 1, i);
					// bacteria[last_living+1] = bacteria[i]; // replace
					++last_living;
				}
				else
				{
					copy_bacteria_to_end(i);
					// bacteria.emplace_back(bacteria[i]); // add copy to end
					++last_living; // to keep with total size
				}
			}
			++i; // move up if either reproduced or did nothing
		} //else possible reproduce
	} // loop once over bacteria vector

	// erase dead bacteria:
	if(last_living == 0) 
		clear_bacteria();
	else
		remove_end(last_living);
}

template<int dim, int numchem>
void
Bacteria<dim, numchem>::clear_bacteria()
{
	locations.clear();

	for(unsigned int c = 0; c < numchem; ++c)
		secretion_rates[c].clear();
}

template<int dim, int numchem>
void 
Bacteria<dim, numchem>::reserve_bacteria(unsigned int number_bacteria)
{
	locations.reserve(number_bacteria);
	for(unsigned int c = 0; c < numchem; ++c)
		secretion_rates[c].reserve(number_bacteria);
}

template<int dim, int numchem>
void
Bacteria<dim,numchem>::remove_end(unsigned int last_living)
{
	// bacteria.erase(bacteria.begin() + last_living, bacteria.end()); // check...
	locations.erase(locations.begin() + last_living, locations.end());

	for(unsigned int c = 0; c < numchem; ++c)
		secretion_rates[c].erase( secretion_rates[c].begin() + last_living,
									secretion_rates[c].end() );
}

template<int dim, int numchem>
void
Bacteria<dim, numchem>::copy_bacteria_to_this(unsigned int id_this, unsigned int id_copy)
{
	locations[id_this] = locations[id_copy]; 
	for(unsigned int c = 0; c < numchem; ++c)
		secretion_rates[c][id_this]  = secretion_rates[c][id_copy];
}

template<int dim, int numchem>
void
Bacteria<dim, numchem>::copy_bacteria_to_end(unsigned int i)
{
	locations.emplace_back(locations[i]);
	for(unsigned int c = 0; c < numchem; ++c)
		secretion_rates[c].emplace_back( secretion_rates[c][i] );
}

template<int dim, int numchem>
void 
Bacteria<dim, numchem>::mutate(unsigned int id, 
	double dt, double mutation_rate, double delta_secretion)
{
	const unsigned int number_bacteria = locations.size();

	for(unsigned int i = 0; i < number_bacteria; ++i)
	{
		double prob = ((double) rand() / (RAND_MAX));
		if(prob < (mutation_rate*dt) )
		{
			const double ds = delta_secretion*(
				2.*((double) rand() / (RAND_MAX)) - 1); 
				// ... use normal distribution instead
			secretion_rates[id][i] += ds;
			if(secretion_rates[id][i] < 0)
				secretion_rates[id][i] = 0; // absorbing lower bc
		} // if mutating
	} // for each bacteria
} // random strength mutation

template<int dim, int numchem>
void 
Bacteria<dim, numchem>::mutate_binary(unsigned int id,
	double dt, double mutation_rate, double original_secretion_rate)
{
	const unsigned int number_bacteria = locations.size();

	for(unsigned int i = 0; i < number_bacteria; ++i)
	{
		double prob = ((double) rand() / (RAND_MAX));
		if(prob < (mutation_rate*dt) )
		{
			if( secretion_rates[id][i] != 0)
				secretion_rates[id][i] = 0;
			else
				secretion_rates[id][i] = original_secretion_rate;
		} // if mutating
	} // for each bacteria
} // binary mutation


template<int dim, int numchem>
void
Bacteria<dim, numchem>::remove_fallen_bacteria(double right_edge)
{
	const double tolerance = 1e-4;
	const unsigned int original_number_bacteria = locations.size();
	unsigned int number_bacteria = locations.size();
	unsigned int number_removed = 0;

	// use swap to end method:
	for(unsigned int i = 0; i < number_bacteria; )
		if( std::fabs( locations[i][0] - right_edge ) < tolerance 
			|| locations[i][0] > right_edge)
		{
			const unsigned int end = original_number_bacteria - number_removed - 1;

			locations[i] = locations[end];
			for(unsigned int c = 0; c < numchem; ++c)
				secretion_rates[c][i] = secretion_rates[c][end];

			++number_removed;
			--number_bacteria;
		}
		else
		{
			++i;
		}

	// remove end:
	locations.erase(locations.end() - number_removed, locations.end());
	for(unsigned int c = 0; c < numchem; ++c)
		secretion_rates[c].erase(secretion_rates[c].end() - number_removed,
									secretion_rates[c].end());
}

template<int dim, int numchem>
void 
Bacteria<dim, numchem>::remove_and_store_fallen(double right_edge, 
											std::vector<double>& pg_rates)
{
	const double tolerance = 1e-4;
	const unsigned int original_number_bacteria = locations.size();
	unsigned int number_bacteria = locations.size();
	unsigned int number_removed = 0;

	// use swap to end method; store public good secretion rate:
	for(unsigned int i = 0; i < number_bacteria; )
		if( locations[i][0] > right_edge - tolerance )
		{
			const unsigned int end = original_number_bacteria - number_removed - 1;

			// if(end < 0 || end >= locations.size() )
			// {
			// 	std::cout << "Number bacteria: " << number_bacteria
			// 		<< " | number_removed: " << number_removed
			// 		<< " | i = " << i << std::endl;
			// 	std::cout << "end index out of range: " << end << std::endl;
			// 	std::cout << "locations size: " << locations.size();
			// 	for(unsigned int c = 0; c < numchem; ++c)
			// 		std::cout << " | secretion_rates " << c 
			// 			<< " size: " << secretion_rates[c].size();
			// 	std::cout << std::endl << std::endl;
			// }

			// store public good secretion rate:
			pg_rates.emplace_back( secretion_rates[0][i] );

			locations[i] = locations[end];
			for(unsigned int c = 0; c < numchem; ++c)
				secretion_rates[c][i] = secretion_rates[c][end];

			++number_removed;
			--number_bacteria;
		}
		else
		{
			++i;
		}

	// remove end:
	locations.erase(locations.end() - number_removed, locations.end());
	for(unsigned int c = 0; c < numchem; ++c)
		secretion_rates[c].erase(secretion_rates[c].end() - number_removed,
									secretion_rates[c].end());
}


// ACCESSORS:
//------------------------------------------------------------------------------------
template<int dim, int numchem>
const std::vector<Point<dim> >& 
Bacteria<dim, numchem>::getLocations() const
{
	return locations;
}

template<int dim, int numchem>
const std::vector<double>& 		
Bacteria<dim, numchem>::getSecretionRates(unsigned int i) const
{
	return secretion_rates[i];
}

template<int dim, int numchem>
std::array<double, numchem> 	
Bacteria<dim, numchem>::getBacteriaSecretionRates(unsigned int b) const
{
	std::array<double, numchem> rates;
	for(unsigned int c = 0; c < numchem; ++c)
		rates[c] = secretion_rates[c][b];

	return rates;
}

template<int dim, int numchem>
unsigned int 
Bacteria<dim, numchem>::getSize() const
{
	return locations.size();
}

template<int dim, int numchem>
bool 
Bacteria<dim, numchem>::isAlive() const
{
	return !locations.empty();
}

template<int dim, int numchem>
void Bacteria<dim, numchem>::print(std::ostream &out) const
{
	for(unsigned int i = 0; i < locations.size(); i++)
	{
		out << locations[i];
		for(unsigned int c = 0; c < numchem; ++c)
			out << " " << secretion_rates[c][i];

		out << std::endl;
	}
} 

template<int dim, int numchem>
void
Bacteria<dim, numchem>::printInfo(std::ostream& out) const
{

    out << "\n\n-----------------------------------------------------" << std::endl;
    out << "\t\tBACTERIA INFO:";
    out << "\n-----------------------------------------------------" << std::endl;

    out << "\nNumber Bacteria: " << locations.size() << std::endl
    	<< "\nDiffusion constant: " << diffusion_constant << std::endl
    	<< "\nNumber chemicals: " << numchem << std::endl
    	<< "\nIntial secretion rates: " << std::endl;
    for(unsigned int j = 0; j < numchem; j++)
    	out << "\t " << j << "th rate: " << secretion_rates[j][0] << std::endl;

    out << "\n-----------------------------------------------------\n\n" << std::endl;

}

} // close namespace
#endif