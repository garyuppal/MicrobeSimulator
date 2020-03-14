#ifndef BACTERIA_H
#define BACTERIA_H


// for testing:
#include <deal.II/base/timer.h>

// #include <deal.II/numerics/vector_tools.h>

#include "single_bacterium.h"

// #include "../discrete_field/FDM_chemical.h"

// in single_bacterium.h:
// #include <deal.II/base/point.h>
// using dealii::Point;
// #include <iostream>
// #include <string>
// #include <sstream>
// #include <array>

namespace MicrobeSimulator{

	template<int dim, int numchem>
	class Bacteria{
	public:
		// constructors:
		Bacteria();
		Bacteria(double db, unsigned int numberBacteria, 
			const std::array<double, numchem>& secretion_rates);

		// Bacteria(double db, unsigned int numberBacteria,
		// 	const std::vector<Point<dim> >& locations, 
		// 	double gs, double ws);

		Bacteria(const Bacteria& b);
	 
		// accessors:
		unsigned int getSize() const;
		Point<dim> getLocation(unsigned int i) const;

		std::array<double, numchem>
			 getSecretionRates(unsigned int bact) const;
		double getSecretionRate(unsigned int bact, unsigned int chem) const;

		// legacy:
		// double getGoodSecretionRate(unsigned int i) const;
		// double getWasteSecretionRate(unsigned int i) const;

		// mutators:
		// void initialize(double db, unsigned int numberBacteria, 
		// 	double gs, double ws, unsigned int number_cheaters = 0); // initialize at origin
		void initialize(double db, unsigned int numberBacteria,
			const std::array<double, numchem>& rates);

		// void initialize(double db, unsigned int numberBacteria,
		// 	const std::vector<Point<dim> >& locations, 
		// 	double gs, double ws, unsigned int number_cheaters = 0); 

		void initialize(double db, unsigned int numberBacteria,
			const std::array<double, numchem>& rates, 
			const std::vector<Point<dim> >& locations);

		// intialize from field:
		void initialize(double db, 
			const std::vector<std::pair<Point<dim>, unsigned int> >& coarse_bacteria,
			const std::array<double, numchem>& rates);

		// void setFitnessConstants(double a1, double a2, double k1, double k2, double b);

		// functions:
		void randomWalk(double time_step,
			const Geometry<dim>& geometry,
			const VelocityInterface<dim>& velocity); 

		// void randomWalk(double timeStep, 
		// 	const Geometry<dim>* const geometry, 
		// 	const AdvectionField<dim>* const advection = NULL);

		// void randomWalk(double timeStep,
		// 	const Geometry<dim>* const geometry,
		// 	const DoFHandler<dim>& stokes_dof,
		// 	const BlockVector<double>& stokes_solution);

		// void randomWalk(double timeStep,
		// 	const Geometry<dim>* const geometry,
		// 	const DoFHandler<dim>& stokes_dof,
		// 	const BlockVector<double>& stokes_solution,
		// 	const PointCellMap<dim>& pcm);
		// legacy:
		// void reproduce(double dt, const DoFHandler<dim>& dofh, 
		// 	const Vector<double>& sol1, const Vector<double>& sol2);

		// void reproduce(double dt, const FDMChemical<dim>& goods,
		// 	const FDMChemical<dim>& waste);

		void reproduce(double dt, const FitnessBase<dim, numchem>& fitness_function);
		void reproduce2(double dt, const FitnessBase<dim, numchem>& fitness_function);

		void remove_fallen_bacteria(double right_edge);
		// void reproduce(Chemicals<dim,NumberGoods> chemicals, fitness function...); 
			// cost will have to be updated internally ... fitness -= costs...
		void mutate(double time_step, 
				double mutation_rate, 
				double deltaSecretion, 
				double original_secretion_rate = 0,
				bool binary_mutation = false);
		 // @todo: write another function --- what about for multiple public goods ? 

		void make_cheaters(unsigned int number_cheaters);


		bool isAlive() const;

		void print(std::ostream &out) const; 
		void printInfo(std::ostream &out) const; 

	private:
		double diffusionConstant;

		std::vector<SingleBacterium<dim, numchem> > bacteria;

		// legacy:
		// double getFitness( const  DoFHandler<dim>& dofh, 
		// 	const Vector<double>& sol1, const Vector<double>& sol2,
		// 	const Point<dim>& location, double sec_rate);

		// double getFitness(const FDMChemical<dim>& goods, const FDMChemical<dim>& waste,
		// 	const Point<dim>& location, double sec_rate);

		// switch to indpendent version:
		// double getFitness(const FitnessBase<dim,numchem>& fitness_function);
	}; // class Bacterium



// IMPLEMENTATION
//-----------------------------------------------------------------------------
	template<int dim, int numchem>
	Bacteria<dim, numchem>::Bacteria() 
	{}


	template<int dim, int numchem>
	Bacteria<dim, numchem>::Bacteria(double db, 
		unsigned int numberBacteria, 
		const std::array<double, numchem>& secretion_rates)
		:
		diffusionConstant(db),
		bacteria(numberBacteria)
	{
		for(unsigned int i = 0; i < numberBacteria; i++)
			bacteria[i].setSecretionRates(secretion_rates);
	}
	
	// template<int dim, int numchem>
	// Bacteria<dim, numchem>::Bacteria(double db, const std::vector<Point<dim> >& locations, 
	// 		double gs, double ws)
	// 	:
	// 	diffusionConstant(db),
	// 	bacteria(locations.size())
	// {
	// 	for(unsigned int i = 0; i < locations.size(); i++)
	// 	{
	// 		bacteria[i].setLocation(locations[i]);
	// 		bacteria[i].setGoodSecretionRate(gs);
	// 		bacteria[i].setWasteSecretionRate(ws);
	// 	}
	// }
		
	template<int dim, int numchem>
	Bacteria<dim, numchem>::Bacteria(const Bacteria& b)
	{
		diffusionConstant = b.diffusionConstant;
		bacteria = b.bacteria; // equal operation for vectors???
	}
	 

	// accessors: 
	template<int dim, int numchem>
	unsigned int 
	Bacteria<dim, numchem>::getSize() const
	{ 
		return bacteria.size();
	}


	template<int dim, int numchem>
	Point<dim> 
	Bacteria<dim, numchem>::getLocation(unsigned int i) const
	{ 
		return bacteria[i].getLocation(); 
	}


	template<int dim, int numchem>
	std::array<double, numchem>
	Bacteria<dim, numchem>::getSecretionRates(unsigned int bact) const
	{ 
		return bacteria[bact].getSecretionRates(); 
	}


	template<int dim, int numchem>
	double
	Bacteria<dim, numchem>::getSecretionRate(unsigned int bact,
			unsigned int chem) const
	{ 
		return bacteria[bact].getSecretionRate(chem); 
	}


	// mutators:

	// legacy:
	// template<int dim, int numchem>
	// void Bacteria<dim, numchem>::initialize(double db, unsigned int numberBacteria, 
	// 	double gs, double ws, unsigned int number_cheaters)
	// {
	// 	bacteria.clear();
	// 	bacteria.reserve(numberBacteria);

	// 	diffusionConstant = db;

	// 	unsigned int i = 0;
	// 	for( ; i < number_cheaters; i++)
	// 		bacteria.push_back( SingleBacterium<dim>(Point<dim>(), 0., ws) );

	// 	for( ; i < numberBacteria; i++)
	// 		bacteria.push_back( SingleBacterium<dim>(Point<dim>() , gs, ws) );
	// }


	// template<int dim, int numchem>
	// void Bacteria<dim, numchem>::initialize(double db, 
	// 	unsigned int numberBacteria,
	// 	const std::vector<Point<dim> >& locations, 
	// 	double gs, double ws, unsigned int number_cheaters)
	// {
	// 	bacteria.clear();
	// 	bacteria.reserve(numberBacteria);

	// 	diffusionConstant = db;

	// 	const unsigned int number_groups = locations.size();

	// 	unsigned int i = 0;
	// 	for( ; i < number_cheaters; i++)
	// 	{
	// 		unsigned int group_index = i % number_groups;	

	// 		bacteria.push_back( SingleBacterium<dim>(locations[group_index], 0., ws));
	// 	} 

	// 	for( ; i < numberBacteria; i++)
	// 	{
	// 		unsigned int group_index = i % number_groups;	

	// 		bacteria.push_back( SingleBacterium<dim>(locations[group_index], gs, ws));
	// 	} // for number of bacteria

	// }


	template<int dim, int numchem>
	void 
	Bacteria<dim, numchem>::initialize(double db, 
		unsigned int numberBacteria,
		const std::array<double, numchem>& rates)
	{
		bacteria.clear();
		bacteria.reserve(numberBacteria);

		diffusionConstant = db;

		for(unsigned int i = 0; i < numberBacteria; i++)
			bacteria.push_back( SingleBacterium<dim, numchem>(Point<dim>() , rates) );
	}


	template<int dim, int numchem>
	void 
	Bacteria<dim, numchem>::initialize(double db,
		unsigned int numberBacteria,
		const std::array<double, numchem>& rates, 
		const std::vector<Point<dim> >& locations)
	{
		bacteria.clear();
		bacteria.reserve(numberBacteria);

		diffusionConstant = db;

		const unsigned int number_groups = locations.size();

		for(unsigned int i=0; i < numberBacteria; i++)
		{
			unsigned int group_index = i % number_groups;	
			bacteria.push_back( SingleBacterium<dim, numchem>(locations[group_index], rates) );
		} // for number of bacteria

	}

	// intialize from field:
	template<int dim, int numchem>
	void 
	Bacteria<dim, numchem>::initialize(double db, 
		const std::vector<std::pair<Point<dim>, unsigned int> >& coarse_bacteria,
		const std::array<double, numchem>& rates)
	{
		bacteria.clear();

		unsigned int number_bacteria = 0;
		const unsigned int n_points = coarse_bacteria.size();
		for(unsigned int i = 0; i < n_points; ++i)
			number_bacteria += coarse_bacteria[i].second;

		std::cout << "... intializing " << number_bacteria << " total bacteria" << std::endl;

		if(number_bacteria == 0)
		{
			std::cout << "No bacteria intialized from continuous solution" << std::endl
				<< "Double check continuous parameters and solution" << std::endl;
			return;
		}

		bacteria.reserve(number_bacteria);

		diffusionConstant = db;

		unsigned int b_id = 0;

		// std::cout << "coarse_bacteria vector size: " << n_points << std::endl;

		for(unsigned int i = 0; i < n_points; ++i)
		{
			const Point<dim> location = coarse_bacteria[i].first;
			const unsigned int n_bact = coarse_bacteria[i].second;
			// std::cout << "assigning " << n_bact << " bacteria to location: " << location << std::endl;
			for(unsigned int b = 0; b < n_bact; ++b)
			{
				++b_id;
				// std::cout << "adding bacteria number: " << b_id << std::endl;
				bacteria.push_back(SingleBacterium<dim, numchem>(location, rates));
			}
		}

		std::cout << "... done intializing bacteria" << std::endl;
	}


	// template<int dim, int numchem>
	// void Bacteria<dim, numchem>::setFitnessConstants(double a1, double a2, 
	// 	double k1, double k2, double b)
	// {
	// 	alpha1 = a1;
	// 	alpha2 = a2;
	// 	saturation_const1 = k1;
	// 	saturation_const2 = k2;
	// 	beta1 = b;
	// }

	// functions:
	template<int dim, int numchem>
	void 
	Bacteria<dim, numchem>::randomWalk(double time_step,
		const Geometry<dim>& geometry,
		const VelocityInterface<dim>& velocity)
	{
		for(unsigned int i = 0; i < bacteria.size(); ++i)
			bacteria[i].randomStep(time_step, diffusionConstant, geometry, velocity);

		if(geometry.getBoundaryConditions()[0] == BoundaryCondition::OPEN)
			remove_fallen_bacteria(geometry.getTopRightPoint()[0]);
	} 

	// legacy:
/*
	template<int dim, int numchem>
	void Bacteria<dim, numchem>::randomWalk(double timeStep, 
			const Geometry<dim>* const geometry, 
			const AdvectionField<dim>* const advection)
	{
		for(unsigned int i = 0; i < bacteria.size(); i++)
			bacteria[i].randomStep(timeStep, diffusionConstant, geometry, advection);

		// for open BC's in x-coordinate:
		if( geometry->getBoundaryConditions()[0] == BoundaryCondition::OPEN )
			remove_fallen_bacteria(geometry->getTopRightPoint()[0]);
	}


	template<int dim, int numchem>
	void 
	Bacteria<dim, numchem>::randomWalk(double timeStep,
			const Geometry<dim>* const geometry,
			const DoFHandler<dim>& stokes_dof,
			const BlockVector<double>& stokes_solution)
	{
		for(unsigned int i = 0; i < bacteria.size(); i++)
			bacteria[i].randomStep(timeStep, diffusionConstant, geometry,
				stokes_dof, stokes_solution);

		// for open BC's in x-coordinate:
		if( geometry->getBoundaryConditions()[0] == BoundaryCondition::OPEN )
			remove_fallen_bacteria(geometry->getTopRightPoint()[0]);
	}


	template<int dim, int numchem>
	void 
	Bacteria<dim, numchem>::randomWalk(double timeStep,
			const Geometry<dim>* const geometry,
			const DoFHandler<dim>& stokes_dof,
			const BlockVector<double>& stokes_solution,
			const PointCellMap<dim>& pcm)
	{
		for(unsigned int i = 0; i < bacteria.size(); i++)
			bacteria[i].randomStep(timeStep, diffusionConstant, geometry,
				stokes_dof, stokes_solution, pcm);

		// for open BC's in x-coordinate:
		if( geometry->getBoundaryConditions()[0] == BoundaryCondition::OPEN )
			remove_fallen_bacteria(geometry->getTopRightPoint()[0]);
	}
*/


	template<int dim, int numchem>
	void Bacteria<dim, numchem>::reproduce(double dt, 
		const FitnessBase<dim, numchem>& fitness_function)
	{
		unsigned int size = bacteria.size();
		std::vector<unsigned int> to_kill, to_rep;

		for(unsigned int i = 0; i < size; i++)
		{
			const double fit = dt*bacteria[i].getFitness(fitness_function); 
			// getFitness(dofh,sol1,sol2,
			// 	bacteria[i].getLocation(),bacteria[i].getGoodSecretionRate());
			if(fit < 0)
				to_kill.push_back(i);
			else
			{
				double prob =  ((double) rand() / (RAND_MAX));
				if(prob<fit)
					to_rep.push_back(i);
			}
		} // find those to kill or reproduce

		// go through cases:
		unsigned int num_to_kill = to_kill.size();
		unsigned int num_to_rep = to_rep.size();

		int new_size = size + (int)num_to_rep - (int)num_to_kill;

		if(num_to_kill <= num_to_rep)
		{
			// new size:
			bacteria.reserve(new_size);

			unsigned int i = 0;
			for( ; i < num_to_kill; i++)
		    	bacteria[ to_kill[i] ] = bacteria[ to_rep[i] ]; 
		    // replace dead with offspring

		    // push back remaining:
		    for(; i < num_to_rep; i++)
		      bacteria.emplace_back( bacteria[ to_rep[i] ] );
		}
		else if( new_size > 0 )
		{
			// first push back offspring: (since deleting changes numbering)
			for(unsigned int i = 0; i < num_to_rep; ++i)
				bacteria.emplace_back( bacteria[ to_rep[i] ] );

			// remove dead:
			for(auto it = bacteria.begin(); it != bacteria.end(); )
			{
				const double fit = dt*(it->getFitness(fitness_function));
				if(fit < 0)
					it = bacteria.erase(it);
				else
					++it;
			}

		}
		else
		{
			bacteria.clear();
		} // empty vector
	}


// try faster:
template<int dim, int numchem>
void 
Bacteria<dim, numchem>::reproduce2(double dt, 
		const FitnessBase<dim, numchem>& fitness_function)
{
	const unsigned int size = bacteria.size();
	unsigned int last_index = size - 1; 
	unsigned int last_living = last_index;

	dealii::Timer timer_fit;
	for(unsigned int i = 0; i < size; ++i)
		const double fit = dt*bacteria[i].getFitness(fitness_function); 
	timer_fit.stop();
	std::cout << "\n--fitness access time: " << timer_fit.cpu_time() 
			<< " s\n" << std::endl; 

	// only need to loop over original n bacteria (not over those moved to end)
	for(unsigned int i = 0; i != last_index; ) 
	{
		const double fit = dt*bacteria[i].getFitness(fitness_function); 
		// decide what to do:
		if(fit < 0)
		{
			// swap to end -- not including those already swapped: 
				//(don't need to copy this to end, just delete)
			if(last_living > last_index) // if new offspring added
			{
				bacteria[i] = bacteria[last_living]; // move last offspring to here
				++i; // since an already reproduced bacteria is replacing this
				--last_living; // since moved dead to end 
			}
			else // should have last index == last_living
			{
				bacteria[i] = bacteria[last_living];
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

				if(last_living < (bacteria.size()-1) ) //if dead already at end
				{
					bacteria[last_living+1] = bacteria[i]; // replace
					++last_living;
				}
				else
				{
					bacteria.emplace_back(bacteria[i]); // add copy to end
					++last_living; // to keep with total size
				}
			}
			++i; // move up if either reproduced or did nothing
		} //else possible reproduce
	} // loop once over bacteria vector

	// erase dead bacteria:
	if(last_living == 0) 
		bacteria.clear();
	else
		bacteria.erase(bacteria.begin() + last_living, bacteria.end()); // check...

} // reproduce2()

	// // legacy (wrong):
	// template<int dim, int numchem>
	// void Bacteria<dim, numchem>::reproduce(double dt, 
	// 	const FitnessBase<dim, numchem>& fitness_function)
	// 	// const  DoFHandler<dim>& dofh, 
	// 	// const Vector<double>& sol1, const Vector<double>& sol2)
	// {
	// 	unsigned int size = bacteria.size();
	// 	std::vector<unsigned int> to_kill, to_rep;

	// 	for(unsigned int i = 0; i < size; i++)
	// 	{
	// 		const double fit = dt*bacteria[i].getFitness(fitness_function); 
	// 		// getFitness(dofh,sol1,sol2,
	// 		// 	bacteria[i].getLocation(),bacteria[i].getGoodSecretionRate());
	// 		if(fit < 0)
	// 			to_kill.push_back(i);
	// 		else
	// 		{
	// 			double prob =  ((double) rand() / (RAND_MAX));
	// 			if(prob<fit)
	// 				to_rep.push_back(i);
	// 		}
	// 	} // find those to kill or reproduce

	// 	// go through cases:
	// 	unsigned int num_to_kill = to_kill.size();
	// 	unsigned int num_to_rep = to_rep.size();

	// 	int new_size = size + num_to_rep - num_to_kill;

	// 	if(num_to_kill <= num_to_rep)
	// 	{
	// 		// new size:
	// 		bacteria.reserve(size + num_to_rep - num_to_kill);

	// 		unsigned int i = 0;
	// 		for( ; i < num_to_kill; i++)
	// 	    	bacteria[ to_kill[i] ] = bacteria[ to_rep[i] ]; 
	// 	    // replace dead with offspring

	// 	    // push back remaining:
	// 	    for(; i < num_to_rep; i++)
	// 	      bacteria.emplace_back( bacteria[ to_rep[i] ] );
	// 	}
	// 	else if( new_size > 0 )
	// 	{
	// 		unsigned int i = 0; 
	// 		auto it = bacteria.begin(); 
	// 		for( ; i < num_to_rep; i++)
	// 		{
	// 		    bacteria[ to_kill[i] ] = bacteria[ to_rep[i] ];
	// 			++i;
	// 			++it;
	// 		}
	// 	    for( ; it != bacteria.end(); )
	// 	 	{

	// 	 	}
	// 	    for(; i < num_to_kill; i++)
	// 	    	// error!
	// 	      bacteria.erase( bacteria.begin() + to_kill[i] ); // could get reverse iterator to end...
	// 	}
	// 	else
	// 	{
	// 		bacteria.clear();
	// 	} // empty vector

	// } // later pass a fitness field object...


	// legacy:
	// template<int dim, int numchem>
	// void Bacteria<dim, numchem>::reproduce(double dt, const FDMChemical<dim>& goods,
	// 	const FDMChemical<dim>& waste)
	// {
	// 	unsigned int size = bacteria.size();
	// 	std::vector<unsigned int> to_kill, to_rep;

	// 	for(unsigned int i = 0; i < size; i++)
	// 	{
	// 		const double fit = dt*getFitness(goods, waste,
	// 			bacteria[i].getLocation(),bacteria[i].getGoodSecretionRate());
	// 		if(fit < 0)
	// 			to_kill.push_back(i);
	// 		else
	// 		{
	// 			double prob =  ((double) rand() / (RAND_MAX));
	// 			if(prob<fit)
	// 				to_rep.push_back(i);
	// 		}
	// 	} // find those to kill or reproduce

	// 	// go through cases:
	// 	unsigned int num_to_kill = to_kill.size();
	// 	unsigned int num_to_rep = to_rep.size();

	// 	int new_size = size + num_to_rep - num_to_kill;

	// 	if(num_to_kill <= num_to_rep)
	// 	{
	// 		// new size:
	// 		bacteria.reserve(size + num_to_rep - num_to_kill);

	// 		unsigned int i = 0;
	// 		for( ; i < num_to_kill; i++)
	// 	    	bacteria[ to_kill[i] ] = bacteria[ to_rep[i] ]; 
	// 	    // replace dead with offspring

	// 	    // push back remaining:
	// 	    for(; i < num_to_rep; i++)
	// 	      bacteria.push_back( bacteria[ to_rep[i] ] );
	// 	}
	// 	else if( new_size > 0 )
	// 	{
	// 	    unsigned int i = 0;
	// 	    for( ; i < num_to_rep; i++)
	// 	      bacteria[ to_kill[i] ] = bacteria[ to_rep[i] ];

	// 	    // delete remaining: // *** could probably speed this part up
	// 	    for(; i < num_to_kill; i++)
	// 	      bacteria.erase( bacteria.begin() + to_kill[i] ); // could get reverse iterator to end...
	// 	}
	// 	else
	// 	{
	// 		bacteria.clear();
	// 	} // empty vector

	// } // later pass a fitness field object...


	template<int dim, int numchem>
	void 
	Bacteria<dim, numchem>::remove_fallen_bacteria(double right_edge)
	{
		const double tolerance = 0.2;

		// std::cout << "right edge at: " << right_edge << std::endl;

		for(auto it = bacteria.begin(); it != bacteria.end();)
		{
			if( std::fabs(it->getLocation()[0] - right_edge) < tolerance )	
				it = bacteria.erase(it);
			else
				++it;
		}

		// for(unsigned int i = 0; i < bacteria.size(); i++)
		// {
		// 	if( std::fabs(bacteria[i].getLocation()[0] - right_edge) < tolerance )
		// 	{
		// 		std::cout << "\t distance from edge: " << std::fabs(bacteria[i].getLocation()[0] - right_edge) << std::endl;
		// 		std::cout << "location at: " << (bacteria.begin() + i)->getLocation() << std::endl;
		// 		bacteria.erase( bacteria.begin() + i );
		// 	}
		// }

		// for (auto it = c.begin(); it != c.end(); ) {
	 //        if (*it % 2 == 0) {
	 //            it = c.erase(it);
	 //        } else {
	 //            ++it;
	 //        }
  //   	}
	}


	// legacy:
	// template<int dim, int numchem>
	// double Bacteria<dim, numchem>::getFitness(const  DoFHandler<dim>& dofh, 
	// 	const Vector<double>& sol1, const Vector<double>& sol2, 
	// 	const Point<dim>& location, double sec_rate)
	// {
	// 	double result = 0;
	// 	double c1, c2;

	// 	c1 = dealii::VectorTools::point_value(dofh, sol1, location);
	// 	c2 = dealii::VectorTools::point_value(dofh, sol2, location);

	// 	result = alpha1*(c1/(c1+saturation_const1)) 
	// 	      - alpha2*(c2/(c2+saturation_const2))
	// 		  - beta1*sec_rate;

	// 	// std::cout << "c1 = " << c1 << std::endl
	// 	// 	<< "c2 = " << c2 << std::endl << std::endl;

	// 	// std::cout << "fitness is: " << result << std::endl << std::endl;
	// 	return result;
	// }


	// template<int dim, int numchem>
	// double Bacteria<dim, numchem>::getFitness(const FDMChemical<dim>& goods, 
	// 	const FDMChemical<dim>& waste,
	// 	const Point<dim>& location, double sec_rate)
	// {
	// 	const double c1 = goods.value(location);
	// 	const double c2 = waste.value(location);

	// 	double result = alpha1*(c1/(c1+saturation_const1)) 
	// 		- alpha2*(c2/(c2+saturation_const2))
	// 		- beta1*sec_rate;

	// 	// std::cout << "Fitness constants are: "
	// 	// 	<< "\t alpha_1: " << alpha1 << std::endl
	// 	// 	<< "\t alpha 2: " << alpha2 << std::endl
	// 	// 	<< "\t beta: " << beta1 << std::endl 
	// 	// 	<< "\t s_good: " << sec_rate << std::endl
	// 	// 	<< "\t k1: " << saturation_const1 << std::endl
	// 	// 	<< "\t k2: " << saturation_const2 << std::endl;

	// 	// std::cout << "c1 = " << c1 << std::endl
	// 	// 	<< "c2 = " << c2 << std::endl << std::endl;

	// 	// std::cout << "fitness is: " << result << std::endl << std::endl;

	// 	return result;
	// }


	template<int dim, int numchem>
	void 
	Bacteria<dim, numchem>::mutate(double time_step, double mutation_rate, 
			double deltaSecretion, double original_secretion_rate, bool binary_mutation)
	{
		for(unsigned int i = 0; i < bacteria.size(); i++)
		{
			// random number between 0 and 1:
			double attempt = ((double) rand() / (RAND_MAX));

			if(attempt < time_step*mutation_rate)
				bacteria[i].mutate(deltaSecretion, original_secretion_rate, binary_mutation);
		}
	}


	template<int dim, int numchem>
	void 
	Bacteria<dim, numchem>::make_cheaters(unsigned int number_cheaters)
	{
		if(number_cheaters > bacteria.size())
			number_cheaters = bacteria.size();

		for(unsigned int b = 0; b < number_cheaters; ++b)
			for(unsigned int c = 0; c < numchem-1; ++c) // keep waste secretion on
				bacteria[b].setSecretionRate(c, 0.);
	}



	template<int dim, int numchem>
	bool 
	Bacteria<dim, numchem>::isAlive() const
	{
		return !bacteria.empty();
	}



	template<int dim, int numchem>
	void Bacteria<dim, numchem>::print(std::ostream &out) const
	{
		for(unsigned int i = 0; i < bacteria.size(); i++)
		{
			bacteria[i].printBacterium(out);
			// out << bacteria[i].getLocation() << " "
			// 	<< bacteria[i].getGoodSecretionRate()
			// 	<< std::endl;
		}
		out << std::endl;
	} 


	template<int dim, int numchem>
	void
	Bacteria<dim, numchem>::printInfo(std::ostream& out) const
	{

	    out << "\n\n-----------------------------------------------------" << std::endl;
	    out << "\t\tBACTERIA INFO:";
	    out << "\n-----------------------------------------------------" << std::endl;

	    out << "\nNumber Bacteria: " << bacteria.size() << std::endl
	    	<< "\nDiffusion constant: " << diffusionConstant << std::endl
	    	<< "\nNumber chemicals: " << numchem << std::endl
	    	<< "\nIntial secretion rates: " << std::endl;
	    for(unsigned int j = 0; j < numchem; j++)
	    	out << "\t " << j << "th rate: " << bacteria[0].getSecretionRate(j) << std::endl;

	    	// << "\nGood secretion rate: " << bacteria[0].getGoodSecretionRate() 
	    	// << "\nWaste secretion rate: " << bacteria[0].getWasteSecretionRate() << std::endl
	    	// << "\nFITNESS: " << std::endl
	    	// << "\t public good benefit: " << alpha1 << std::endl
	    	// << "\t waste harm: " << alpha2 << std::endl
	    	// << "\t goods saturation constant: " << saturation_const1 << std::endl
	    	// << "\t waste saturation constant: " << saturation_const2 << std::endl
	    	// << "\t secretion cost: " << beta1 << std::endl;


        out << "\n-----------------------------------------------------\n\n" << std::endl;

	}
}

#endif  // bacterium.h
