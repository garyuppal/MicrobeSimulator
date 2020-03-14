#ifndef MICROBESIMULATOR_AGING_CELLS_H
#define MICROBESIMULATOR_AGING_CELLS_H

#include "./bacterium_base.h"
#include <vector>

namespace MicrobeSimulator{ namespace MultiBacteria{

// // perhaps this is a design flaw and we should
// // restruture things, but for now, delete randomstep:
// void randomStep(double time_step, double diffusion_constant,
// 	const Geometry<dim>& geometry, 
// 	const Velocity::AdvectionHandler<dim>& velocity) = delete;

// // similarly with waste chemicals .. this class only uses one chemical!!


// also we want an array of cells ... bacterium base will do fine,
	// it has more funcionality than we need ... and fitness will have to be
	// for 2 chem, even though we interface with only one chemical...

template<int dim>
class AgingCells{
public:
	AgingCells();

	void init(unsigned int n_cells, double sr, // use the public good rate of bacterium base 
		const Point<dim>& bottom_left, const Point<dim>& top_right);

	// can non-randomly, but uniformly setup cells:
	void init(double sr, // use the public good rate of bacterium base 
		const Point<dim>& bottom_left, const Point<dim>& top_right,
		double density, double resolution);

	void init(double sr, const std::vector<Point<dim> >& locations);


	// the following uses two chems because of bacteria base... should redesign this eventually...
	void random_death(double dt, const FitnessBase<dim, 2>& fitness_function); 	

	// for chemical interface:
	std::vector<Point<dim> > 	getLocations() const; // will change as cells die
	double 						getSecretionRate() const; // use base::pg -- does not mutate
	std::vector<double>			getSecretionRateVector() const; 
										// also store locally
	bool 						isAlive() const;

	void 						print(std::ostream& out) const;
	void 						printInfo(std::ostream& out) const;
private:
	std::vector<BacteriumBase<dim> > cells;

	double 	secretion_rate;

	// this method could be defined outside of this class ...
	Point<dim> getRandomLocationInDomain(const Point<dim>& bottom_left,
		const Point<dim>& top_right) const;
};

// IMPLEMENTATION:
// ------------------------------------------------------------------------------
template<int dim>
AgingCells<dim>::AgingCells()
	:
	secretion_rate(0)
{}

template<int dim>
void 
AgingCells<dim>::init(unsigned int n_cells, double sr, 
	const Point<dim>& bottom_left, const Point<dim>& top_right)
{
	cells.clear();
	cells.reserve(n_cells);
	secretion_rate = sr;
	
	for(unsigned int i = 0; i < n_cells; ++i)
	{
		Point<dim> location = getRandomLocationInDomain(bottom_left, top_right);

		cells.emplace_back( BacteriumBase<dim>(location, 0, 0) ); 
			// NOTE: this suggests a more basic `` particles '' base class
	}
}

// setup cells uniformly with proper (rescaled) secretion rate based on coarse grained area
template<>
void 
AgingCells<2>::init(double sr, 
	const Point<2>& bottom_left, const Point<2>& top_right,
	double density, double resolution) 
{
	cells.clear();

	const double density_to_mass_factor = resolution*resolution; 
	// unit resolution implies single cell per area

	double total_volume = 1;
	unsigned int n_cells = 1;
	for(unsigned int i = 0; i < 2; ++i)
	{
		const double width = top_right[i]-bottom_left[i];
		total_volume *= width;
		n_cells *= (unsigned int)(density*width/resolution);
	}
			
	cells.reserve(n_cells);

	// from chemical:
	// const double scale_factor = 0.0016;
	secretion_rate = 625*(density_to_mass_factor*sr);
	// 		/scale_factor;

	for(double y = 0.5*resolution + bottom_left[1]; y < top_right[1]; y+= resolution)
		for(double x = 0.5*resolution + bottom_left[0]; x < top_right[0]; x += resolution)
			cells.emplace_back( BacteriumBase<2>( Point<2>(x,y) , 0, 0) );
}

// for 3 dimensions:
// template<>
// void 
// AgingCells<3>::init(double secretion_rate, 
// 	const Point<3>& bottom_left, const Point<3>& top_right,
// 	double density, double coarse_grained_volume) 
// {
// 	const double density_to_mass_factor = coarse_grained_volume; // times scale?
// 	const double resolution = std::pow(coarse_grained_volume, 1./3. );

// 	for(double z = bottom_left[2]; z < top_right[2]; z += resolution)
// 		for(double y = bottom_left[1]; y < top_right[1]; y += resolution)
// 			for(double x = bottom_left[0]; x < top_right[0]; x += resolution)
// 				cells.emplace_back( BacteriumBase<3>( Point<3>(x,y,z) , 0, 0) );
// }

template<int dim>
void
AgingCells<dim>::init(double sr, const std::vector<Point<dim> >& locations)
{
	secretion_rate = 625*sr;

	const unsigned int n = locations.size();
	cells.clear();
	cells.reserve(n);

	for(unsigned int i = 0; i < n; ++i)
		cells.emplace_back( BacteriumBase<dim>(locations[i],0,0) );
}


template<int dim>
Point<dim> 
AgingCells<dim>::getRandomLocationInDomain(const Point<dim>& bottom_left,
		const Point<dim>& top_right) const
{
	Point<dim> temp_point;

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		temp_point[dim_itr] = (top_right[dim_itr]-bottom_left[dim_itr])
			*((double) rand() / (RAND_MAX)) + bottom_left[dim_itr]; 

	return temp_point;
}


// the following uses two chems because of bacteria base... should redesign this eventually...
template<int dim>
void 
AgingCells<dim>::random_death(double dt, const FitnessBase<dim, 2>& fitness_function)
{
	for(auto it = cells.begin(); 
			it != cells.end(); )
	{
		const double p_death = it->getFitness(fitness_function)*dt;
		double prob = ((double) rand() / (RAND_MAX));
		if(prob < p_death)
			it = cells.erase(it);
		else
			++it;
	}
}

// for chemical interface:
template<int dim>
std::vector<Point<dim> > 	
AgingCells<dim>::getLocations() const
{
	std::vector<Point<dim> > locations;
	locations.reserve(cells.size());

	for(unsigned int i = 0; i < cells.size(); ++i)
		locations.emplace_back(cells[i].getLocation());

	return locations;
}

template<int dim>
double 		
AgingCells<dim>::getSecretionRate() const
{
	return secretion_rate;
}

template<int dim>
std::vector<double>			
AgingCells<dim>::getSecretionRateVector() const
{
	return std::vector<double>(cells.size(), secretion_rate);
}

						

template<int dim>						
bool 						
AgingCells<dim>::isAlive() const
{
	return !cells.empty();
}

template<int dim>
void 						
AgingCells<dim>::print(std::ostream& out) const
{
	for(unsigned int i = 0; i < cells.size(); ++i)
		out << cells[i].getLocation() << std::endl;
}

template<int dim>
void 						
AgingCells<dim>::printInfo(std::ostream& out) const
{
	out << "\n\n-----------------------------------------------------" << std::endl
		<< "AGING CELLS:" << std::endl
		<< "-----------------------------------------------------" << std::endl
		<< "\t number: " << cells.size() << std::endl
		<< "\t secretion rate: " << secretion_rate 
		<< "\n-----------------------------------------------------\n\n" << std::endl;
}

}} // close namespace
#endif