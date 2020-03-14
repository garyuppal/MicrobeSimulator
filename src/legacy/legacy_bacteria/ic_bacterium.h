#ifndef MICROBESIMULATOR_IC_BACTERIUM_H
#define MICROBESIMULATOR_IC_BACTERIUM_H

#include "./bacterium_base.h"
#include "../utility/bool_functions.h"

namespace MicrobeSimulator{ namespace MultiBacteria{

template<int dim>
class ICBacterium : public BacteriumBase<dim>
{
public:
	ICBacterium();
	ICBacterium(double pg, double w);
	ICBacterium(const Point<dim>& p, double pg, double w);

	// intermittent cheating functions:
	void updateTimeAndSecretion(double dt, 
		const BoolFunctions::BoolFunction& switch_function);

	void setFullSecretionRate(double pg);
	void setTime(double t);

	double getTime() const;
	double getFullSecretionRate() const;

	void print(std::ostream& out) const override;

private:
	// inherits location
	// inherits public_good_secretion -- use this for fitness, current rate
	// inherits waste_secretion

	double time; // internal clock for switching strategies
	double full_secretion_rate; // this is specifically for a 2pg bacterium...
	// double switching_rate; 

	// unsigned int number_switches;

	void update_secretion(const BoolFunctions::BoolFunction& switch_function);
};

// IMPLEMENTATION:
// --------------------------------------------------------------------------------
template<int dim>
ICBacterium<dim>::ICBacterium()
	:
	time(0),
	full_secretion_rate(0)
	// switching_rate(0),
	// number_switches(0)
{}

template<int dim>
ICBacterium<dim>::ICBacterium(double pg, double w)
	:
	BacteriumBase<dim>(pg,w),
	time(0),
	full_secretion_rate(pg)
	// switching_rate(sr),
	// number_switches(0)
{}

template<int dim>
ICBacterium<dim>::ICBacterium(const Point<dim>& p, double pg, double w)
	:
	BacteriumBase<dim>(p,pg,w),
	time(0),
	full_secretion_rate(pg)
	// switching_rate(sr),
	// number_switches(0)
{}


// intermittent cheating functions:
template<int dim>
void 
ICBacterium<dim>::updateTimeAndSecretion(double dt,
	const BoolFunctions::BoolFunction& switch_function)
{
	time += dt;
	update_secretion(switch_function);
}

template<int dim>
void 
ICBacterium<dim>::update_secretion(const BoolFunctions::BoolFunction& switch_function)
{
	if( switch_function.value(time) )
		this->public_good_secretion = full_secretion_rate;
	else
		this->public_good_secretion = 0;

	// if( std::floor(time/switching_rate) > number_switches )
	// {
	// 	// std::cout << "updating IC secretion!!" << std::endl;
	// 	if(this->public_good_secretion != 0)
	// 		this->public_good_secretion = 0;
	// 	else
	// 		this->public_good_secretion = full_secretion_rate;

	// 	++number_switches;
	// }
}

template<int dim>
void 
ICBacterium<dim>::setFullSecretionRate(double pg)
{
	full_secretion_rate = pg;
}

// template<int dim>
// void 
// ICBacterium<dim>::setSwitchingRate(double sr)
// {
// 	switching_rate = sr;
// }

template<int dim>
void 
ICBacterium<dim>::setTime(double t)
{
	time = t;
}

template<int dim>
double 
ICBacterium<dim>::getTime() const
{
	return time;
}

template<int dim>
double 
ICBacterium<dim>::getFullSecretionRate() const
{
	return full_secretion_rate;
}

// template<int dim>
// double 
// ICBacterium<dim>::getSwitchingRate() const
// {
// 	return switching_rate;
// }

// template<int dim>
// unsigned int 
// ICBacterium<dim>::getNumberSwitches() const
// {
// 	return number_switches;
// }

template<int dim>
void 
ICBacterium<dim>::print(std::ostream& out) const
{
    out << this->location << " "
    	<< this->public_good_secretion << " "
    	<< this->waste_secretion << " "
    	<< time << " "
    	<< full_secretion_rate << std::endl;
    	// " "
    	// << switching_rate << " "
    	// << number_switches << std::endl;
}


}} // close namespaces
#endif