#ifndef MICROBESIMULATOR_REG_BACTERIUM_H
#define MICROBESIMULATOR_REG_BACTERIUM_H

/**
*  as is, regBacterium is no different from base ...
*/
#include "./bacterium_base.h"

namespace MicrobeSimulator{ namespace MultiBacteria{

template<int dim>
class RegBacterium : public BacteriumBase<dim>
{
public:
	RegBacterium();
	RegBacterium(double pg, double w);
	RegBacterium(const Point<dim>& p, double pg, double w);
	
	// void print(std::ostream& out) const override;
};

template<int dim>
RegBacterium<dim>::RegBacterium()
{}

template<int dim>
RegBacterium<dim>::RegBacterium(double pg, double w)
	:
	BacteriumBase<dim>(pg,w)
{}

template<int dim>
RegBacterium<dim>::RegBacterium(const Point<dim>& p, double pg, double w)
	:
	BacteriumBase<dim>(p,pg,w)
{}

// member of base class:

// template<int dim>
// void
// RegBacterium<dim>::print(std::ostream& out) const
// {
//     out << this->location << " "
//     	<< this->public_good_secretion << " "
//     	<< this->waste_secretion << std::endl;
// }

}} //close namespace
#endif