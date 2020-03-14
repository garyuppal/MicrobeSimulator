#ifndef MICROBESIMULATOR_AS_BACTERIUM_H
#define MICROBESIMULATOR_AS_BACTERIUM_H

#include "./bacterium_base.h"

namespace MicrobeSimulator{ namespace MultiBacteria{ 

template<int dim>
class ASBacterium : public BacteriumBase<dim>
{
public:
	ASBacterium();

private:
	// inherits location
	// inherits public_good_secretion
	// inherits waste_secretion

	// suicide strategy... when around cheaters, thus when 
	// ratio of waste to public good is too large

	double threshold_ratio; // for determining suicide strategy
	// this is essentially just having a larger susceptiblity to waste
};

// IMPLEMENTATION:
// -------------------------------------------------------------
template<int dim>
ASBacterium<dim>::ASBacterium()
	:
	threshold_ratio(0)
{}



}} // close namespace
#endif