#ifndef MICROBESIMULATOR_BACTERIA_PARAMETERS_H
#define MICROBESIMULATOR_BACTERIA_PARAMETERS_H

// for points:
#include <deal.II/base/point.h>

#include "./meta_argparser.h"

namespace MicrobeSimulator{ namespace Parameters{

class BacteriaParameters{
public: 
	BacteriaParameters();

	void parse(const std::string& parameter_file);
	void checkParameters();

private:
    	// bacteria:
    	unsigned int number_bacteria;
    	unsigned int number_groups;
    	double bacteria_diffusion_constant;
    	double waste_secretion_rate;
    	double good_secretion_rate;
    	double mutation_rate;
    	unsigned int initial_cheaters;
    		// deterministic mutation:
    	unsigned int number_mutate;
    	double deterministic_mutation_time;
    	std::vector<dealii::Point<2> > initial_locations;
    	bool field_intialization;
    	double left_subdomain_length;
    	bool adding_new_groups;
    	double reintroduction_period;

    	// fitness:
		double alpha_good;
		double alpha_waste;
		double secretion_cost;
		double good_saturation;
		double waste_saturation;

		double inefficiency_penalty;
};

// IMPLEMENTATION:
// --------------------------------------------------------------------

}} // close namespaces
#endif