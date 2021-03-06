#pragma once

#include "./aging_fdm_field.h"

#include "../utility/parameter_handler.h"

#include "./aging_chemical_interface.h"



namespace MicrobeSimulator{ namespace Aging{

/** \brief Finite difference method implementation of evolving chemical field
* using an upwinded forward euler scheme for dynamics 
*/
template<int dim>
class AgingChemical : public AgingChemicalInterface<dim>{
public:
	AgingChemical(const ParameterHandler& prm, 
				// const Geometry<dim>& geo, 
				// const Velocity::AdvectionHandler<dim>& velocity_function,
				double dt, 
				unsigned int id);

	// static void declare_parameters(ParameterHandler& prm); // chemical handler takes care of this for now

	// interface functions:
	double value(const Point<dim>& p) const override;  
	std::vector<double> 
		value_list(const std::vector<Point<dim> >& points) const override;

	// project continuous function onto discrete field
	void project_function(const Function<dim>& initial_condition) override;

	// updating methods:

	// // no secretion:
	// void update() override;

	// // secretion, but no control:
	// void update(const std::vector<Point<dim> >& locations, 
	// 					const std::vector<double>& amounts) override;

	// // secretion and function: (could do zero amount to do just function)
	// void update(const std::vector<Point<dim> >& locations, 
	// 					const std::vector<double>& amounts,
	// 					const Function<dim>& control_function) override;

	void update(const std::vector<Point<dim> >& source_locations, 
					double sources,
					const std::vector<Point<dim> >& sink_locations, 
					double sinks) override;


	// integrate field over total volume 
	double getMass() const override;

	// output
	void print(std::ostream& out) const override;
	void printInfo(std::ostream& out) const override;

	// inherited accessors:
	// double getDiffusionConstant() const;
	// double getDecayConstant() const;
	// double getTimeStep() const;
	// unsigned int getID() const;

	// void process_solution(dealii::ConvergenceTable& convergence_table, // can replace by just vector of values later
	// 				const Function<dim>& exact_solution,
	// 				const unsigned int cycle) const override;

private:
	Aging::Aging_FDM_Field<dim> chemical;
	Aging::Aging_FDM_Field<dim> auxiliary; // for updating using Forward Euler Method

	Aging::Aging_FDM_Field<dim, Tensor<1,dim> > flow_field; // for updating using Forward Euler Method

	// void setFlowField(const Velocity::AdvectionHandler<dim>& velocity_function);

	// inherited constants:
	// double diffusion_constant;
	// double decay_constant;
	// double time_step;
	// unsigned int chemical_id;

	void updateDiffusion();
	void updateAdvection();
	void updateSources(const std::vector<Point<dim> >& locations, 
						double amounts); 
	void updateSinks(const std::vector<Point<dim> >& locations, 
						double amounts); 
	void updateControls(const Function<dim>& control_function); 

	void addAux();

	void ensure_positive(); 
};

//IMPL:
//--------------------------------------------------------------------------------

/** \brief Constructor for finite difference method implemented chemical field */
template<int dim>
AgingChemical<dim>::AgingChemical(const ParameterHandler& prm, 
							// const Geometry<dim>& geo, 
							// const Velocity::AdvectionHandler<dim>& velocity_function,
							double dt, 
							unsigned int id) 
	:
	AgingChemicalInterface<dim>(prm, dt, id)
{
	const std::string section = "Chemicals";

	// get bounds:
	std::vector<Point<dim> > bounds = prm.get_oneD_point_list(section, "Bounds");
	Point<dim> lower = bounds[0];
	Point<dim> upper = bounds[1];
	
	// boundary conditions:
	std::array<BoundaryCondition, dim> bcs;
	std::vector<std::string> strbcs = prm.get_list(section, "Boundary conditions");

	for(unsigned int i = 0; i < dim; ++i)
		bcs[i] = stringToBoundaryCondition(strbcs[i]);
	
	// intialize with zeros:
	const double init_value = 0.;

	// get discretization:
	std::array<unsigned int, dim> disc;
	std::vector<double> res = prm.get_double_vector(section, "FDM discretization"); // gives resolution

	for(unsigned int i = 0; i < dim; ++i)
	{
		double dx = res[id]; // same for all dimensions but can vary between chemicals
		disc[i] = std::ceil( (upper[i]-lower[i])/dx );
	}

	// intialize chemicals and auxillary fields:
	Tensor<1, dim> flow;
	flow[0] = prm.get_double(section,"Flow rate");
	chemical.init(lower, upper, disc, bcs, init_value);
	auxiliary.init(lower, upper, disc, bcs, init_value);
	flow_field.init(lower, upper, disc, bcs, flow );

	// set local flow field:
	// setFlowField(velocity_function);
}

/** \brief Set local flow field for future advection updates (2D) */
// template<>
// void 
// AgingChemical<2>::setFlowField(const Velocity::AdvectionHandler<2>& velocity_function)
// {
// 	for(unsigned int i = 0; i < flow_field.getDiscretization(0); ++i)
// 		for(unsigned int j = 0; j < flow_field.getDiscretization(1); ++j)
// 			flow_field.at(i,j) = velocity_function.value( flow_field.getCellCenterPoint(i,j) );
// }

// /** \brief Set local flow field for future advection updates (3D) */
// template<>
// void 
// AgingChemical<3>::setFlowField(const Velocity::AdvectionHandler<3>& velocity_function)
// {
// 	for(unsigned int i = 0; i < flow_field.getDiscretization(0); ++i)
// 		for(unsigned int j = 0; j < flow_field.getDiscretization(1); ++j)
// 			for(unsigned int k = 0; k < flow_field.getDiscretization(2); ++k)
// 				flow_field.at(i,j,k) = velocity_function.value( flow_field.getCellCenterPoint(i,j,k) );
// }

// interface functions:

/** \brief Return value of field at given point */
template<int dim>
double 
AgingChemical<dim>::value(const Point<dim>& p) const
{
	return chemical.at(p);
}

/** \brief Return vector of values at given points */
template<int dim>
std::vector<double> 
AgingChemical<dim>::value_list(const std::vector<Point<dim> >& points) const
{
	const unsigned int n_pts = points.size();
	std::vector<double> result;
	result.reserve(n_pts);

	for(unsigned int i = 0; i < n_pts; ++i)
		result.emplace_back(this->value(points[i]));

	return result;
}

/** \brief Project continuous function onto discrete field (1D) */
/** Replaces current chemical field with discrete representation of given function */
template<>
void 
AgingChemical<1>::project_function(const Function<1>& function)
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		chemical.at(i) = function.value( chemical.getCellCenterPoint(i) );
}

/** \brief Project continuous function onto discrete field (2D) */
/** Replaces current chemical field with discrete representation of given function */
template<>
void 
AgingChemical<2>::project_function(const Function<2>& function)
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < chemical.getDiscretization(1); ++j)
			chemical.at(i,j) = function.value( chemical.getCellCenterPoint(i,j) );
}

/** \brief Project continuous function onto discrete field (3D) */
/** Replaces current chemical field with discrete representation of given function */
template<>
void 
AgingChemical<3>::project_function(const Function<3>& function)
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < chemical.getDiscretization(1); ++j)
			for(unsigned int k = 0; k < chemical.getDiscretization(2); ++k)
				chemical.at(i,j,k) = function.value( chemical.getCellCenterPoint(i,j,k) );
}

// updating methods:
// -------------------------------------------------------------

/** \brief One dimensional update of diffusion */
/** Replaces current aux field with laplacian of current chemical field */
template<>
void 
AgingChemical<1>::updateDiffusion()
{
	const double inv_dx_sqr = chemical.getInverseCellWidthSquared(0);

	for(int i = 0; i < (int)chemical.getDiscretization(0); ++i)
		auxiliary.at(i) = diffusion_constant*(
			inv_dx_sqr*(chemical.at_bc(i+1) - 2*chemical.at(i)
				+chemical.at_bc(i-1))
			);
}

/** \brief Two dimensional update of diffusion */
/** Replaces current aux field with laplacian of current chemical field */
template<>
void 
AgingChemical<2>::updateDiffusion()
{
	const double inv_dx_sqr = chemical.getInverseCellWidthSquared(0);
	const double inv_dy_sqr = chemical.getInverseCellWidthSquared(1);

	for(int i = 0; i < (int)chemical.getDiscretization(0); ++i)
		for(int j = 0; j < (int)chemical.getDiscretization(1); ++j)
			auxiliary.at(i,j) = diffusion_constant*(
				inv_dx_sqr*(chemical.at_bc(i+1,j) - 2*chemical.at(i,j)
					+chemical.at_bc(i-1,j))
				+inv_dy_sqr*(chemical.at_bc(i,j+1) - 2*chemical.at(i,j)
					+chemical.at_bc(i,j-1))
				);
}

/** \brief Three dimensional update of diffusion */
/** Replaces current aux field with laplacian of current chemical field */
template<>
void 
AgingChemical<3>::updateDiffusion()
{
	const double inv_dx_sqr = chemical.getInverseCellWidthSquared(0);
	const double inv_dy_sqr = chemical.getInverseCellWidthSquared(1);
	const double inv_dz_sqr = chemical.getInverseCellWidthSquared(2);

	for(int i = 0; i < (int)chemical.getDiscretization(0); ++i)
		for(int j = 0; j < (int)chemical.getDiscretization(1); ++j)
			for(int k = 0; k < (int)chemical.getDiscretization(2); ++k)
				auxiliary.at(i,j,k) = diffusion_constant*(
					inv_dx_sqr*( chemical.at_bc(i+1,j,k) - 2*chemical.at(i,j,k)
						+chemical.at_bc(i-1,j,k) )
					+inv_dy_sqr*( chemical.at_bc(i,j+1,k) - 2*chemical.at(i,j,k)
						+chemical.at_bc(i,j-1,k) )
					+inv_dz_sqr*( chemical.at_bc(i,j,k+1) - 2*chemical.at(i,j,k)
						+chemical.at_bc(i,j,k-1) )
					);
}


/** \brief Update auxillary field to include advection effect (1D) */
/** Adds to current aux field */
template<>
void 
AgingChemical<1>::updateAdvection() 
{
	const double inv_dx = auxiliary.getInverseCellWidth(0);

	for(int i = 0; i < (int)auxiliary.getDiscretization(0); ++i)
	{
		const Tensor<1, 1> velocity = flow_field.at(i);

		double v_grad_x = 0;

		// x
		if( velocity[0] > 0 )
			v_grad_x = velocity[0]*inv_dx*(chemical.at_bc(i) - chemical.at_bc(i-1));
		else
			v_grad_x = velocity[0]*inv_dx*(chemical.at_bc(i+1) - chemical.at_bc(i));

		auxiliary.at(i) += - v_grad_x; 
	}
}	

/** \brief Update auxillary field to include advection effect (2D) */
/** Adds to current aux field */
template<>
void 
AgingChemical<2>::updateAdvection() // could pass in advection and just ignore for fem, dg impl's
{
	const double inv_dx = auxiliary.getInverseCellWidth(0);
	const double inv_dy = auxiliary.getInverseCellWidth(1);

	for(int i = 0; i < (int)auxiliary.getDiscretization(0); ++i)
		for(int j = 0; j < (int)auxiliary.getDiscretization(1); ++j)
		{
			const Tensor<1, 2> velocity = flow_field.at(i,j);

			double v_grad_x = 0;
			double v_grad_y = 0;

			// x
			if( velocity[0] > 0 )
				v_grad_x = velocity[0]*inv_dx*(chemical.at_bc(i,j) - chemical.at_bc(i-1,j));
			else
				v_grad_x = velocity[0]*inv_dx*(chemical.at_bc(i+1,j) - chemical.at_bc(i,j));

			// y
			if( velocity[1] > 0 )
				v_grad_y = velocity[1]*inv_dy*(chemical.at_bc(i,j) - chemical.at_bc(i,j-1));
			else
				v_grad_y = velocity[1]*inv_dy*(chemical.at_bc(i,j+1) - chemical.at_bc(i,j));

			auxiliary.at(i,j) += - v_grad_x - v_grad_y; 
		}
}	

/** \brief Update auxillary field to include advection effect (3D) */
/** Adds to current aux field */
template<>
void 
AgingChemical<3>::updateAdvection() 
{
	const double inv_dx = auxiliary.getInverseCellWidth(0);
	const double inv_dy = auxiliary.getInverseCellWidth(1);
	const double inv_dz = auxiliary.getInverseCellWidth(2);

	for(int i = 0; i < (int)auxiliary.getDiscretization(0); ++i) // cast as int to avoid warnings, should be safe as long as discretization is not too large
		for(int j = 0; j < (int)auxiliary.getDiscretization(1); ++j)
			for(int k = 0; k < (int)auxiliary.getDiscretization(2); ++k)
			{
				const Tensor<1, 3> velocity = flow_field.at(i,j,k);

				double v_grad_x = 0;
				double v_grad_y = 0;
				double v_grad_z = 0;

				// x
				if( velocity[0] > 0 )
					v_grad_x = velocity[0]*inv_dx*(chemical.at_bc(i,j,k) - chemical.at_bc(i-1,j,k));
				else
					v_grad_x = velocity[0]*inv_dx*(chemical.at_bc(i+1,j,k) - chemical.at_bc(i,j,k));

				// y
				if( velocity[1] > 0 )
					v_grad_y = velocity[1]*inv_dy*(chemical.at_bc(i,j,k) - chemical.at_bc(i,j-1,k));
				else
					v_grad_y = velocity[1]*inv_dy*(chemical.at_bc(i,j+1,k) - chemical.at_bc(i,j,k));

				// z
				if( velocity[2] > 0 )
					v_grad_z = velocity[2]*inv_dz*(chemical.at_bc(i,j,k) - chemical.at_bc(i,j,k-1));
				else
					v_grad_z = velocity[2]*inv_dz*(chemical.at_bc(i,j,k+1) - chemical.at_bc(i,j,k));

				auxiliary.at(i,j,k) += - v_grad_x - v_grad_y - v_grad_z; 
			}
}

/** \brief Update auxiliary field by sources */
/** Adds to current aux field */
template<int dim>
void 
AgingChemical<dim>::updateSources(const std::vector<Point<dim> >& locations, 
						double amounts)
{
	double idv = 1.; // divide by cell volume to scale with resolution

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		idv *= auxiliary.getInverseCellWidth(dim_itr);

	for(unsigned int i = 0; i < locations.size(); ++i)
		auxiliary.at(locations[i]) += idv*amounts; 
}

template<int dim>
void 
AgingChemical<dim>::updateSinks(const std::vector<Point<dim> >& locations, 
						double amounts)
{
	double idv = 1.; // divide by cell volume to get density

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		idv *= auxiliary.getInverseCellWidth(dim_itr);

	// subtract from auxilliary consumnption amount times local concentration
	for(unsigned int i = 0; i < locations.size(); ++i)
		auxiliary.at(locations[i]) += -idv*amounts*(chemical.at(locations[i]));
									// {  gamma * n  }* {   phi[i]   } 
}

/** \brief Update auxiliary field by control function (2D) */
/** Adds to current aux field */
template<>
void 
AgingChemical<2>::updateControls(const Function<2>& control_function)
{
	for(unsigned int i = 0; i < auxiliary.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < auxiliary.getDiscretization(1); ++j)
			auxiliary.at(i,j) += control_function.value( auxiliary.getCellCenterPoint(i,j) );
}

/** \brief Update auxiliary field by control function (3D) */
/** Adds to current aux field */
template<>
void 
AgingChemical<3>::updateControls(const Function<3>& control_function)
{
	for(unsigned int i = 0; i < auxiliary.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < auxiliary.getDiscretization(1); ++j)
			for(unsigned int k = 0; k < auxiliary.getDiscretization(2); ++k)
				auxiliary.at(i,j,k) += control_function.value( auxiliary.getCellCenterPoint(i,j,k) );
}

/** \brief One dimensional implementation of internal update */
template<>
void 
AgingChemical<1>::addAux()
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		chemical.at(i) = (1.0-time_step*decay_constant)*chemical.at(i)
			+ (this->time_step)*auxiliary.at(i);
}

/** \brief Two dimensional implementation of internal update */
template<>
void 
AgingChemical<2>::addAux()
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < chemical.getDiscretization(1); ++j)
			chemical.at(i,j) = (1.0-time_step*decay_constant)*chemical.at(i,j)
				+ (this->time_step)*auxiliary.at(i,j);
}

/** \brief Three dimensional implementation of internal update */
template<>
void 
AgingChemical<3>::addAux()
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < chemical.getDiscretization(1); ++j)
			for(unsigned int k = 0; k < chemical.getDiscretization(2); ++k)
				chemical.at(i,j,k) = (1.0-time_step*decay_constant)*chemical.at(i,j,k) 
					+ (this->time_step)*auxiliary.at(i,j,k);
}

template<int dim>
void
AgingChemical<dim>::ensure_positive()
{
	const unsigned int n = chemical.getSize();
	
	for(unsigned int i = 0; i < n; ++i)
		chemical[i] = (chemical[i]<0)? 0 : chemical[i];
}

// /** \brief Update field by its time step */
// template<int dim>
// void 
// AgingChemical<dim>::update()
// {
// 	updateDiffusion();
// 	updateAdvection();
// 	addAux();
// }

// /** \brief Update field by its time step 
// * and add sources to given locations by given amounts */
// template<int dim>
// void 
// AgingChemical<dim>::update(const std::vector<Point<dim> >& locations, 
// 					const std::vector<double>& amounts)
// {
// 	updateDiffusion();
// 	updateAdvection();
// 	updateSources(locations, amounts);
// 	addAux();
// }

// /** \brief Update field by its time step,
// * add sources to given locations by given amounts 
// * and add values of control functions */
// template<int dim>
// void 
// AgingChemical<dim>::update(const std::vector<Point<dim> >& locations, 
// 					const std::vector<double>& amounts,
// 					const Function<dim>& control_function)
// {
// 	updateDiffusion();
// 	updateAdvection();
// 	updateSources(locations, amounts);
// 	updateControls(control_function);
// 	addAux();	
// }

template<int dim>
void 
AgingChemical<dim>::update(const std::vector<Point<dim> >& source_locations, 
				double sources,
				const std::vector<Point<dim> >& sink_locations, 
				double sinks)
{
	updateDiffusion();
	updateAdvection();
	updateSources(source_locations, sources);
	updateSinks(sink_locations, sinks);
	// updateControls(control_function);
	addAux();	

	ensure_positive(); // in case consume too much
}




/** \brief integrate field over total volume */
template<int dim>
double 
AgingChemical<dim>::getMass() const
{
	double dv = 1.0; 
	for(unsigned int i = 0; i < dim; ++i)
		dv *= chemical.getCellWidth(i);

	double total_mass = 0.;

	// iterate over all elements of chemical:
	for(unsigned int i = 0; i < chemical.getSize(); ++i)
		total_mass += chemical[i]; 

	return total_mass*dv;
}

/** \brief Output field to ostream */
template<int dim>
void 
AgingChemical<dim>::print(std::ostream& out) const
{
	chemical.print(out);
}

/** \brief Output chemical info */
template<int dim>
void 
AgingChemical<dim>::printInfo(std::ostream& out) const
{
	out << "\n\n" << Utility::medium_line << std::endl
		<< "\t\t FINITE DIFFERENCE CHEMICAL: " << std::endl
		<<  Utility::medium_line << std::endl
		<< "Diffusion constant: " << this->diffusion_constant << std::endl
		<< "Decay constant: " << this->decay_constant << std::endl
		<< "Time step: " << this->time_step << std::endl
		<< std::endl << Utility::medium_line 
		<< std::endl << std::endl << std::endl;	
 	/** @todo print out field info too */
			// chemical.printInfo(out);
}

// template<int dim>
// void 
// AgingChemical<dim>::process_solution(dealii::ConvergenceTable& /* convergence_table */, 
// 				const Function<dim>& /* exact_solution */,
// 				const unsigned int /* cycle */) const
// {
// 	std::cout << "Not implemented" << std::endl;
// 	assert(false);
// }

}} // CLOSE NAMESPACES
/* aging_chemical.h */