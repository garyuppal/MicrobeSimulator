#ifndef MICROBESIMULATOR_REFACTORED_FDM_CHEMICAL_H
#define MICROBESIMULATOR_REFACTORED_FDM_CHEMICAL_H

#include "./fdm_field.h"

#include "../utility/parameter_handler.h"
#include "../geometry/geometry.h" 
#include "../advection/advection_handler.h"

#include "./chemical_interface.h"



namespace MicrobeSimulator{ namespace RefactoredChemicals{

/** \brief Finite difference method implementation of evolving chemical field
* using an upwinded forward euler scheme for dynamics 
*/
template<int dim>
class FDMChemical : public ChemicalInterface<dim>{
public:
	FDMChemical(const ParameterHandler& prm, 
				const Geometry<dim>& geo, 
				const Velocity::AdvectionHandler<dim>& velocity_function,
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

	// no secretion:
	void update() override;

	// secretion, but no control:
	void update(const std::vector<Point<dim> >& locations, 
						const std::vector<double>& amounts) override;

	// secretion and function: (could do zero amount to do just function)
	void update(const std::vector<Point<dim> >& locations, 
						const std::vector<double>& amounts,
						const Function<dim>& control_function) override;


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

private:
	FDM_Field<dim> chemical;
	FDM_Field<dim> auxiliary; // for updating using Forward Euler Method

	FDM_Field<dim, Tensor<1,dim> > flow_field; // for updating using Forward Euler Method

	void setFlowField(const Velocity::AdvectionHandler<dim>& velocity_function);

	// inherited constants:
	// double diffusion_constant;
	// double decay_constant;
	// double time_step;
	// unsigned int chemical_id;

	void updateDiffusion();
	void updateAdvection();
	void updateSources(const std::vector<Point<dim> >& locations, 
						const std::vector<double>& amounts); 
	void updateControls(const Function<dim>& control_function); 

	void addAux();
};

//IMPL:
//--------------------------------------------------------------------------------

/** \brief Constructor for finite difference method implemented chemical field */
template<int dim>
FDMChemical<dim>::FDMChemical(const ParameterHandler& prm, 
							const Geometry<dim>& geo, 
							const Velocity::AdvectionHandler<dim>& velocity_function,
							double dt, 
							unsigned int id) 
	:
	ChemicalInterface<dim>(prm, dt, id)
{
	// get bounds:
	const Point<dim> lower = geo.getBottomLeftPoint();
	const Point<dim> upper = geo.getTopRightPoint();
	
	// boundary conditions:
	const std::array<BoundaryCondition, dim> bcs = geo.getBoundaryConditions();
	
	// intialize with zeros:
	const double init_value = 0.;

	// get discretization:
	std::array<unsigned int, dim> disc;
	/** @todo this only works for 2 dim, need to generalize parameter handler or
	* declare parameters in a different manner */
	std::vector<Point<2> > res = prm.get_point_list("Chemicals", "FDM discretization"); // gives resolution

	for(unsigned int i = 0; i < dim; ++i)
	{
		double dx = res[id][i];
		disc[i] = std::ceil( (upper[i]-lower[i])/dx );
	}

	// intialize chemicals and auxillary fields:
	chemical.init(lower, upper, disc, bcs, init_value);
	auxiliary.init(lower, upper, disc, bcs, init_value);
	flow_field.init(lower, upper, disc, bcs, Tensor<1, dim>() );

	// set local flow field:
	setFlowField(velocity_function);
}

/** \brief Set local flow field for future advection updates (2D) */
template<>
void 
FDMChemical<2>::setFlowField(const Velocity::AdvectionHandler<2>& velocity_function)
{
	for(unsigned int i = 0; i < flow_field.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < flow_field.getDiscretization(1); ++j)
			flow_field.at(i,j) = velocity_function.value( flow_field.getCellCenterPoint(i,j) );
}

/** \brief Set local flow field for future advection updates (3D) */
template<>
void 
FDMChemical<3>::setFlowField(const Velocity::AdvectionHandler<3>& velocity_function)
{
	for(unsigned int i = 0; i < flow_field.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < flow_field.getDiscretization(1); ++j)
			for(unsigned int k = 0; k < flow_field.getDiscretization(2); ++k)
				flow_field.at(i,j,k) = velocity_function.value( flow_field.getCellCenterPoint(i,j,k) );
}

// interface functions:

/** \brief Return value of field at given point */
template<int dim>
double 
FDMChemical<dim>::value(const Point<dim>& p) const
{
	return chemical.at(p);
}

/** \brief Return vector of values at given points */
template<int dim>
std::vector<double> 
FDMChemical<dim>::value_list(const std::vector<Point<dim> >& points) const
{
	const unsigned int n_pts = points.size();
	std::vector<double> result;
	result.reserve(n_pts);

	for(unsigned int i = 0; i < n_pts; ++i)
		result.emplace_back(this->value(points[i]));

	return result;
}

/** \brief Project continuous function onto discrete field (2D) */
/** Replaces current chemical field with discrete representation of given function */
template<>
void 
FDMChemical<2>::project_function(const Function<2>& function)
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < chemical.getDiscretization(1); ++j)
			chemical.at(i,j) = function.value( chemical.getCellCenterPoint(i,j) );
}

/** \brief Project continuous function onto discrete field (3D) */
/** Replaces current chemical field with discrete representation of given function */
template<>
void 
FDMChemical<3>::project_function(const Function<3>& function)
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < chemical.getDiscretization(1); ++j)
			for(unsigned int k = 0; k < chemical.getDiscretization(2); ++k)
				chemical.at(i,j,k) = function.value( chemical.getCellCenterPoint(i,j,k) );
}

// updating methods:
// -------------------------------------------------------------

/** \brief Two dimensional update of diffusion */
/** Replaces current aux field with laplacian of current chemical field */
template<>
void 
FDMChemical<2>::updateDiffusion()
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
FDMChemical<3>::updateDiffusion()
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

/** \brief Update auxillary field to include advection effect (2D) */
/** Adds to current aux field */
template<>
void 
FDMChemical<2>::updateAdvection() // could pass in advection and just ignore for fem, dg impl's
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
FDMChemical<3>::updateAdvection() 
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
FDMChemical<dim>::updateSources(const std::vector<Point<dim> >& locations, 
						const std::vector<double>& amounts)
{
	double idv = 1.; // divide by cell volume to scale with resolution

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		idv *= auxiliary.getInverseCellWidth(dim_itr);

	for(unsigned int i = 0; i < locations.size(); ++i)
		auxiliary.at(locations[i]) += idv*amounts[i]; 
}

/** \brief Update auxiliary field by control function (2D) */
/** Adds to current aux field */
template<>
void 
FDMChemical<2>::updateControls(const Function<2>& control_function)
{
	for(unsigned int i = 0; i < auxiliary.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < auxiliary.getDiscretization(1); ++j)
			auxiliary.at(i,j) += control_function.value( auxiliary.getCellCenterPoint(i,j) );
}

/** \brief Update auxiliary field by control function (3D) */
/** Adds to current aux field */
template<>
void 
FDMChemical<3>::updateControls(const Function<3>& control_function)
{
	for(unsigned int i = 0; i < auxiliary.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < auxiliary.getDiscretization(1); ++j)
			for(unsigned int k = 0; k < auxiliary.getDiscretization(2); ++k)
				auxiliary.at(i,j,k) += control_function.value( auxiliary.getCellCenterPoint(i,j,k) );
}

/** \brief Two dimensional implementation of internal update */
template<>
void 
FDMChemical<2>::addAux()
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < chemical.getDiscretization(1); ++j)
			chemical.at(i,j) = (1.0-time_step*decay_constant)*chemical.at(i,j)
				+ (this->time_step)*auxiliary.at(i,j);
}

/** \brief Three dimensional implementation of internal update */
template<>
void 
FDMChemical<3>::addAux()
{
	for(unsigned int i = 0; i < chemical.getDiscretization(0); ++i)
		for(unsigned int j = 0; j < chemical.getDiscretization(1); ++j)
			for(unsigned int k = 0; k < chemical.getDiscretization(2); ++k)
				chemical.at(i,j,k) = (1.0-time_step*decay_constant)*chemical.at(i,j,k) 
					+ (this->time_step)*auxiliary.at(i,j,k);
}

/** \brief Update field by its time step */
template<int dim>
void 
FDMChemical<dim>::update()
{
	updateDiffusion();
	updateAdvection();
	addAux();
}

/** \brief Update field by its time step 
* and add sources to given locations by given amounts */
template<int dim>
void 
FDMChemical<dim>::update(const std::vector<Point<dim> >& locations, 
					const std::vector<double>& amounts)
{
	updateDiffusion();
	updateAdvection();
	updateSources(locations, amounts);
	addAux();
}

/** \brief Update field by its time step,
* add sources to given locations by given amounts 
* and add values of control functions */
template<int dim>
void 
FDMChemical<dim>::update(const std::vector<Point<dim> >& locations, 
					const std::vector<double>& amounts,
					const Function<dim>& control_function)
{
	updateDiffusion();
	updateAdvection();
	updateSources(locations, amounts);
	updateControls(control_function);
	addAux();	
}

/** \brief integrate field over total volume */
template<int dim>
double 
FDMChemical<dim>::getMass() const
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
FDMChemical<dim>::print(std::ostream& out) const
{
	chemical.print(out);
}

/** \brief Output chemical info */
template<int dim>
void 
FDMChemical<dim>::printInfo(std::ostream& out) const
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

}} // CLOSE NAMESPACES
#endif