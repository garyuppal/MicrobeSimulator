#ifndef MICROBESIMULATOR_FDM_CHEMICAL_H
#define MICROBESIMULATOR_FDM_CHEMICAL_H


#include "../advection/advection.h"
#include "./discrete_field.h"

namespace MicrobeSimulator{

	template<int dim>
	class FDMChemical{
	public:
		FDMChemical();
		FDMChemical(const Point<dim>& lower, const Point<dim>& upper,
			const std::array<unsigned int, dim>& disc);
		FDMChemical(const FDMChemical& chem);

		void initialize(const Point<dim>& lower, const Point<dim>& upper,
			const std::array<unsigned int, dim>& disc,
			const std::array<BoundaryCondition, dim>& bcs);

		// access
		double value(const Point<dim>& p) const;
		double& at(const Point<dim> &p);

		double getDiffusionConstant() const;
		double getDecayConstant() const;

		// mutators:
		void setDiffusionConstant(double diff);
		void setDecayConstant(double decay);

		void projectFunction(const DiscreteFunction<dim>& function);

		// evolve:
		void secreteTo(const Point<dim> p, double amount);
		void update(double time_step, const AdvectionField<dim>* const advection = NULL);

		void print(std::ostream& out) const;
		void printAux(std::ostream& out) const;

		void printInfo(std::ostream& out) const;

	private:
		double diffusion_constant;
		double decay_constant;

		DiscreteField<dim> chemical;
		DiscreteField<dim> auxiliary;

		void updateDiffusion();
		void updateAdvection( const AdvectionField<dim>* const advection);
	};

//IMPLEMENTATION:
//--------------------------------------------------------------------------------

	template<int dim>
	FDMChemical<dim>::FDMChemical()
	{}

	template<int dim>
	FDMChemical<dim>::FDMChemical(const Point<dim>& lower, const Point<dim>& upper,
		const std::array<unsigned int, dim>& disc)
		:
		chemical(lower,upper,disc),
		auxiliary(lower,upper,disc)
	{}

	template<int dim>
	FDMChemical<dim>::FDMChemical(const FDMChemical& chem)
	{
		diffusion_constant = chem.diffusion_constant;
		decay_constant = chem.decay_constant;
		chemical = chem.chemical; // does this work?? do we need to overload = for field???
		auxiliary = chem.auxiliary;
	}

	template<int dim>
	void FDMChemical<dim>::initialize(const Point<dim>& lower, const Point<dim>& upper,
		const std::array<unsigned int, dim>& disc, 
		const std::array<BoundaryCondition, dim>& bcs)
	{
		chemical.initialize(lower,upper,disc, bcs);
		auxiliary.initialize(lower,upper,disc, bcs);
	}


	// access
	template<int dim>
	double FDMChemical<dim>::value(const Point<dim>& p) const
	{
		return chemical.at(p);
	}

	template<int dim>
	double& FDMChemical<dim>::at(const Point<dim> &p)
	{
		return chemical.at(p);
	}

	template<int dim>
	double FDMChemical<dim>::getDiffusionConstant() const
	{ 
		return diffusion_constant; 
	}

	template<int dim>
	double FDMChemical<dim>::getDecayConstant() const
	{ 
		return decay_constant; 
	}

	// mutators:
	template<int dim>
	void FDMChemical<dim>::setDiffusionConstant(double diff)
	{
		diffusion_constant = diff;
	}

	template<int dim>
	void FDMChemical<dim>::setDecayConstant(double decay)
	{
		decay_constant = decay;
	}


	template<int dim>
	void FDMChemical<dim>::projectFunction(const DiscreteFunction<dim>& function)
	{
		chemical.projectFunction(function);
	}


	// evolve:
	template<int dim>
	void FDMChemical<dim>::secreteTo(const Point<dim> p, double amount)
	{
		chemical.at(p) = chemical.at(p) + amount;
	}

	template<int dim>
	void FDMChemical<dim>::update(double time_step, const AdvectionField<dim>* const advection)
	{
		updateDiffusion();

		// std::cout << "Decay constant: " << decay_constant << std::endl;

		if(advection != NULL)
			updateAdvection(advection);

		if(dim == 2)
		{
			// std::cout << "CHEMICAL: " << std::endl;
			// chemical.print(std::cout);
			// std::cout << "\n\nauxiliary:\n\n" << std::endl;
			// auxiliary.print(std::cout);

			// std::cout << "updating chemical:" << std::endl;
			for(unsigned int i = 0; i < chemical.getDiscretization(0); i++)
				for(unsigned int j = 0; j < chemical.getDiscretization(1); j++)
					chemical.at(i,j) = (1.0-time_step*decay_constant)*chemical.at(i,j)
						+ time_step*auxiliary.at(i,j);
		}
		else if(dim == 3)
		{
			for(unsigned int i = 0; i < chemical.getDiscretization(0); i++)
				for(unsigned int j = 0; j < chemical.getDiscretization(1); j++)
					for(unsigned int k = 0; k < chemical.getDiscretization(2); k++)
						chemical.at(i,j,k) = (1.0-time_step*decay_constant)*chemical.at(i,j,k) 
							+ time_step*auxiliary.at(i,j,k);
		}
		// chemical = chemical + time_step*auxiliary // can we overload these operations???
	}

	template<int dim>
	void FDMChemical<dim>::print(std::ostream& out) const
	{
		chemical.print(out);
	}

	template<int dim>
	void FDMChemical<dim>::printAux(std::ostream& out) const
	{
		auxiliary.print(out);
	}

	// private:
	template<int dim>
	void FDMChemical<dim>::updateDiffusion()
	{
		// clear field:
		// auxiliary.reinitialize();

		if(dim == 2)
		{
			const double inv_dx_sqr = chemical.getInverseCellWidthSquared(0);
			const double inv_dy_sqr = chemical.getInverseCellWidthSquared(1);

			// std::cout <<" diffusion_constant: " << diffusion_constant << std::endl;
			// std::cout << " inv_dx_sqr: " << inv_dx_sqr << std::endl;
			// std::cout << " inv_dy_sqr: " << inv_dy_sqr << std::endl << std::endl;

			for(int i = 0; i < chemical.getDiscretization(0); i++)
				for(int j = 0; j < chemical.getDiscretization(1); j++)
					auxiliary.at(i,j) = diffusion_constant*(
						inv_dx_sqr*(chemical.at_bc(i+1,j) - 2*chemical.at(i,j)
							+chemical.at_bc(i-1,j))
						+inv_dy_sqr*(chemical.at_bc(i,j+1) - 2*chemical.at(i,j)
							+chemical.at_bc(i,j-1))
						);
		}
		else if(dim == 3)
		{
			const double inv_dx_sqr = chemical.getInverseCellWidthSquared(0);
			const double inv_dy_sqr = chemical.getInverseCellWidthSquared(1);
			const double inv_dz_sqr = chemical.getInverseCellWidthSquared(2);

			for(int i = 0; i < chemical.getDiscretization(0); i++)
				for(int j = 0; j < chemical.getDiscretization(1); j++)
					for(int k = 0; k < chemical.getDiscretization(2); k++)
						auxiliary.at(i,j,k) = diffusion_constant*(
							inv_dx_sqr*( chemical.at_bc(i+1,j,k) - 2*chemical.at(i,j,k)
								+chemical.at_bc(i-1,j,k) )
							+inv_dy_sqr*( chemical.at_bc(i,j+1,k) - 2*chemical.at(i,j,k)
								+chemical.at_bc(i,j-1,k) )
							+inv_dz_sqr*( chemical.at_bc(i,j,k+1) - 2*chemical.at(i,j,k)
								+chemical.at_bc(i,j,k-1) )
							);
		}
		else
		{
			throw std::invalid_argument("Chemical diffusion only implementted for 2 or 3 dimensions.");
		}

	}

	template<int dim>
	void FDMChemical<dim>::updateAdvection( const AdvectionField<dim>* const advection)
	{
		// upwind scheme
		// performed after diffusion, add to auxiliary field
		if(dim == 2)
		{
			const double inv_dx = auxiliary.getInverseCellWidth(0);
			const double inv_dy = auxiliary.getInverseCellWidth(1);

			for(int i = 0; i < auxiliary.getDiscretization(0); i++)
				for(int j = 0; j < auxiliary.getDiscretization(1); j++)
				{
					const Tensor<1,dim> velocity = advection->value( auxiliary.getCellCenterPoint(i,j) );

					double v_grad_x = 0;
					double v_grad_y = 0;

					// x
					if( velocity[0] > 0 )
					{
						v_grad_x = velocity[0]*inv_dx*( 
							chemical.at_bc(i,j) - chemical.at_bc(i-1,j)
							);
					}
					else
					{
						v_grad_x = velocity[0]*inv_dx*( 
							chemical.at_bc(i+1,j) - chemical.at_bc(i,j)
							);
					}
					// y
					if( velocity[1] > 0 )
					{
						v_grad_y = velocity[1]*inv_dy*( 
							chemical.at_bc(i,j) - chemical.at_bc(i,j-1)
							);
					}
					else
					{
						v_grad_y = velocity[1]*inv_dy*( 
							chemical.at_bc(i,j+1) - chemical.at_bc(i,j)
							);
					}
					auxiliary.at(i,j) += - v_grad_x - v_grad_y; 
				}
		}
		else if(dim == 3)
		{
			const double inv_dx = auxiliary.getInverseCellWidth(0);
			const double inv_dy = auxiliary.getInverseCellWidth(1);
			const double inv_dz = auxiliary.getInverseCellWidth(2);

			for(int i = 0; i < auxiliary.getDiscretization(0); i++)
				for(int j = 0; j < auxiliary.getDiscretization(1); j++)
					for(int k = 0; k < auxiliary.getDiscretization(2); k++)
					{
						const Tensor<1,dim> velocity = advection->value( auxiliary.getCellCenterPoint(i,j,k) );

						double v_grad_x = 0;
						double v_grad_y = 0;
						double v_grad_z = 0;

						// x
						if( velocity[0] > 0 )
						{
							v_grad_x = velocity[0]*inv_dx*( 
								chemical.at_bc(i,j,k) - chemical.at_bc(i-1,j,k)
								);
						}
						else
						{
							v_grad_x = velocity[0]*inv_dx*( 
								chemical.at_bc(i+1,j,k) - chemical.at_bc(i,j,k)
								);
						}
						// y
						if( velocity[1] > 0 )
						{
							v_grad_y = velocity[1]*inv_dy*( 
								chemical.at_bc(i,j,k) - chemical.at_bc(i,j-1,k)
								);
						}
						else
						{
							v_grad_y = velocity[1]*inv_dy*( 
								chemical.at_bc(i,j+1,k) - chemical.at_bc(i,j,k)
								);
						}
						// z
						if( velocity[2] > 0 )
						{
							v_grad_z = velocity[2]*inv_dz*( 
								chemical.at_bc(i,j,k) - chemical.at_bc(i,j,k-1)
								);
						}
						else
						{
							v_grad_z = velocity[2]*inv_dz*( 
								chemical.at_bc(i,j,k+1) - chemical.at_bc(i,j,k)
								);
						}
						auxiliary.at(i,j,k) += - v_grad_x - v_grad_y - v_grad_z; 
					}
		}
		else
		{
			throw std::invalid_argument("Chemical advection only implementted for 2 or 3 dimensions.");
		}

	}

	template<int dim>
	void 
	FDMChemical<dim>::printInfo(std::ostream& out) const
	{
		out << "\n\n-----------------------------------------------------" << std::endl
			<< "\t\t FINITE DIFFERENCE CHEMICAL: "
			<< "\n-----------------------------------------------------" << std::endl
			<< "Diffusion constant: " <<  diffusion_constant << std::endl
			<< "Decay constant: " << decay_constant << std::endl;

		chemical.printInfo(out);
	}
}


#endif


