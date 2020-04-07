#ifndef MICROBE_SIMULATOR_REFACTORED_CHEMICAL_FE_BASE_H
#define MICROBE_SIMULATOR_REFACTORED_CHEMICAL_FE_BASE_H


/** necessary header files...
* @todo clean up, remove those not used...
*/

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
// can probably remove cg...
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <deal.II/lac/precondition.h>

// #include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/base/tensor_function.h>
#include <deal.II/grid/grid_in.h>

#include <fstream>
#include <sstream>
#include <limits>
#include <memory>

#include "../advection/velocity_interface.h"
#include "../advection/stokes_solver.h"

// remove these ...
#include "../utility/cell_iterator_map.h"
#include "../utility/fe_tools.h"

namespace MicrobeSimulator{ namespace RefactoredChemicals{
	using namespace dealii;

template<int dim>
class Chemical_FE_Base{  
public:
	// Chemical_FE_Base() 	: triangulation(NULL),
	// chemical_fe_degree(0),
	// chemical_fe( 0 ),
	// chemical_dof_handler(*triangulation),
	// using_velocity(false) {} // remove, to check for now

	Chemical_FE_Base(const Triangulation<dim>& tria,
					const unsigned int degree = 1); // done

	// void setup(const VelocityInterface<dim>& velocity_function); // still add velocity setup...
	// void setup(const Velocity::AdvectionHandler<dim>& velocity_function); // done
	void setup(const Velocity::AdvectionHandler<dim>& velocity_function,
				const std::array<BoundaryCondition, dim>& bcs); // done

	// legacy:
	void attach_point_cell_map(const FETools::PointCellMap<dim>& pcm); // done
	
	void attach_point_cell_map(std::shared_ptr<FETools::PointCellMap<dim> > pcm); // done


	/// ACCESSORS:
	// const FETools::PointCellMap<dim>* getPointCellMap() const; // done
	std::shared_ptr<FETools::PointCellMap<dim> > getPointCellMap() const; // done // ***not acutally const

	const Triangulation<dim>& get_triangulation() const;
	const SparsityPattern& getSparsityPattern() const; // done
	unsigned int get_n_dofs() const; // done
	const SparseMatrix<double>& get_mass_matrix() const; // done
	const SparseMatrix<double>& get_stiffness_matrix() const; // done
	unsigned int get_fe_degree() const; // done
	const FE_Q<dim>& get_fe() const; // done
	const DoFHandler<dim>& get_dof_handler() const; // done
	const ConstraintMatrix& get_constraints() const; // done -- change to affine...
	const Vector<double>& get_mass_check_vector() const; // done

	// const Velocity::AdvectionHandler<dim>* getStokesSolution() const;

	// velocity values: -- add velocity function support ...
	std::vector<Tensor<1, dim> >	
	get_cell_velocity_values(unsigned int cell_id) const; // done

	Tensor<1, dim>		
	get_cell_quad_velocity_value(unsigned int cell_id, unsigned int q_point) const; // done

private:
	// common to all chemicals:

	const Triangulation<dim>* const		triangulation;
	const unsigned int      			chemical_fe_degree;
	FE_Q<dim>               			chemical_fe;
	DoFHandler<dim>         			chemical_dof_handler;

	// const FETools::PointCellMap<dim>*	point_cell_map;
	std::shared_ptr<FETools::PointCellMap<dim> > point_cell_map;

	ConstraintMatrix        			chemical_constraints;
	SparsityPattern         			chemical_sparsity_pattern;

	SparseMatrix<double>    			chemical_mass_matrix;
	SparseMatrix<double>    			chemical_diffusion_matrix;

	Vector<double>						mass_check_vector;

	// velocity: (extracted for chemical matrix assembly)
	std::vector<std::vector<Tensor<1, dim> > >	velocity_values; //(n_q_points);
	bool using_velocity;

	// initialization functions
	// void setup_dofs(); // done
	void setup_dofs(const std::array<BoundaryCondition, dim>& bcs); // done
	void assemble_chemical_matrices(); // done -- makes mass and diffusion matrices...
	// void setup_velocity_values(const VelocityInterface<dim>& velocity_function);
	void setup_velocity_values(const Velocity::AdvectionHandler<dim>& velocity_function); // done
	void create_mass_check_vector(); // done
};


// IMPLEMENTATION:
//-------------------------------------------------------------------------
template<int dim>
Chemical_FE_Base<dim>::Chemical_FE_Base(const Triangulation<dim>& tria,
										const unsigned int degree)
	:
	triangulation(&tria),
	chemical_fe_degree(degree),
	chemical_fe( chemical_fe_degree ),
	chemical_dof_handler(*triangulation),
	using_velocity(false)
{}


// template<int dim>
// void 
// Chemical_FE_Base<dim>::setup(const VelocityInterface<dim>& velocity_function)
// {
// 	std::cout << "...setting up chemical base with velocity function" << std::endl;

// 	setup_dofs();
// 	assemble_chemical_matrices();
// 	create_mass_check_vector();

// 	setup_velocity_values(velocity_function);
// }


// remove this and use one below in general:
// template<int dim>
// void 
// Chemical_FE_Base<dim>::setup(const Velocity::AdvectionHandler<dim>& velocity_function)
// {
// 	std::cout << "...setting up chemical base with stokes velocity" << std::endl;

// 	setup_dofs();
// 	assemble_chemical_matrices();
// 	create_mass_check_vector();

// 	// stokes = &velocity_function;

// 	// if(velocity_function.isActive())
// 		setup_velocity_values(velocity_function);
// }

// to take over setup method above ***
template<int dim>
void 
Chemical_FE_Base<dim>::setup(const Velocity::AdvectionHandler<dim>& velocity_function,
			const std::array<BoundaryCondition, dim>& bcs)
{
	std::cout << "...setting up chemical base with stokes velocity" << std::endl;

	setup_dofs(bcs); // also to replace old setup_dofs()
	assemble_chemical_matrices();
	create_mass_check_vector();

	// stokes = &velocity_function;

	// if(velocity_function.isActive())
	setup_velocity_values(velocity_function);
}

template<int dim>
void 
Chemical_FE_Base<dim>::attach_point_cell_map(const FETools::PointCellMap<dim>& pcm)
{
	point_cell_map = &pcm;
}

template<int dim>
void 
Chemical_FE_Base<dim>::attach_point_cell_map(
	std::shared_ptr<FETools::PointCellMap<dim> > pcm)
{
	point_cell_map = pcm;
}

/** setup vector velocity_values for adding advection to chemicals
*/

// template<int dim>
// void
// Chemical_FE_Base<dim>::setup_velocity_values(const VelocityInterface<dim>& velocity_function)
// {
// 	velocity_values = velocity_function.get_fe_velocity_values(chemical_fe_degree,
// 		chemical_fe, chemical_dof_handler);
// }

template<int dim>
void 
Chemical_FE_Base<dim>::setup_velocity_values(const Velocity::AdvectionHandler<dim>& velocity_function)
{
	std::cout << std::endl << "SETTING UP CHEMICAL VELOCITY VALUES " << std::endl;
	velocity_values = velocity_function.get_fe_velocity_values(chemical_fe_degree,
		chemical_fe, chemical_dof_handler);

	using_velocity = true; 
}


// template<int dim>
// void 
// Chemical_FE_Base<dim>::setup_dofs()
// {
// 	{
// 		chemical_dof_handler.distribute_dofs(chemical_fe);
// 		chemical_constraints.clear();
// 		DoFTools::make_hanging_node_constraints(chemical_dof_handler,
// 		                                	    chemical_constraints);

// 		// **** option to make periodic *****

// 		chemical_constraints.close();
// 	}
// 	{
// 		DynamicSparsityPattern dsp(chemical_dof_handler.n_dofs());
// 		DoFTools::make_sparsity_pattern(chemical_dof_handler,
// 			                            dsp,
// 			                            chemical_constraints,
// 			/*keep constrained_dofs = */ false); 
// 		chemical_sparsity_pattern.copy_from(dsp);
// 	}

// 	chemical_mass_matrix.reinit(chemical_sparsity_pattern);
// 	chemical_diffusion_matrix.reinit(chemical_sparsity_pattern);
// }

template<int dim>
void 
Chemical_FE_Base<dim>::setup_dofs(const std::array<BoundaryCondition, dim>& bcs)
{
	{
		chemical_dof_handler.distribute_dofs(chemical_fe);
		chemical_constraints.clear();
		DoFTools::make_hanging_node_constraints(chemical_dof_handler,
		                                	    chemical_constraints);

		// **** option to make periodic *****:
		for(unsigned int i = 0; i < dim; ++i)
		{
			// std::cout << "CHECK HERE IF BC IS PERIODIC" << std::endl;
			if(bcs[i] == BoundaryCondition::WRAP) 
			{
				std::cout << "\n...USING PERIODIC " << i <<"th BOUNDARY" << std::endl;
				// ADD PERIODICITY:
				std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> >
				periodicity_vector;

				const unsigned int direction = i;
				unsigned int bid_one = 0 + 2*i;
				unsigned int bid_two = 1 + 2*i;

				GridTools::collect_periodic_faces(chemical_dof_handler, bid_one, bid_two, direction,
				                            periodicity_vector); //, offset, rotation_matrix);

				DoFTools::make_periodicity_constraints<DoFHandler<dim> > 
				(periodicity_vector, chemical_constraints); //, fe.component_mask(velocities), first_vector_components);
			} // if periodic 
		} // for each dimension

		chemical_constraints.close();
	}
	{
		DynamicSparsityPattern dsp(chemical_dof_handler.n_dofs());
		DoFTools::make_sparsity_pattern(chemical_dof_handler,
			                            dsp,
			                            chemical_constraints,
			/*keep constrained_dofs = */ false); 
		chemical_sparsity_pattern.copy_from(dsp);
	}

	chemical_mass_matrix.reinit(chemical_sparsity_pattern);
	chemical_diffusion_matrix.reinit(chemical_sparsity_pattern);
}


template<int dim>
void 
Chemical_FE_Base<dim>::assemble_chemical_matrices()
{
	std::cout << "...assembling chemical matrices" << std::endl;
	chemical_mass_matrix = 0;
	chemical_diffusion_matrix = 0;

	QGauss<dim> quadrature_formula(chemical_fe_degree + 2);
	FEValues<dim> fe_values(chemical_fe, quadrature_formula,
	      update_values | update_gradients | update_JxW_values);

	const unsigned int dofs_per_cell = chemical_fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();

	FullMatrix<double> local_mass_matrix(dofs_per_cell, dofs_per_cell);
	FullMatrix<double> local_diffusion_matrix(dofs_per_cell, dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	std::vector<double>           phi_T (dofs_per_cell);
	std::vector<Tensor<1, dim> >  grad_phi_T(dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator
		cell = chemical_dof_handler.begin_active(),
		endc = chemical_dof_handler.end();

	for(; cell != endc; ++cell)
	{
		local_mass_matrix = 0;
		local_diffusion_matrix = 0;

		fe_values.reinit(cell);

		for(unsigned int q = 0; q < n_q_points; ++q)
		{
			for(unsigned int k = 0; k < dofs_per_cell; ++k)
			{
				grad_phi_T[k] = fe_values.shape_grad(k,q);
				phi_T[k] = fe_values.shape_value(k,q);
			}

			for(unsigned int i = 0; i < dofs_per_cell; ++i)
				for (unsigned int j = 0; j < dofs_per_cell; ++j)
				{
				local_mass_matrix(i,j) += (phi_T[i] * phi_T[j]
				*fe_values.JxW(q));

				local_diffusion_matrix(i,j) += (grad_phi_T[i] * grad_phi_T[j]
				*fe_values.JxW(q));
				}
		} // for quadrature points

		cell->get_dof_indices(local_dof_indices);
		chemical_constraints.distribute_local_to_global(local_mass_matrix,
						                                local_dof_indices,
						                                chemical_mass_matrix);
		chemical_constraints.distribute_local_to_global(local_diffusion_matrix,
						                                local_dof_indices,
						                                chemical_diffusion_matrix);
	} // for cells
}


template<int dim>
void
Chemical_FE_Base<dim>::create_mass_check_vector()
{
	mass_check_vector.reinit(chemical_dof_handler.n_dofs());

	const QGauss<dim>  quadrature_formula((chemical_fe.degree+1));
	FEValues<dim> fe_values (chemical_fe, quadrature_formula,
			update_values | update_gradients | update_JxW_values);
	const unsigned int dofs_per_cell = chemical_fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();

	Vector<double>       cell_rhs (dofs_per_cell);
	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

	for (const auto &cell: chemical_dof_handler.active_cell_iterators())
	{
		fe_values.reinit (cell);

		cell_rhs = 0;
		for (unsigned int q_index=0; q_index < n_q_points; ++q_index)
		{
			for (unsigned int i=0; i < dofs_per_cell; ++i)
				cell_rhs(i) += (fe_values.shape_value (i, q_index) *
				1 *
				fe_values.JxW (q_index));
		}
		cell->get_dof_indices (local_dof_indices);
		for (unsigned int i=0; i < dofs_per_cell; ++i)
			mass_check_vector(local_dof_indices[i]) += cell_rhs(i);
	}
}


/* ACCESSORS:
* these will be needed to let chemical classes make use of the base class
*/
// ---------------------------------------------------------------------------------

template<int dim>
std::shared_ptr<FETools::PointCellMap<dim> > //const FETools::PointCellMap<dim>* 
Chemical_FE_Base<dim>::getPointCellMap() const
{
	return point_cell_map;
}

template<int dim>
const Triangulation<dim>& 
Chemical_FE_Base<dim>::get_triangulation() const 
{
	return *triangulation;
}

template<int dim>
const SparsityPattern& 
Chemical_FE_Base<dim>::getSparsityPattern() const
{
	return chemical_sparsity_pattern;
}


template<int dim>
unsigned int 
Chemical_FE_Base<dim>::get_n_dofs() const
{
	return chemical_dof_handler.n_dofs();
}


template<int dim>
const SparseMatrix<double>& 
Chemical_FE_Base<dim>::get_mass_matrix() const
{
	return chemical_mass_matrix;
}


template<int dim>
const SparseMatrix<double>& 
Chemical_FE_Base<dim>::get_stiffness_matrix() const
{
	return chemical_diffusion_matrix;
}


template<int dim>
unsigned int 
Chemical_FE_Base<dim>::get_fe_degree() const
{
	return chemical_fe_degree;
}


template<int dim>
const FE_Q<dim>& 
Chemical_FE_Base<dim>::get_fe() const
{
	return chemical_fe;
}


template<int dim>
const DoFHandler<dim>& 
Chemical_FE_Base<dim>::get_dof_handler() const
{
	return chemical_dof_handler;
}


template<int dim>
const ConstraintMatrix& 
Chemical_FE_Base<dim>::get_constraints() const
{
	return chemical_constraints;
}

// velocity values:
template<int dim>
std::vector<Tensor<1, dim> >	
Chemical_FE_Base<dim>::get_cell_velocity_values(unsigned int cell_id) const
{
	if(using_velocity == false)
		return std::vector<Tensor<1,dim> >{Tensor<1, dim>()};

	return velocity_values[cell_id]; 
}


template<int dim>
Tensor<1, dim>		
Chemical_FE_Base<dim>::get_cell_quad_velocity_value(unsigned int cell_id, 
													unsigned int q_point) const
{
	if(using_velocity == false)
		return Tensor<1,dim>(); // zero by default

	// std::cout << "size of velocity values at cell (n_q_points): " 
	// 	<< velocity_values[0].size()
	// 	<< " | cell-qpoint access attempt at: " << cell_id 
	// 	<< ", " << q_point << std::endl;

	return velocity_values[cell_id][q_point]; //.at(cell_id).at(q_point); //[cell_id][q_point]; 
}


template<int dim>
const Vector<double>& 
Chemical_FE_Base<dim>::get_mass_check_vector() const
{
	return mass_check_vector;
}


}} // close namespace
#endif