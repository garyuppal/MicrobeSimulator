#ifndef MICROBESIMULATOR_STOKES_SOLVER_H
#define MICROBESIMULATOR_STOKES_SOLVER_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
// #include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>

#include <iostream>
#include <fstream>
#include <memory>

#include "./velocity_interface.h"
#include "../utility/fe_tools.h"

namespace MicrobeSimulator{ namespace Velocity{
	using namespace dealii;

/** Define seperate preconditioner types for 2 vs 3d
* construct using smart pointer
*/

template <int dim>
struct InnerPreconditioner;

template <>
struct InnerPreconditioner<2>
{
	typedef SparseDirectUMFPACK type;
};

template <>
struct InnerPreconditioner<3>
{
	typedef SparseILU<double> type;
};


//BOUNDARY CONDITIONS
//--------------------------------------------------------------------------------------------

  // INLET BOUNDARY:
  template <int dim>
  class InletBoundary : public Function<dim>
  {
  public:
    InletBoundary (double inlet_vel) 
    	: 
    	Function<dim>(dim+1),
    	inlet_velocity(inlet_vel)
	{}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
private:
	double inlet_velocity;
  };


  template <int dim>
  double
  InletBoundary<dim>::value (const Point<dim>& /* p */,
                              const unsigned int component) const
  {
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));
    if (component == 0)
      return inlet_velocity; // along x direction
    return 0;
  }


  template <int dim>
  void
  InletBoundary<dim>::vector_value (const Point<dim> &p,
                                     Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = InletBoundary<dim>::value (p, c);
  }



  // NO SLIP:
  template <int dim>
  class NoSlipBoundary : public Function<dim>
  {
  public:
    NoSlipBoundary () : Function<dim>(dim+1) {}
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
  };


  template <int dim>
  double
  NoSlipBoundary<dim>::value (const Point<dim>& /* p */,
                              const unsigned int  /* component */ ) const
  {
    // Assert (component < this->n_components,
    //         ExcIndexRange (component, 0, this->n_components));
    // if (component == 0)
    //   return 1; // (p[0] < 0 ? -1 : (p[0] > 0 ? 1 : 0));
    return 0;
  }

  template <int dim>
  void
  NoSlipBoundary<dim>::vector_value (const Point<dim> &p,
                                     Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = NoSlipBoundary<dim>::value (p, c);
  }

//RIGHT HAND SIDE
//--------------------------------------------------------------------------------------------

  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>(dim+1) {}
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
  };

  template <int dim>
  double
  RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                             const unsigned int /*component*/) const
  {
    return 0;
  }

  template <int dim>
  void
  RightHandSide<dim>::vector_value (const Point<dim> &p,
                                    Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = RightHandSide<dim>::value (p, c);
  }


// INVERSE MATRIX
//--------------------------------------------------------------------------------------------
  template <class MatrixType, class PreconditionerType>
  class InverseMatrix : public Subscriptor
  {
  public:
    InverseMatrix (const MatrixType         &m,
                   const PreconditionerType &preconditioner);
    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;
  private:
    const SmartPointer<const MatrixType> matrix;
    const SmartPointer<const PreconditionerType> preconditioner;
  };


  template <class MatrixType, class PreconditionerType>
  InverseMatrix<MatrixType,PreconditionerType>::InverseMatrix
  (const MatrixType         &m,
   const PreconditionerType &preconditioner)
    :
    matrix (&m),
    preconditioner (&preconditioner)
  {}


  template <class MatrixType, class PreconditionerType>
  void InverseMatrix<MatrixType,PreconditionerType>::vmult
  (Vector<double>       &dst,
   const Vector<double> &src) const
  {
    SolverControl solver_control (src.size(), 1e-6*src.l2_norm()); // *** large tolerance for speed
    SolverCG<>    cg (solver_control);
    dst = 0;
    cg.solve (*matrix, dst, src, *preconditioner);
  }



// SCHUR COMPLEMENT
//--------------------------------------------------------------------------------------------
  template <class PreconditionerType>
  class SchurComplement : public Subscriptor
  {
  public:
    SchurComplement (const BlockSparseMatrix<double> &stokes_system_matrix,
                     const InverseMatrix<SparseMatrix<double>, PreconditionerType> &A_inverse);
    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;
  private:
    const SmartPointer<const BlockSparseMatrix<double> > stokes_system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>, PreconditionerType> > A_inverse;
    mutable Vector<double> tmp1, tmp2;
  };


  template <class PreconditionerType>
  SchurComplement<PreconditionerType>::SchurComplement
  (const BlockSparseMatrix<double>                              &stokes_system_matrix,
   const InverseMatrix<SparseMatrix<double>,PreconditionerType> &A_inverse)
    :
    stokes_system_matrix (&stokes_system_matrix),
    A_inverse (&A_inverse),
    tmp1 (stokes_system_matrix.block(0,0).m()),
    tmp2 (stokes_system_matrix.block(0,0).m())
  {}


  template <class PreconditionerType>
  void SchurComplement<PreconditionerType>::vmult (Vector<double>       &dst,
                                                   const Vector<double> &src) const
  {
    stokes_system_matrix->block(0,1).vmult (tmp1, src);
    A_inverse->vmult (tmp2, tmp1);
    stokes_system_matrix->block(1,0).vmult (dst, tmp2);
  }

//------------------------------------------------------------------------------------------------

/** STOKES SOLVER:
* This class solves the stokes equation for the simulation
* and provides access to the flow field
* This class inherits from the VelocityInterface class. This requires we provide
* implementation for the value() method which provides access to the solution.
* For of FE_Chemical however, we'd also like to have access to the block vector
* solution to construct a more accurate system matrix without needing to interpolate the 
* solution. That is, we use the same grid and take the Finite Element solution given here
*/


template<int dim>
class StokesSolver : public VelocityInterface<dim>{
public:
	StokesSolver(Triangulation<dim>& tria, unsigned int fe_degree = 1);

	// constructor to use for advection_handler:
	// since base class does not have other setup methods
	StokesSolver(Triangulation<dim>& tria, 
				const std::vector<unsigned int>& noslip,
				unsigned int n_ref_cycles, 
				double inlet_vel, 
				const std::string& output_directory);

	virtual Tensor<1, dim> value(const Point<dim>& location) const;

	// provide boundary ids for no slip condition:
	void setNoSlipIDs(const std::vector<unsigned int>& ids); 

	void setNumberRefinementSteps(unsigned int n_ref);

	// void attach_point_cell_map(const FETools::PointCellMap<dim>& pcm);

	void solve(const double inlet_velocity,
				const std::string& output_directory = ".",
				const bool output_solution = true);

	// other ACCESSORS:
	double get_maximum_velocity(double max_coordinate) const override; 
	std::vector<std::vector<Tensor<1, dim> > >
		get_fe_velocity_values(unsigned int fe_degree, 
			const FiniteElement<dim>& fe,
			const DoFHandler<dim>& dof) const override;
	
	// not included in base:
	const FESystem<dim>&		get_fe() const;
	unsigned int 				get_fe_degree() const;
	const BlockVector<double>&	get_solution() const;
	const DoFHandler<dim>&		get_dof_handler() const;
	bool						isActive() const;

	void printInfo(std::ostream& out) const override;
private:
	Triangulation<dim>* 		triangulation; // can be modified
	const unsigned int 			stokes_degree;
	unsigned int 				number_refinement_cycles;

	// const FETools::PointCellMap<dim>*	point_cell_map;

	std::vector<unsigned int>	no_slip_boundaries;
	FESystem<dim>				stokes_fe;
	DoFHandler<dim>				stokes_dof_handler;

	AffineConstraints<double> 	stokes_constraints;
	BlockSparsityPattern 		stokes_sparsity_pattern;
	BlockSparseMatrix<double>	stokes_system_matrix;

	BlockSparsityPattern 		preconditioner_sparsity_pattern;
	BlockSparseMatrix<double>	preconditioner_matrix;

	BlockVector<double> 		stokes_solution;
	BlockVector<double> 		stokes_rhs;

	std::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
	bool						active; 			

	// helper methods:
	void 	setup_dofs(const double inlet_velocity);
	void 	assemble_stokes_system();
	void 	solve_stokes();
	void	output_results(const std::string& output_directory,
							const unsigned int refinement_cycle) const;
	void 	refine_mesh();
};

// IMPLEMENTATION:
// ------------------------------------------------------------------------------------
template<int dim>
StokesSolver<dim>::StokesSolver(Triangulation<dim>& tria, unsigned int fe_degree)
	:
	triangulation(&tria),
	stokes_degree(fe_degree),
	number_refinement_cycles(0),
	stokes_fe(FE_Q<dim>(stokes_degree+1), dim, FE_Q<dim>(stokes_degree), 1), // taylor-hood elements
	stokes_dof_handler(*triangulation),
	active(false)
{}

template<int dim>
StokesSolver<dim>::StokesSolver(Triangulation<dim>& tria, 
								const std::vector<unsigned int>& noslip,
								unsigned int n_ref_cycles, 
								double inlet_vel, 
								const std::string& output_directory)
	:
	triangulation(&tria),
	stokes_degree(1),
	number_refinement_cycles(n_ref_cycles),
	no_slip_boundaries(noslip),
	stokes_fe(FE_Q<dim>(stokes_degree+1), dim, FE_Q<dim>(stokes_degree), 1), // taylor-hood elements
	stokes_dof_handler(*triangulation),
	active(false)
{
	// std::cout << "checking no slip boundaries..." << std::endl;
	// for(unsigned int i = 0; i < no_slip_boundaries.size(); ++i)
	// 	std::cout << "no slip id: " << no_slip_boundaries[i] << std::endl;
	// std::cout << ".........................................." << std::endl;
	
	solve(inlet_vel, output_directory, true);
}
								

/** Interfacing value function
* This can be sped up if provided with a cell mapping, from which we can cache solutions for a 
* geometry.... can also make use of a discrete field object ...
*/
template<int dim>
Tensor<1, dim>
StokesSolver<dim>::value(const Point<dim>& location) const
{
	// double warning_fix_stokes_quick_access;
	
	if(active == false)
		return Point<dim>(); //Tensor<1,dim>(); // zero default constructor

	Vector<double> velocity;
	velocity.reinit(dim + 1); // +1 for pressure

	VectorTools::point_value(stokes_dof_handler,
							stokes_solution,
							location,
							velocity);

	// FETools::vector_point_value(point_cell_map,
	// 							stokes_dof_handler,
	// 							stokes_solution,
	// 							location,
	// 							velocity);

	// Tensor<1,dim> return_value;
	Tensor<1, dim> return_value;
	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
		return_value[dim_itr] = velocity[dim_itr];

	return return_value;
}

/** OVERRIDE method for extracting velocity values at quadrature points
* for chemical matrix construction
*/
template<int dim>
std::vector<std::vector<dealii::Tensor<1, dim> > >
StokesSolver<dim>::get_fe_velocity_values(unsigned int fe_degree, 
									const FiniteElement<dim>& /* fe */,
									const DoFHandler<dim>& dof) const
{

	const QGauss<dim> quadrature_formula(fe_degree + 2);
	FEValues<dim> stokes_fe_values(stokes_fe, quadrature_formula, update_values);

	const unsigned int n_q_points = quadrature_formula.size();

	std::vector<std::vector<dealii::Tensor<1, dim> > > velocity_values;
	velocity_values.reserve(dof.n_dofs());
	std::vector<Tensor<1, dim> > 		cell_velocity_values(n_q_points);

	const FEValuesExtractors::Vector velocities(0);

	auto cell 			= dof.begin_active();
	const auto endc		= dof.end();
	auto stokes_cell	= stokes_dof_handler.begin_active();
	unsigned int cell_id = 0;

	for(; cell != endc; ++cell, ++stokes_cell, ++cell_id)
	{
		stokes_fe_values.reinit(stokes_cell);
		stokes_fe_values[velocities].get_function_values(stokes_solution, 
														cell_velocity_values);
		velocity_values.emplace_back(cell_velocity_values);
	}

	return velocity_values;
}

// template<int dim>
// void 
// StokesSolver<dim>::attach_point_cell_map(const FETools::PointCellMap<dim>& pcm)
// {
// 	point_cell_map = &pcm;
// }

template<int dim>
void 
StokesSolver<dim>::setNoSlipIDs(const std::vector<unsigned int>& ids)
{
	no_slip_boundaries = ids;
}


template<int dim>
void 
StokesSolver<dim>::setNumberRefinementSteps(unsigned int n_ref)
{
	number_refinement_cycles = n_ref;
}

/** Setup system degrees of freedom:  
* make sure no slip boundaries are set before calling
*/

template<int dim>
void 
StokesSolver<dim>::setup_dofs(const double inlet_velocity)
{
	A_preconditioner.reset ();
	stokes_system_matrix.clear ();
	preconditioner_matrix.clear ();
	stokes_dof_handler.distribute_dofs (stokes_fe);

	/// renumber DOFS to help ILU ***
	DoFRenumbering::Cuthill_McKee (stokes_dof_handler); 

	std::vector<unsigned int> block_component (dim+1,0);
	block_component[dim] = 1;
	DoFRenumbering::component_wise (stokes_dof_handler, block_component);

	// DIRECHLET BOUNDARY CONDITIONS:
	{
		const unsigned int inlet_boundary_id = (inlet_velocity >=0)? 0: 1; // 1 is open, flip if negative velocity

		stokes_constraints.clear ();
		FEValuesExtractors::Vector velocities(0);
		DoFTools::make_hanging_node_constraints (stokes_dof_handler,
		                                     stokes_constraints);

		/// constant velocity at inlet:
		VectorTools::interpolate_boundary_values (stokes_dof_handler,
		                                      inlet_boundary_id,  
		                                      InletBoundary<dim>(inlet_velocity),
		                                      stokes_constraints,
		                                      stokes_fe.component_mask(velocities)); 

		// no slip boundary conditions on walls:
		for(unsigned int i = 0; i < no_slip_boundaries.size(); ++i)
			VectorTools::interpolate_boundary_values (stokes_dof_handler,
			                                        no_slip_boundaries[i],  
			                                        NoSlipBoundary<dim>(),
			                                        stokes_constraints,
			                                        stokes_fe.component_mask(velocities)); 

	}
	stokes_constraints.close ();

	std::vector<types::global_dof_index> dofs_per_block (2);
	DoFTools::count_dofs_per_block (stokes_dof_handler, 
	                              dofs_per_block, 
	                              block_component);
	const unsigned int n_u = dofs_per_block[0],
	                 n_p = dofs_per_block[1];
	std::cout << "   Number of active cells: "
	        << triangulation->n_active_cells()
	        << std::endl
	        << "   Number of degrees of freedom: "
	        << stokes_dof_handler.n_dofs()
	        << " (" << n_u << '+' << n_p << ')'
	        << std::endl;

	/// scoping helps to release memory after operation complete
	{
		BlockDynamicSparsityPattern dsp (2,2);
		dsp.block(0,0).reinit (n_u, n_u);
		dsp.block(1,0).reinit (n_p, n_u);
		dsp.block(0,1).reinit (n_u, n_p);
		dsp.block(1,1).reinit (n_p, n_p);
		dsp.collect_sizes();
		Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
		for (unsigned int c=0; c<dim+1; ++c)
		  for (unsigned int d=0; d<dim+1; ++d)
		    if (! ((c==dim) && (d==dim)))
		      coupling[c][d] = DoFTools::always;
		    else
		      coupling[c][d] = DoFTools::none;
		DoFTools::make_sparsity_pattern (stokes_dof_handler, 
		                                coupling, dsp, 
		                                stokes_constraints, 
		                                false);
		stokes_sparsity_pattern.copy_from (dsp);
	}

	{
		BlockDynamicSparsityPattern preconditioner_dsp (2,2);
		preconditioner_dsp.block(0,0).reinit (n_u, n_u);
		preconditioner_dsp.block(1,0).reinit (n_p, n_u);
		preconditioner_dsp.block(0,1).reinit (n_u, n_p);
		preconditioner_dsp.block(1,1).reinit (n_p, n_p);
		preconditioner_dsp.collect_sizes();
		Table<2,DoFTools::Coupling> preconditioner_coupling (dim+1, dim+1);
		for (unsigned int c=0; c<dim+1; ++c)
		  for (unsigned int d=0; d<dim+1; ++d)
		    if (((c==dim) && (d==dim)))
		      preconditioner_coupling[c][d] = DoFTools::always;
		    else
		      preconditioner_coupling[c][d] = DoFTools::none;
		DoFTools::make_sparsity_pattern (stokes_dof_handler, preconditioner_coupling,
		                                 preconditioner_dsp, stokes_constraints, false);
		preconditioner_sparsity_pattern.copy_from (preconditioner_dsp);
	}

	// initialize objects:
	stokes_system_matrix.reinit (stokes_sparsity_pattern);
	preconditioner_matrix.reinit (preconditioner_sparsity_pattern);

	stokes_solution.reinit (2);
	stokes_solution.block(0).reinit (n_u);
	stokes_solution.block(1).reinit (n_p);
	stokes_solution.collect_sizes ();

	stokes_rhs.reinit (2);
	stokes_rhs.block(0).reinit (n_u);
	stokes_rhs.block(1).reinit (n_p);
	stokes_rhs.collect_sizes ();
} /// setup_dofs() 


template<int dim>
void 
StokesSolver<dim>::assemble_stokes_system()
{
	stokes_system_matrix=0;
	stokes_rhs=0;
	preconditioner_matrix = 0;

	QGauss<dim>   quadrature_formula(stokes_degree+2);
	FEValues<dim> stokes_fe_values (stokes_fe, quadrature_formula,
									update_values    |
									update_quadrature_points  |
									update_JxW_values |
									update_gradients);

	const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
	const unsigned int   n_q_points      = quadrature_formula.size();

	FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
	FullMatrix<double>   local_preconditioner_matrix (dofs_per_cell, dofs_per_cell);

	Vector<double>       local_rhs (dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

	const RightHandSide<dim>          right_hand_side; // *** defined in program above
	std::vector<Vector<double> >      rhs_values (n_q_points,
	                                            Vector<double>(dim+1));

	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);

	std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
	std::vector<double>                  div_phi_u   (dofs_per_cell);
	std::vector<double>                  phi_p       (dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator
	cell = stokes_dof_handler.begin_active(),
	endc = stokes_dof_handler.end();

	for (; cell!=endc; ++cell)
	{
	  stokes_fe_values.reinit (cell);
	  local_matrix = 0;
	  local_preconditioner_matrix = 0;
	  local_rhs = 0;
	  right_hand_side.vector_value_list(stokes_fe_values.get_quadrature_points(),
	                                    rhs_values);
	  for (unsigned int q=0; q<n_q_points; ++q)
	    {
	      for (unsigned int k=0; k<dofs_per_cell; ++k)
	        {
	          symgrad_phi_u[k] = stokes_fe_values[velocities].symmetric_gradient (k, q);
	          div_phi_u[k]     = stokes_fe_values[velocities].divergence (k, q);
	          phi_p[k]         = stokes_fe_values[pressure].value (k, q);
	        }
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
	        { // *** BILINEAR FORMS ***
	          for (unsigned int j=0; j<=i; ++j) // *** loop upto i only, because symmetric ***
	            {
	              local_matrix(i,j) += (2 * (symgrad_phi_u[i] * symgrad_phi_u[j])
	                                    - div_phi_u[i] * phi_p[j]
	                                    - phi_p[i] * div_phi_u[j])
	                                   * stokes_fe_values.JxW(q);
	              local_preconditioner_matrix(i,j) += (phi_p[i] * phi_p[j])
	                                                  * stokes_fe_values.JxW(q);
	            }
	          const unsigned int component_i =
	            stokes_fe.system_to_component_index(i).first;
	          local_rhs(i) += stokes_fe_values.shape_value(i,q) *
	                          rhs_values[q](component_i) *
	                          stokes_fe_values.JxW(q);
	        }
	    }
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=i+1; j<dofs_per_cell; ++j)
	      {
	        local_matrix(i,j) = local_matrix(j,i);
	        local_preconditioner_matrix(i,j) = local_preconditioner_matrix(j,i); 
	        // flip indices to fill in rest ***
	      }
	  cell->get_dof_indices (local_dof_indices);
	  stokes_constraints.distribute_local_to_global (local_matrix, local_rhs,
	                                          local_dof_indices,
	                                          stokes_system_matrix, stokes_rhs);
	  stokes_constraints.distribute_local_to_global (local_preconditioner_matrix,
	                                          local_dof_indices,
	                                          preconditioner_matrix);
	}
	// *** make preconditioner  
	std::cout << "   Computing preconditioner..." << std::endl << std::flush;

	A_preconditioner
	= std::make_shared<typename InnerPreconditioner<dim>::type>();
	A_preconditioner->initialize (stokes_system_matrix.block(0,0),
	                            typename InnerPreconditioner<dim>::type::AdditionalData());
}


template<int dim>
void 
StokesSolver<dim>::solve_stokes()
{
	const InverseMatrix<SparseMatrix<double>,
	    typename InnerPreconditioner<dim>::type>
	    A_inverse (stokes_system_matrix.block(0,0), *A_preconditioner);

	Vector<double> tmp (stokes_solution.block(0).size());

	{
		Vector<double> schur_rhs (stokes_solution.block(1).size());
		A_inverse.vmult (tmp, stokes_rhs.block(0));
		stokes_system_matrix.block(1,0).vmult (schur_rhs, tmp);
		schur_rhs -= stokes_rhs.block(1);
		SchurComplement<typename InnerPreconditioner<dim>::type>
		schur_complement (stokes_system_matrix, A_inverse);
		SolverControl solver_control (stokes_solution.block(1).size(),
		                              1e-6*schur_rhs.l2_norm());
		SolverCG<>    cg (solver_control);
		SparseILU<double> preconditioner;
		preconditioner.initialize (preconditioner_matrix.block(1,1),
		                           SparseILU<double>::AdditionalData());
		InverseMatrix<SparseMatrix<double>,SparseILU<double> >
		m_inverse (preconditioner_matrix.block(1,1), preconditioner);
		cg.solve (schur_complement, stokes_solution.block(1), schur_rhs,
		          m_inverse);
		stokes_constraints.distribute (stokes_solution);
		std::cout << "  "
		          << solver_control.last_step()
		          << " outer CG Schur complement iterations for pressure"
		          << std::endl;
	}
	{
		stokes_system_matrix.block(0,1).vmult (tmp, stokes_solution.block(1));
		tmp *= -1;
		tmp += stokes_rhs.block(0);
		A_inverse.vmult (stokes_solution.block(0), tmp);
		stokes_constraints.distribute (stokes_solution);
	}
}


template<int dim>
void 
StokesSolver<dim>::output_results(const std::string& output_directory,
									const unsigned int refinement_cycle) const
{
	std::vector<std::string> stokes_solution_names (dim, "velocity");
	stokes_solution_names.emplace_back("pressure");

	std::vector<DataComponentInterpretation::DataComponentInterpretation>
		data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
	data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);

	DataOut<dim> data_out;
	data_out.attach_dof_handler (stokes_dof_handler);
	data_out.add_data_vector (stokes_solution, 
							stokes_solution_names,
	                        DataOut<dim>::type_dof_data,
	                        data_component_interpretation);
	data_out.build_patches ();

	std::string file_name = output_directory
						+ "/stokes_solution-"
	                    + Utilities::int_to_string(refinement_cycle, 2)
	                    + ".vtk";
    std::cout << "file: " << file_name << std::endl;
	std::ofstream output (file_name);
	data_out.write_vtk (output);
}


/** use KellYErrorEstimator, only for pressure to flag refinement
* refine regions with large pressure gradients
*/

template<int dim>
void 
StokesSolver<dim>::refine_mesh()
{
	Vector<float> estimated_error_per_cell (triangulation->n_active_cells());
	FEValuesExtractors::Scalar pressure(dim);
	KellyErrorEstimator<dim>::estimate (stokes_dof_handler,
	                                  QGauss<dim-1>(stokes_degree+1),
	                                  typename FunctionMap<dim>::type(),
	                                  stokes_solution,
	                                  estimated_error_per_cell,
	                                  stokes_fe.component_mask(pressure));
	GridRefinement::refine_and_coarsen_fixed_number (*triangulation,
	                                               estimated_error_per_cell,
	                                               0.3, 0.0); /// where do these numbers come from ???
	triangulation->execute_coarsening_and_refinement ();
}


template<int dim>
void 
StokesSolver<dim>::solve(const double inlet_velocity,
						const std::string& output_directory,
						const bool output_solution)
{
	std::cout << "\n\nSOLVING STOKES..." << std::endl << std::endl;

	for (unsigned int refinement_cycle = 0; 
		refinement_cycle < number_refinement_cycles;
		++refinement_cycle)
    {
      std::cout << "Refinement cycle: " << refinement_cycle << std::endl;
      if (refinement_cycle > 0)
        refine_mesh ();
      setup_dofs (inlet_velocity);
      std::cout << "\t Assembling..." << std::endl << std::flush;
      assemble_stokes_system ();
      std::cout << "\t Solving..." << std::flush;
      solve_stokes ();

      if(output_solution == true)
	      output_results (output_directory, refinement_cycle);
      std::cout << std::endl << std::endl;
    }

    active = true; // to inform others (accessing objects) a solution has been obtained
}


/** ACCESSORS FOR STOKES SOLUTION
*/
template <int dim>
double 
StokesSolver<dim>::get_maximum_velocity(double /* max_coordinate */) const
{
	if(!active)
		return 0;

	const QIterated<dim> quadrature_formula(QTrapez<1>(), stokes_degree + 1);
	const unsigned int   n_q_points = quadrature_formula.size();
	FEValues<dim> fe_values(stokes_fe, quadrature_formula, update_values);
	std::vector<Tensor<1, dim> > velocity_values(n_q_points);

	double max_velocity = 0;

	const FEValuesExtractors::Vector velocities(0);
	for (const auto &cell : stokes_dof_handler.active_cell_iterators())
	{
		fe_values.reinit(cell);
		fe_values[velocities].get_function_values(stokes_solution,
		                                        velocity_values);
		for (unsigned int q = 0; q < n_q_points; ++q)
			max_velocity = std::max(max_velocity, velocity_values[q].norm());
	}
	return max_velocity;
}


template<int dim>
const FESystem<dim>& 		
StokesSolver<dim>::get_fe() const
{
	return stokes_fe;
}


template<int dim>
const BlockVector<double>& 		
StokesSolver<dim>::get_solution() const
{
	return stokes_solution;
}


template<int dim>
const DoFHandler<dim>&
StokesSolver<dim>::get_dof_handler() const
{
	return stokes_dof_handler;
}


template<int dim>
bool
StokesSolver<dim>::isActive() const
{
	return active;
}

template<int dim>
void 
StokesSolver<dim>::printInfo(std::ostream& out) const
{
	out << Utility::short_line << std::endl
		<< "\t General Stokes flow:" << std::endl
		<< Utility::short_line << std::endl
		<< "..." << std::endl
		<< Utility::short_line << std::endl << std::endl;
}


}} // close namespace
#endif