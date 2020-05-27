#ifndef MICROBE_SIMULATOR_FE_CHEMICAL_H
#define MICROBE_SIMULATOR_FE_CHEMICAL_H

#include <deal.II/base/timer.h>

// direct solver:
#include <deal.II/lac/sparse_direct.h>

#include "./chemical_fe_base.h"
#include "./chemical_interface.h"
#include "./control_functions.h"
#include "../utility/fe_tools.h"

namespace MicrobeSimulator{ namespace Chemicals{
	using namespace dealii;

template<int dim>
class FE_Chemical : public ChemicalInterface<dim>{
public:	
	FE_Chemical(); // done ...

	FE_Chemical(std::shared_ptr<Chemical_FE_Base<dim> > chem_base, //const Chemical_FE_Base<dim>& chem_base, 
		double diffusion, double decay, double b, double dt);
	
	// setup:
	void reinit(const Chemical_FE_Base<dim>& chem_base,
		double diffusion, double decay, double b, double dt); // done

	// interface & update:
	double value(const Point<dim>& p) const override;

	std::vector<double> 
		value_list(const std::vector<Point<dim> >& points) const override;
		
	void project_function(const Function<dim>& initial_condition) override;

	// for evolving without sources (initial condition)
	void update() override;
	void update(const std::vector<Point<dim> >& locations, 
						const std::vector<double>& amounts) override;

	void update(const std::vector<Point<dim> >& locations, 
							const std::vector<double>& amounts,
							const Function<dim>& control_function) override;

	// accessors:
	double 	getDiffusionConstant() const; // done
	double 	getDecayConstant() const; // done
	double 	getViscosityBeta() const; // done
	double 	getTimeStep() const override;  // done

	// methods:
	void 	project_initial_condition(const Function<dim>& initial_condition); // done

	double 	getMass() const; // done 
	double 	getMin() const; // done // for debugging -- using min/max of vector for now...
	double 	getMax() const; // done // for debugging -- using min/max of vector for now...

	void 	output_solution(const std::string& output_directory,
						 const unsigned int chem_id,
						 const unsigned int save_step) const; // done

	void 	output_solution(const std::string& output_directory,
						 const unsigned int chem_id,
						 const unsigned int save_step,
						 const unsigned int cycle) const; 

	void 	print(std::ostream& out, unsigned int chem_id) const override;
	void 	printInfo(std::ostream& out) const override;

private:
	// const Chemical_FE_Base<dim>*	chemical_base;
	std::shared_ptr<Chemical_FE_Base<dim> >					chemical_base; 

	double 							diffusion_constant;
	double 							decay_constant;
	double 							viscosity_beta;
	double 							time_step;

	// chemical objects:
	SparseMatrix<double>			system_matrix;
	Vector<double>					solution;
	Vector<double>					old_solution;
// could still use...
		// Vector<double>					old_old_solution; // for second order time schemes

	Vector<double> 					right_hand_side; 
	Vector<double>					source; 
	Vector<double>					control_function_fem_vector;
	
	SparseDirectUMFPACK			 	A_direct;

	// methods:
	// setup:
	void 	reinit(); // done
	void	setup_solver(); // assemble system matrix and A_direct
    double 	compute_viscosity(const std::vector<Tensor<1,dim> >& velocity_values,
                const double cell_diameter); // done -- run tests to find beta ...

    // update and solve:
	void	update_rhs(); // done
	void	update_rhs(const std::vector<Point<dim> >& locations, 
					const std::vector<double>& amounts); // done
	void	update_rhs(const std::vector<Point<dim> >& locations, 
					const std::vector<double>& amounts,
					const Function<dim>& control_function);
	void 	update_source_vector(const std::vector<Point<dim> >& locations, 
										const std::vector<double>& amounts); // done
	void 	solve(); // done
};


/** IMPLEMENTATION
*/
//-----------------------------------------------------------------------------
template<int dim>
FE_Chemical<dim>::FE_Chemical()
	:
	chemical_base(NULL)
	// use_bdf2_scheme(false)
	// viscosity_beta( 0. /* 0.017*dim */ ) // need to run tests to figure 
{}

template<int dim> //const Chemical_FE_Base<dim>&
FE_Chemical<dim>::FE_Chemical(std::shared_ptr<Chemical_FE_Base<dim> >	 chem_base, 
	double diffusion, double decay, double b, double dt)
	:
	chemical_base(chem_base),
	diffusion_constant(diffusion),
	decay_constant(decay),
	viscosity_beta(b),
	time_step(dt)
{
	reinit();
}

// INTERFACE:
// --------------------------------------------------------------------------------
template<int dim>
double 
FE_Chemical<dim>::value(const Point<dim>& p) const
{
	// return VectorTools::point_value(chemical_base->get_dof_handler(), solution, p);

    // Vector<double> value(1);
    // FETools::vector_point_value(chemical_base->getPointCellMap(), 
    // 							chemical_base->get_dof_handler(),
    // 	 						solution, 
    // 	 						p, 
    // 	 						value);

    // return value(0);

    return dealii::VectorTools::point_value(chemical_base->get_dof_handler(),
    										solution,
    										p);
}

template<int dim>
std::vector<double> 
FE_Chemical<dim>::value_list(const std::vector<Point<dim> >& points) const
{
	std::vector<double> result;
	result.reserve(points.size());

	for(unsigned int i = 0; i < points.size(); ++i)
		result.emplace_back(this->value(points[i]));

	return result;
}
	
template<int dim>
void 
FE_Chemical<dim>::project_function(const Function<dim>& initial_condition)
{
	project_initial_condition(initial_condition);	
}



template<int dim>
void 
FE_Chemical<dim>::update()
{
	update_rhs();
	solve();

	old_solution = solution;
}

template<int dim>
void 
FE_Chemical<dim>::update(const std::vector<Point<dim> >& locations, 
					const std::vector<double>& amounts)
{
	update_rhs(locations, amounts);
	solve();

	old_solution = solution;
}

template<int dim>
void 
FE_Chemical<dim>::update(const std::vector<Point<dim> >& locations, 
						const std::vector<double>& amounts,
						const Function<dim>& control_function)
{
	update_rhs(locations, amounts, control_function);
	solve();

	old_solution = solution;	
}


// ACCESSORS:
// --------------------------------------------------------------------------------
template<int dim>
double 
FE_Chemical<dim>::getDiffusionConstant() const
{
	return diffusion_constant;
}


template<int dim>
double 
FE_Chemical<dim>::getDecayConstant() const
{
	return decay_constant;
}

template<int dim>
double 
FE_Chemical<dim>::getViscosityBeta() const
{
	return viscosity_beta;
}

template<int dim>
double 	
FE_Chemical<dim>::getTimeStep() const
{
	return time_step;
}


// INTIALIZATION:
// --------------------------------------------------------------------------------------
template<int dim>
void 
FE_Chemical<dim>::reinit(const Chemical_FE_Base<dim>& chem_base, 
	double diffusion, double decay, double b, double dt)
{
	chemical_base = &chem_base;
	diffusion_constant = diffusion;
	decay_constant = decay;
	viscosity_beta = b;
	time_step = dt;

	reinit();
}


template<int dim>
void 
FE_Chemical<dim>::reinit()
{
	system_matrix.reinit(chemical_base->getSparsityPattern());

	solution.reinit(chemical_base->get_n_dofs());
	old_solution.reinit(chemical_base->get_n_dofs());

	right_hand_side.reinit(chemical_base->get_n_dofs());
	source.reinit(chemical_base->get_n_dofs());

	control_function_fem_vector.reinit(chemical_base->get_n_dofs());

	setup_solver(); 
}


template<int dim>
void 
FE_Chemical<dim>::project_initial_condition(const Function<dim>& initial_condition)
{
	VectorTools::project(chemical_base->get_dof_handler(),
						chemical_base->get_constraints(),
						QGauss<dim>(chemical_base->get_fe_degree() + 2 ),
						initial_condition,
						old_solution);
	solution = old_solution;
}


// UPDATE AND SOLVE:
// -------------------------------------------------------------------------------------------
template<int dim>
void	
FE_Chemical<dim>::update_rhs()
{
	chemical_base->get_mass_matrix().vmult(right_hand_side, old_solution); // RHS = MU^{k-1}
}

template<int dim>
void	
FE_Chemical<dim>::update_rhs(const std::vector<Point<dim> >& locations, 
				const std::vector<double>& amounts)
{
	update_rhs(); // RHS = MU^{k-1}
	update_source_vector(locations, amounts);
	right_hand_side.add(time_step, source); // RHS += k*S
}

template<int dim>
void	
FE_Chemical<dim>::update_rhs(const std::vector<Point<dim> >& locations, 
				const std::vector<double>& amounts,
				const Function<dim>& control_function)
{
	update_rhs(); // RHS = MU^{k-1}
	// ***FEM SCHEME FOR CONTROL FUNCTION: add control to right hand side
	update_source_vector(locations, amounts);
	right_hand_side.add(time_step, source); // RHS += k*S
	// project onto FEM vector:
	VectorTools::project(chemical_base->get_dof_handler(),
					chemical_base->get_constraints(),
					QGauss<dim>(chemical_base->get_fe_degree() + 2 ),
					control_function,
					control_function_fem_vector);

	// control_function_fem_vector.print(std::cout); 
	// add control to rhs:
	right_hand_side.add(time_step, control_function_fem_vector);
}


template<int dim>
void
FE_Chemical<dim>::update_source_vector(const std::vector<Point<dim> >& locations, 
										const std::vector<double>& amounts)
{
	source = 0;

	const double scale_factor = 0.0016;
	const unsigned int number_bacteria = locations.size();

	const unsigned int dofs_per_cell = chemical_base->get_fe().dofs_per_cell;
	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

	// loop over bacteria:
	for(unsigned int i = 0; i < number_bacteria; ++i)
	{
		std::pair<typename DoFHandler<dim>::active_cell_iterator, Point<dim> >
		cell_point = 
			chemical_base->getPointCellMap()->get_cell_point_pair(locations[i]);

		Quadrature<dim> q(GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

		FEValues<dim> fe_values(StaticMappingQ1<dim>::mapping, 
				                chemical_base->get_dof_handler().get_fe(),
				                q, 
				                UpdateFlags(update_values));

		fe_values.reinit(cell_point.first);

		cell_point.first->get_dof_indices (local_dof_indices);

		for (unsigned int j=0; j<dofs_per_cell; ++j)
			source(local_dof_indices[j]) += //check ***
				scale_factor*amounts[i]*fe_values.shape_value(j,0);
		// chemical_base->get_constraints().distribute_local_to_global(cell_matrix,
		// 					                                local_dof_indices,
		// 					                                system_matrix);

	}
}


template<int dim>
void
FE_Chemical<dim>::setup_solver()
{
	std::cout << "\n... assembling chemical system matrix " << std::endl;
	
	system_matrix = 0;

	const QGauss<dim> quadrature_formula(chemical_base->get_fe_degree() + 2);

	FEValues<dim>	chemical_fe_values(chemical_base->get_fe(),
			                            quadrature_formula,                                           
			                            update_values    |
			                            update_gradients |
			                            // update_hessians  |
			                            update_quadrature_points  |
			                            update_JxW_values);

	const QGauss<dim-1> face_quad(chemical_base->get_fe_degree() + 2);
	// FEFaceValues<dim> fe_face_values(chemical_base->get_fe(), face_quad,
	//         update_values | update_normal_vectors | update_quadrature_points | update_JxW_values);
	// const unsigned int n_face_q_points = face_quad.size();


	const unsigned int dofs_per_cell = chemical_base->get_fe().dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();

	FullMatrix<double> 	cell_matrix(dofs_per_cell, dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	// shape functions:
	std::vector<double> 			phi(dofs_per_cell);
	std::vector<Tensor<1, dim> >	grad_phi(dofs_per_cell);

	//typename DoFHandler<dim>::active_cell_iterator 
	auto cell         = chemical_base->get_dof_handler().begin_active();
	const auto endc   = chemical_base->get_dof_handler().end();
	unsigned int cell_id = 0; // for velocity access
	// auto stokes_cell  = stokes_dof_handler.begin_active();

	for(; cell!=endc; ++cell , ++cell_id)
	{
		cell_matrix = 0;

		chemical_fe_values.reinit(cell);

		// local artificial viscosity: // could treat advection explicitly ... need ``RHS matrix''
		const double nu = compute_viscosity(
							chemical_base->get_cell_velocity_values(cell_id),
							cell->diameter());

		// loop over quadrature points:
		for(unsigned int q = 0; q < n_q_points; ++q)
		{
			for(unsigned int k = 0; k < dofs_per_cell; ++k)
			{
				grad_phi[k] = chemical_fe_values.shape_grad(k,q);
				phi[k] = chemical_fe_values.shape_value(k,q);
			}

			// velocity values:
			const Tensor<1, dim> velocity = 
					chemical_base->get_cell_quad_velocity_value(cell_id, q);

			for(unsigned int i = 0; i < dofs_per_cell; ++i)
			{
				for(unsigned int j = 0; j < dofs_per_cell; ++j)
				{
					cell_matrix(i,j) += (
						// mass terms:
						(1. + decay_constant*time_step)*phi[i]*phi[j]
				
						// diffusion + artificial viscosity:
						+ grad_phi[i]*( (diffusion_constant + nu)*time_step )*grad_phi[j]
				
						// advection: ( maybe see what happens if integrating by parts here ...)
						+ phi[i]*time_step*velocity*grad_phi[j]

						// robin boundary (no flux) terms:
						// +  ...
						// left-off-here
						)*chemical_fe_values.JxW(q);									
				} // for matrix j
			} // for matrix i

		} // for quadrature points

		// // face terms for robin boundary conditions:
		// for(unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f)
		// {
		// 	if( (cell ->face(f)->boundary_id() == 0) ) // // can do all but outlet ...|| (cell ->face(f)->boundary_id() == 6) )
		// 	{
		// 		// std::cout << " at boundary 5" << std::endl;
		// 		fe_face_values.reinit(cell,f);

		// 		for(unsigned int q_index = 0; q_index < n_face_q_points; ++q_index)
		// 		{
		// 			const Tensor<1,2> velVal; // = chemical_base->getStokesSolution()->value(fe_face_values.quadrature_point(q_index)); // advection field value *****
		// 			const double neumann_value  = 1.0*(velVal * fe_face_values.normal_vector(q_index));

		// 			// if( neumann_value > 0)
		// 			  std::cout << "Boundary " << cell->face(f)->boundary_id()
		// 			    << " Neumann value: " << neumann_value << std::endl;
		// 			// PRINT VALUE AND LOCATION TO FILE!!!

		// 			for(unsigned int i=0; i < dofs_per_cell; ++i)
		// 			{
		// 				for(unsigned int j=0; j < dofs_per_cell; ++j)
		// 				{
		// 					// robin_
		// 					cell_matrix(i,j) += ( 
		// 							fe_face_values.shape_value(i,q_index) *
		// 							neumann_value * fe_face_values.shape_value(j, q_index) 
		// 						)* fe_face_values.JxW(q_index);         // may need a minus sign *****
		// 				}
		// 			} // for cell indicies
		// 		} // for quadrature points
		// 	} // if at a circle boundary
		// } // for cell faces (for boundary condit

		cell->get_dof_indices(local_dof_indices);
		chemical_base->get_constraints().distribute_local_to_global(cell_matrix,
							                                local_dof_indices,
							                                system_matrix);
	} // for cells

	// intialize direct solver:
	A_direct.initialize(system_matrix);
} // setup_solver()


template<int dim>
double 
FE_Chemical<dim>::compute_viscosity(const std::vector<Tensor<1,dim> >& velocity_values,
            						const double cell_diameter)
{
	double max_velocity = 0;

	const unsigned int n_q_points = velocity_values.size();

	for(unsigned int q = 0; q < n_q_points; ++q)
	{
		const Tensor<1,dim> u = velocity_values[q];

		max_velocity = std::max(std::sqrt(u*u), max_velocity);
	}

	// nu = beta ||w(x)|| h(x)
	return viscosity_beta*max_velocity*cell_diameter;
}


template<int dim>
void 
FE_Chemical<dim>::solve()
{
	A_direct.vmult(solution, right_hand_side);
    chemical_base->get_constraints().distribute(solution);
}


// MORE ACCESSORS (For debugging):
// ---------------------------------------------------------------------------------------
template<int dim>
double 
FE_Chemical<dim>::getMass() const
{
	return solution*chemical_base->get_mass_check_vector();
}


template<int dim>
double 
FE_Chemical<dim>::getMin() const
{
	double min_value = solution[0];

	for(unsigned int i = 0; i < solution.size(); ++i)
		min_value = std::min(min_value, solution[i]);

	return min_value;
}


template<int dim>
double 
FE_Chemical<dim>::getMax() const
{
	double max_value = solution[0];

	for(unsigned int i = 0; i < solution.size(); ++i)
		max_value = std::max(max_value, solution[i]);

	return max_value;
}


template<int dim>
void 
FE_Chemical<dim>::output_solution(const std::string& output_directory,
						 const unsigned int chem_id,
						 const unsigned int save_step) const
{
	DataOut<dim> data_out;
	data_out.attach_dof_handler(chemical_base->get_dof_handler());

	std::string chemical_name = "Chemical" + Utilities::int_to_string(chem_id,3);

	data_out.add_data_vector(solution, chemical_name);
	data_out.build_patches();
	const std::string filename = output_directory
				  + "/"
	              + chemical_name
	              + "_"
	              + Utilities::int_to_string(save_step,4)
	              + ".vtk";
	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
}


template<int dim>
void 
FE_Chemical<dim>::output_solution(const std::string& output_directory,
						 const unsigned int chem_id,
						 const unsigned int save_step,
						 const unsigned int cycle) const
{
	DataOut<dim> data_out;
	data_out.attach_dof_handler(chemical_base->get_dof_handler());

	std::string chemical_name = "Chemical" 
				+ Utilities::int_to_string(chem_id,3)
				+ "_"
	            + Utilities::int_to_string(cycle,2);

	data_out.add_data_vector(solution, chemical_name);
	data_out.build_patches();
	const std::string filename = output_directory
				  + "/"
	              + chemical_name
	              + "-"
	              + Utilities::int_to_string(save_step,4)
	              + ".vtk";
	std::ofstream output(filename.c_str());
	data_out.write_vtk(output);
}

template<int dim>
void 	
FE_Chemical<dim>::print(std::ostream& out, unsigned int chem_id) const
{
	DataOut<dim> data_out;
	data_out.attach_dof_handler(chemical_base->get_dof_handler());

	std::string chemical_name = "Chemical" 
				+ Utilities::int_to_string(chem_id,3);

	data_out.add_data_vector(solution, chemical_name);
	data_out.build_patches();

	data_out.write_vtk(out);	
}

template<int dim>
void 	
FE_Chemical<dim>::printInfo(std::ostream& out) const
{
	out << "\n\n-----------------------------------------------------" << std::endl
		<< "\t\t FE CHEMICAL INFO:"
		<< "\n-----------------------------------------------------" << std::endl;

   	out << "Diffusion constant: " << diffusion_constant << std::endl
   		<< "Decay constant: " << decay_constant << std::endl
   		<< "Viscosity beta: " << viscosity_beta << std::endl
   		<< "Time step: " << time_step << std::endl;

	out << "\n-----------------------------------------------------\n\n" << std::endl;
}

}} // close namespace
#endif