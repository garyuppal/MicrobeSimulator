#ifndef MICROBESIMULATOR_VELOCITY_INTERFACE_H
#define MICROBESIMULATOR_VELOCITY_INTERFACE_H

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

// for fe values:
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
// quadrature formula:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe.h>

namespace MicrobeSimulator{ namespace Velocity{
	using namespace dealii;

template<int dim>
class VelocityInterface{
public:
	VelocityInterface() {}
	virtual ~VelocityInterface() {}

	// pure virtual methods:
	virtual Tensor<1, dim> value(const Point<dim>& location) const = 0;
	virtual double get_maximum_velocity(double max_coordinate) const = 0;

	// virtual with default implementation:
	virtual std::vector<std::vector<Tensor<1, dim> > >
		get_fe_velocity_values(unsigned int fe_degree, 
			const FiniteElement<dim>& fe,
			const DoFHandler<dim>& dof) const;

	// value list: (shouldn't need to be overridden):
	void value_list(const std::vector<Point<dim> >& points,
					std::vector<Tensor<1, dim> >& values) const;

	virtual void printInfo(std::ostream& out) const = 0;
};

// default implementation for value_list
template<int dim>
void
VelocityInterface<dim>::value_list(const std::vector<Point<dim> >& points,
					std::vector<Tensor<1, dim> >& values) const
{
	if( values.size() != points.size() )
		values.resize( points.size() );

	for(unsigned int i = 0; i < points.size(); ++i)
		values[i] = this->value(points[i]);
}

// default implementation (mainly used for all classes except stokes_solver):
template<int dim>
std::vector<std::vector<Tensor<1, dim> > >
VelocityInterface<dim>::get_fe_velocity_values(unsigned int fe_degree, 
			const FiniteElement<dim>& fe,
			const DoFHandler<dim>& dof) const 
{
	const QGauss<dim> quadrature_formula(fe_degree + 2);
	FEValues<dim>  fe_values(fe, quadrature_formula,
	      update_values | update_quadrature_points | update_gradients | update_JxW_values);

	const unsigned int n_q_points = quadrature_formula.size();

	std::cout << "n_q_points: " << n_q_points << std::endl;

	std::vector<std::vector<dealii::Tensor<1, dim> > > velocity_values;
	velocity_values.reserve(dof.n_dofs());
	std::vector<Tensor<1, dim> > 		cell_velocity_values(n_q_points); 
				//, Tensor<1,dim>()); // temporary ''fix''

	// std::cout << "STILL NEED TO FIX THIS" << std::endl;

	auto cell 			= dof.begin_active();
	const auto endc 	= dof.end();
	for(; cell != endc; ++cell)
	{
		fe_values.reinit(cell);

		this->value_list(fe_values.get_quadrature_points(), cell_velocity_values);

		// std::cout << "QUADRATURE POINTS: (" 
		// 	<< fe_values.get_quadrature_points().size() << ") " << std::endl;
		// for(unsigned int i = 0; i < fe_values.get_quadrature_points().size(); ++i)
		// 	std::cout << fe_values.get_quadrature_points()[i] << " ";
		// std::cout << std::endl;

		velocity_values.emplace_back(cell_velocity_values);
	}

	// std::cout << "Velocity values for (" << velocity_values.size()
	// 	<< ") cells" << std::endl << std::endl;

	return velocity_values;
}

}}
#endif