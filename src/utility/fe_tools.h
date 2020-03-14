#ifndef MICROBESIMULATOR_FE_TOOLS_H
#define MICROBESIMULATOR_FE_TOOLS_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

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

#include "../utility/cell_iterator_map.h"

#include <vector>

// probably dont actually need all above headers

namespace MicrobeSimulator{ namespace FETools{
	using namespace dealii;

// @ note: may want to include spacedim for later 3d use ...
template<int dim, typename VectorType>
void
vector_point_value(const PointCellMap<dim>* 		point_cell_map,
				const DoFHandler<dim> &    	  		dof,
				const VectorType &                	fe_function,
				const Point<dim> &                	point,
				Vector<typename VectorType::value_type> &value)
{
	using Number = typename VectorType::value_type;
	const FiniteElement<dim> &fe = dof.get_fe();

	std::pair<typename DoFHandler<dim>::active_cell_iterator, Point<dim> >
	cell_point = 
		point_cell_map->get_cell_point_pair(point);

	// rest from deal ii code for point_value:

	const Quadrature<dim> quadrature(GeometryInfo<dim>::project_to_unit_cell(
															cell_point.second));

	FEValues<dim> fe_values(StaticMappingQ1<dim>::mapping,
										 fe, quadrature, update_values);

	fe_values.reinit(cell_point.first);

	std::vector<Vector<Number> > u_value(1, Vector<Number>(fe.n_components()));
	fe_values.get_function_values(fe_function, u_value);

	value = u_value[0];
}


template<int dim, typename VectorType>
void
vector_point_value(std::shared_ptr<FETools::PointCellMap<dim> >		point_cell_map,
				const DoFHandler<dim> &    	  		dof,
				const VectorType &                	fe_function,
				const Point<dim> &                	point,
				Vector<typename VectorType::value_type> &value)
{
	using Number = typename VectorType::value_type;
	const FiniteElement<dim> &fe = dof.get_fe();

	std::pair<typename DoFHandler<dim>::active_cell_iterator, Point<dim> >
	cell_point = 
		point_cell_map->get_cell_point_pair(point);

	// rest from deal ii code for point_value:

	const Quadrature<dim> quadrature(GeometryInfo<dim>::project_to_unit_cell(
															cell_point.second));

	FEValues<dim> fe_values(StaticMappingQ1<dim>::mapping,
										 fe, quadrature, update_values);

	fe_values.reinit(cell_point.first);

	std::vector<Vector<Number> > u_value(1, Vector<Number>(fe.n_components()));
	fe_values.get_function_values(fe_function, u_value);

	value = u_value[0];
}


}} // close namespace
#endif