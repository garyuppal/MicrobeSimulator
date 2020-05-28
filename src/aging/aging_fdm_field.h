#pragma once

#include "../refactored_chemicals/field_base.h"
#include "../utility/utility.h"  // for BoundaryCondition enum

#include <array>

// TO DO Need to also add ghost cells.... (increase bounds by +/- dx?), should be a better way..., or could have two sets of bounds
// with fdm_chem keeping ''real'' bounds and the underlying fields with a ghost buffer ....
namespace MicrobeSimulator{ namespace Aging{

template<int dim, typename Number = double>
class Aging_FDM_Field : public RefactoredChemicals::FieldBase<dim,Number>{
public:
	Aging_FDM_Field();

	// Aging_FDM_Field(const Point<dim>& lower, 
	// 		const Point<dim>& upper,
	// 		const std::array<unsigned int, dim>& disc, 
	// 		const Number& init_value);

	void init(const Point<dim>& lower, 
			const Point<dim>& upper,
			const std::array<unsigned int, dim>& disc, 
			const std::array<BoundaryCondition, dim>& bcs,
			const Number& init_value);


	// inherited:
	// virtual void reinit(const Point<dim>& lower, const Point<dim> upper,
	// 	const std::array<unsigned int, dim>& disc, const T& init_value);

	// convert read access to interpolation:
	// virtual T at(const Point<dim> point) const;

	// read access with boundary conditions:
	double at_bc(int i) const;
	double at_bc(int i, int j) const;
	double at_bc(int i, int j, int k) const;

	void print(std::ostream& out) const override;

	double getInverseCellWidthSquared(unsigned int dimension) const;

private:
	std::array<double, dim> invsqr_d_cells;	
	// store locally for easier diffusion implementation

	std::array<BoundaryCondition, dim> boundary_conditions;

	void 	setInverseSquareCellWidth();
};

// IMPL
// -------------------------------------------------------------------------------
template<int dim, typename Number>
Aging_FDM_Field<dim,Number>::Aging_FDM_Field()
{}

// template<int dim, typename Number>
// Aging_FDM_Field<dim, Number>::Aging_FDM_Field(const Point<dim>& lower, const Point<dim>& upper,
// 		const std::array<unsigned int, dim>& disc, const Number& init_value)
// 	:
// 	FieldBase<dim,Number>(lower, upper, disc, init_value)
// {
// 	setInverseSquareCellWidth();	
// }

template<int dim, typename Number>
void 
Aging_FDM_Field<dim, Number>::init(const Point<dim>& lower, const Point<dim>& upper,
		const std::array<unsigned int, dim>& disc, const std::array<BoundaryCondition, dim>& bcs, 
		const Number& init_value)
{
	this->reinit(lower, upper, disc, init_value);
	setInverseSquareCellWidth();
	boundary_conditions = bcs;
}

template<int dim, typename Number>
void 	
Aging_FDM_Field<dim,Number>::setInverseSquareCellWidth()
{
	for(unsigned int i = 0; i < dim; ++i)
		invsqr_d_cells[i] = std::pow(this->inv_d_cells[i],2);
}

template<int dim, typename Number>
double 
Aging_FDM_Field<dim,Number>::getInverseCellWidthSquared(unsigned int dimension) const
{
	return invsqr_d_cells[dimension];
}


// read access:
// ----------------------------------------------------------------------------

template<>
double 
Aging_FDM_Field<1, double>::at_bc(int i) const
{
	if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::WRAP))
		i = (int)(this->n_cells)[0]-1;
	else if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
		i = 0;
	else if( (i == (int)(this->n_cells)[0]) && (boundary_conditions[0] == BoundaryCondition::WRAP))
		i = 0;
	else if( (i == (int)(this->n_cells)[0]) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
		i = (int)(this->n_cells)[0]-1;

	return field[i];
}

template<>
double 
Aging_FDM_Field<2, double>::at_bc(int /*i*/) const
{
	std::cout << "not implemented" << std::endl;
	assert(false);
}

template<>
double 
Aging_FDM_Field<3, double>::at_bc(int /*i*/) const
{
	std::cout << "not implemented" << std::endl;
	assert(false);
}

template<>
double 
Aging_FDM_Field<1, double>::at_bc(int /*i*/, int /*j*/) const
{
	std::cout << "not implemented" << std::endl;
	assert(false);
}

template<>
double 
Aging_FDM_Field<2, double>::at_bc(int i, int j) const
{
	if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::WRAP))
		i = (int)(this->n_cells)[0]-1;
	else if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
		i = 0;
	else if( (i == (int)(this->n_cells)[0]) && (boundary_conditions[0] == BoundaryCondition::WRAP))
		i = 0;
	else if( (i == (int)(this->n_cells)[0]) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
		i = (int)(this->n_cells)[0]-1;

	if( (j == -1) && (boundary_conditions[1] == BoundaryCondition::WRAP))
		j = (int)(this->n_cells)[1]-1;
	else if( (j == -1) && (boundary_conditions[1] == BoundaryCondition::REFLECT))
		j = 0;
	else if( (j == (int)(this->n_cells)[1]) && (boundary_conditions[1] == BoundaryCondition::WRAP))
		j = 0;
	else if( (j == (int)(this->n_cells)[1]) && (boundary_conditions[1] == BoundaryCondition::REFLECT))
		j = (int)(this->n_cells)[1]-1;

	return field[ i + j*(this->n_cells)[0] ];
}

template<>
double 
Aging_FDM_Field<3, double>::at_bc(int /*i*/, int /*j*/) const
{
	std::cout << "not implemented" << std::endl;
	assert(false);
}


template<>
double
Aging_FDM_Field<1, double>::at_bc(int /*i*/, int /*j*/, int /*k*/) const
{
	std::cout << "not implemented" << std::endl;
	assert(false);
}

template<>
double
Aging_FDM_Field<2, double>::at_bc(int /*i*/, int /*j*/, int /*k*/) const
{
	std::cout << "not implemented" << std::endl;
	assert(false);
}

template<>
double
Aging_FDM_Field<3, double>::at_bc(int i, int j, int k) const
{
	if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::WRAP))
		i = (this->n_cells)[0]-1;
	else if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
		i = 0;
	else if( (i == (int)(this->n_cells)[0]) && (boundary_conditions[0] == BoundaryCondition::WRAP))
		i = 0;
	else if( (i == (int)(this->n_cells)[0]) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
		i = (int)(this->n_cells)[0]-1;

	if( (j == -1) && (boundary_conditions[1] == BoundaryCondition::WRAP))
		j = (int)(this->n_cells)[1]-1;
	else if( (j == -1) && (boundary_conditions[1] == BoundaryCondition::REFLECT))
		j = 0;
	else if( (j == (int)(this->n_cells)[1]) && (boundary_conditions[1] == BoundaryCondition::WRAP))
		j = 0;
	else if( (j == (int)(this->n_cells)[1]) && (boundary_conditions[1] == BoundaryCondition::REFLECT))
		j = (int)(this->n_cells)[1]-1;

	if( (k == -1) && (boundary_conditions[2] == BoundaryCondition::WRAP))
		k = (int)(this->n_cells)[2]-1;
	else if( (k == -1) && (boundary_conditions[2] == BoundaryCondition::REFLECT))
		k = 0;
	else if( (k == (int)(this->n_cells)[2]) && (boundary_conditions[2] == BoundaryCondition::WRAP))
		k = 0;
	else if( (k == (int)(this->n_cells)[2]) && (boundary_conditions[2] == BoundaryCondition::REFLECT))
		k = (int)(this->n_cells)[2]-1;


	return field[ i + j*(this->n_cells)[0] + k*(this->n_cells)[0]*(this->n_cells)[1] ];
}

// print methods:
// ----------------------------------------------------------------------------

/** \brief Specialized field print for 1D field of doubles */
template<>
void 
Aging_FDM_Field<1, double>::print(std::ostream& out) const
{
	for(unsigned int i = 0; i < this->field.size(); ++i)
		out << this->field[i] << " ";
	out << std::endl;
}

/** \brief Specialized field print for 2D field of doubles */
template<>
void 
Aging_FDM_Field<2, double>::print(std::ostream& out) const
{
	for(unsigned int j = 0; j < this->n_cells[1]; ++j)
	{
		for(unsigned int i = 0; i < this->n_cells[0]; ++i)
		{
			out << this->at(i,j) << " ";
		}
		out << std::endl;
	}
}

/** \brief Specialized field print for 3D field of doubles */
template<>
void 
Aging_FDM_Field<3, double>::print(std::ostream& out) const
{
	for(unsigned int k = 0; k < this->n_cells[2]; ++k)
	{
		for(unsigned int j = 0; j < this->n_cells[1]; ++j)
		{
			for(unsigned int i = 0; i < this->n_cells[0]; ++i)
			{
				out << this->at(i,j,k) << " ";
			}
			out << std::endl;
		}
		out << std::endl;
	} 
}

/** \brief Default to print all out */
template<int dim, typename Number>
void
Aging_FDM_Field<dim, Number>::print(std::ostream& out) const
{
	for(unsigned int i = 0; i < this->field.size(); ++i)
		out << this->field[i] << " ";
	out << std::endl;
}

}} // close namespaces
/* aging_Aging_FDM_Field.h */