#ifndef MATH_FDM_FIELD_H
#define MATH_FDM_FIELD_H

// #include "./function.h"
#include "./field_base.h"

namespace Chemotaxis{ namespace Math{

template<int dim, typename Number = double>
class FDM_Field : public FieldBase<dim,Number>{
public:
	FDM_Field();

	FDM_Field(const Vect<dim>& lower, const Vect<dim>& upper,
		const std::array<unsigned int, dim>& disc, const Number& init_value);
	// virtual void reinit(const Vect<dim>& lower, const Vect<dim> upper,
	// 	const std::array<unsigned int, dim>& disc, const T& init_value) override;

	// convert read access to interpolation:
	// virtual T at(const Vect<dim> point) const;

	virtual void print(std::ostream& out) const;

	double getInverseCellWidthSquared(unsigned int dimension) const;

private:
	std::array<double, dim> invsqr_d_cells;	

	void 	setInverseSquareCellWidth();
};

// IMPL
// -------------------------------------------------------------------------------
template<int dim, typename Number>
FDM_Field<dim,Number>::FDM_Field()
{}

template<int dim, typename Number>
FDM_Field<dim, Number>::FDM_Field(const Vect<dim>& lower, const Vect<dim>& upper,
		const std::array<unsigned int, dim>& disc, const Number& init_value)
	:
	FieldBase<dim,Number>(lower, upper, disc, init_value)
{
	setInverseSquareCellWidth();	
}

template<int dim, typename Number>
void 	
FDM_Field<dim,Number>::setInverseSquareCellWidth()
{
	for(unsigned int i = 0; i < dim; ++i)
		invsqr_d_cells[i] = std::pow(this->inv_d_cells[i],2);
}
// ----------------------------------------------------------------------------

template<>
void 
FDM_Field<1, double>::print(std::ostream& out) const
{
	for(unsigned int i = 0; i < this->field.size(); ++i)
		out << this->field[i] << " ";
	out << std::endl;
}

template<>
void 
FDM_Field<2, double>::print(std::ostream& out) const
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

template<>
void 
FDM_Field<3, double>::print(std::ostream& out) const
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


}} // close namespaces
#endif