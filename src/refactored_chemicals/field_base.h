#ifndef MICROBESIMULATOR_REFACTORED_FIELD_BASE_H
#define MICROBESIMULATOR_REFACTORED_FIELD_BASE_H

#include <deal.II/base/point.h>
using dealii::Point;

#include "../utility/utility.h"

#include <iostream>
#include <vector>
#include <array>


namespace MicrobeSimulator{ namespace RefactoredChemicals{

/** \brief basic field class
* used most generally as a spatial multi-dimenstional container
* fields that need calculus or other operations derive from this class
* another example is the used of cell-linked-lists where we store a field
* of lists of neighbors of interacting particles
*/
template<int dim, class T>
class FieldBase{
public:
	FieldBase();
	virtual ~FieldBase() {}

	// other constructors:
	FieldBase(const Point<dim>& lower, const Point<dim>& upper,
		const std::array<unsigned int, dim>& disc, const T& init_value);
	// FieldBase(const FDM_Field<dim,Number>& df);

	// re-initialization:
	virtual void reinit(const Point<dim>& lower, const Point<dim> upper,
		const std::array<unsigned int, dim>& disc, const T& init_value); // give initial value T ()

	// read access by grid location:
	const T& at(unsigned int i) const; // one-dimensional field
	const T& at(unsigned int i, unsigned int j) const; // two-dimensional field
	const T& at(unsigned int i, unsigned int j, unsigned int k) const; // three-dimensional field
	// write access
	T& at(unsigned int i); 
	T& at(unsigned int i, unsigned int j); 
	T& at(unsigned int i, unsigned int j, unsigned int k); 

	// read access by points:
	virtual T at(const Point<dim> point) const;
		// option to change to interpolation for numerical fields
	// write access by points (find containing cell):
	T& at(const Point<dim> point); 

	// read only direct access:
	const T& operator[](unsigned int i) const;
	T& operator[](unsigned int i);

	unsigned int getSize() const;

	// OTHER ACCESSORS:
	unsigned int getDiscretization(unsigned int dimension) const;
	double getCellWidth(unsigned int dimension) const;
	double getInverseCellWidth(unsigned int dimension) const;

	Point<dim> getCellCenterPoint(unsigned int i) const;
	Point<dim> getCellCenterPoint(unsigned int i, unsigned int j) const;
	Point<dim> getCellCenterPoint(unsigned int i, unsigned int j,unsigned int k) const;

	// MUTATORS:
	// void setBounds(const Point<dim>& lower, const Point<dim> upper); 
		// would need to modify discretization and possbily data

	// move to fdm_field
	// void setBoundaryConditions(const std::array<BoundaryCondition, dim>& bcs);

	virtual void print(std::ostream& out) const;

	unsigned int indexFromPoint(const Point<dim>& point) const;
protected:
	std::vector<T> field; // is there a better structure to use than vector?

	// spatial boundaries:
	Point<dim> bottom_left;
	Point<dim> top_right;

	// discretization:
	// these are not all independent but we store them for easy access
	// constructors and intializers are responsible for making these consistent
	// and initializing one from the other
	std::array<unsigned int, dim> n_cells; // nx, ny, nz
	std::array<double, dim> d_cells; // dx, dy, dz  
	std::array<double, dim> inv_d_cells; // (1/dx), etc...

	// std::array<BoundaryCondition, dim> boundary_conditions;

	unsigned int total_size;

	// initialization helpers:
	void setCellWidths();
	void setTotalSize();
	void init_field(const T& init_value);
};


// IMPL
// --------------------------------------------------------------------
template<int dim, class T>
FieldBase<dim,T>::FieldBase()
	:
	total_size(0)
{}

template<int dim, class T>
void 
FieldBase<dim,T>::setCellWidths()
{
	for(unsigned int i = 0; i < dim; ++i){
		d_cells[i] = (top_right - bottom_left)[i]/((double)n_cells[i]);
		inv_d_cells[i] = 1./d_cells[i];
	}
}

template<int dim, class T>
void
FieldBase<dim,T>::setTotalSize()
{
	total_size = 1;
	for(unsigned int i = 0; i < dim; ++i)
		total_size *= n_cells[i];
}

template<int dim, class T>
void 
FieldBase<dim,T>::init_field(const T& init_value)
{
	setTotalSize();
	field.assign(total_size, init_value);
}

template<int dim, class T>
unsigned int 
FieldBase<dim,T>::indexFromPoint(const Point<dim>& point) const
{
	unsigned int i = 0;
	unsigned int discProd = 1;

	for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
	{
		double reset_value = point[dim_itr] - bottom_left[dim_itr];

		assert(reset_value>=0); // may not happen due to round off errors

		i += discProd*floor( reset_value*inv_d_cells[dim_itr] );
		discProd *= n_cells[dim_itr];
	}

	assert(i < total_size);
	return i;
}

// --------------------------------------------------------------------

template<int dim, class T>
FieldBase<dim,T>::FieldBase(const Point<dim>& lower, const Point<dim>& upper,
	const std::array<unsigned int, dim>& disc, const T& init_value)
	:
	bottom_left(lower),
	top_right(upper),
	n_cells(disc),
	total_size(0)
{
	setCellWidths();
	init_field(init_value);
} 

// re-initialization:
template<int dim, class T>
void 
FieldBase<dim,T>::reinit(const Point<dim>& lower, const Point<dim> upper,
	const std::array<unsigned int, dim>& disc, const T& init_value)
{
	bottom_left = lower;
	top_right = upper;
	n_cells = disc;
	setCellWidths();
	init_field(init_value);
}

template<int dim, class T>
unsigned int 
FieldBase<dim, T>::getSize() const
{
	return field.size();
}

template<int dim, class T>
const T& 
FieldBase<dim, T>::operator[](unsigned int i) const
{
	return field[i];
}

template<int dim, class T>
T& 
FieldBase<dim, T>::operator[](unsigned int i)
{
	return field[i];
}

// read access by grid location:
template<int dim, class T>
const T& 
FieldBase<dim,T>::at(unsigned int i) const	
{
	assert(dim==1);

	return field[i];
}

template<int dim, class T>
const T& 
FieldBase<dim,T>::at(unsigned int i, unsigned int j) const
{
	assert(dim==2);

	return field[i + (j*n_cells[0])];
}

template<int dim, class T>
const T& 
FieldBase<dim,T>::at(unsigned int i, unsigned int j, unsigned int k) const
{
	assert(dim==3);

	return field[i + (j*n_cells[0]) + (k*n_cells[0]*n_cells[1])];
}

// write access
template<int dim, class T>
T& 
FieldBase<dim,T>::at(unsigned int i) 	
{
	assert(dim==1);

	return field[i];
}

template<int dim, class T>
T&
FieldBase<dim,T>::at(unsigned int i, unsigned int j) 
{
	assert(dim==2);

	return field[i + (j*n_cells[0])];
}

template<int dim, class T>
T&
FieldBase<dim,T>::at(unsigned int i, unsigned int j, unsigned int k) 
{
	assert(dim==3);

	return field[i + (j*n_cells[0]) + (k*n_cells[0]*n_cells[1])];
}

template<int dim, class T>
T 
FieldBase<dim,T>::at(const Point<dim> point) const 
{ 
	return field[indexFromPoint(point)];
}

template<int dim, class T>
T& 
FieldBase<dim,T>::at(const Point<dim> point)
{
	return field[indexFromPoint(point)];
}

template<int dim, class T>
unsigned int 
FieldBase<dim,T>::getDiscretization(unsigned int dimension) const
{
	return n_cells[dimension];
}

template<int dim, class T>
double 
FieldBase<dim,T>::getCellWidth(unsigned int dimension) const
{
	return d_cells[dimension];
}

template<int dim, class T>
double 
FieldBase<dim,T>::getInverseCellWidth(unsigned int dimension) const
{
	return inv_d_cells[dimension];
}

template<int dim, class T>
Point<dim> 
FieldBase<dim,T>::getCellCenterPoint(unsigned int i) const
{
	assert(dim==1);

	return Point<dim>(bottom_left[0] + (i+0.5)*getCellWidth(0));
}

template<int dim, class T>
Point<dim> 
FieldBase<dim,T>::getCellCenterPoint(unsigned int i, unsigned int j) const
{
	assert(dim==2);

	const double x = bottom_left[0] + (i+ 0.5)*getCellWidth(0);
	const double y = bottom_left[1] + (j + 0.5)*getCellWidth(1); 

	return Point<dim>(x,y);
}

template<int dim, class T>
Point<dim> 
FieldBase<dim,T>::getCellCenterPoint(unsigned int i, unsigned int j,unsigned int k) const
{
	assert(dim==3);

	const double x = bottom_left[0] + (i+ 0.5)*getCellWidth(0);
	const double y = bottom_left[1] + (j + 0.5)*getCellWidth(1); 
	const double z = bottom_left[2] + (k + 0.5)*getCellWidth(2);

	return Point<dim>(x,y,z);	
}

// template<int dim, class T>
// void 
// FieldBase<dim,T>::setBoundaryConditions(const std::array<BoundaryCondition, dim>& bcs)
// {
// 	boundary_conditions = bcs;
// }

template<int dim, class T>
void 
FieldBase<dim,T>::print(std::ostream& out) const
{
	out << "Field Base: can override for specific field type" << std::endl;
}

}} // CLOSE NAMESPACE
#endif