#ifndef MICROBESIMULATOR_DISCRETE_FIELD_H
#define MICROBESIMULATOR_DISCRETE_FIELD_H

#include <deal.II/base/point.h>
using dealii::Point;

#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "../utility/enum_types.h"
#include "./discrete_function.h"

/*
	WITH FILE INTIALIZATION, CAN SAVE CONFIGURATIONS OR RUNS, AND RESTART THEM
	LATER BY LOADING LAST SAVED SYSTEM!!!!
*/

namespace MicrobeSimulator{

	template<int dim>
	class DiscreteField{
	public:
		DiscreteField();
		DiscreteField(const Point<dim>& lower, const Point<dim>& upper,
			const std::array<unsigned int, dim>& disc);
		DiscreteField(const DiscreteField<dim>& df);

		// initialization:
		void initialize(const Point<dim>& lower, const Point<dim>& upper,
			const std::array<unsigned int, dim>& disc,
			const std::array<BoundaryCondition, dim>& bcs); // not all discrete fields have bcs??...
		void initialize(std::ifstream& infile, const Point<dim> &lower, 
			double velocity_scale, double geometry_scale);
		void reinitialize();

		// overload [] for regular indexing

		// accessors:
		unsigned int getDiscretization(unsigned int dimension) const; 
		BoundaryCondition getBoundaryCondition(unsigned int dimension) const;
		Point<dim> getLowerBound() const;
		Point<dim> getUpperBound() const;

		double getCellWidth(unsigned int dimension) const;
		double getInverseCellWidth(unsigned int dimension) const;
		double getInverseCellWidthSquared(unsigned int dimension) const;

		// dimension dependent accessors:
		double at(unsigned int i) const;
		double at(unsigned int i, unsigned int j) const;
		double at(unsigned int i, unsigned int j, unsigned int k) const; 

		// with boundary condition checking:
		double at_bc(int i) const;
		double at_bc(int i, int j) const;
		double at_bc(int i, int j, int k) const; 

		// read and write access:
		double& at(unsigned int i);
		double& at(unsigned int i, unsigned int j);
		double& at(unsigned int i, unsigned int j, unsigned int k); 

		// access by point:
		double at(const Point<dim>& p) const;
		double& at(const Point<dim>& p);


		// mutators: // Warning: these reinitialize fields to zero!!!
		void setBounds(const Point<dim>& lower, const Point<dim>& upper);
		void setDiscretization(const std::array<unsigned int, dim>& disc);	
		void setBoundaryConditions(const std::array<BoundaryCondition, dim>& bcs);

		// project function:
		void projectFunction(const DiscreteFunction<dim>& function);
		Point<dim> getCellCenterPoint(unsigned int i, unsigned int j) const;
		Point<dim> getCellCenterPoint(unsigned int i, unsigned int j,unsigned int k) const;


		void print(std::ostream& out) const;
		void printInfo(std::ostream& out) const;

	private:
		std::vector<double> field;

		Point<dim> bottom_left;
		Point<dim> top_right;
		std::array<unsigned int, dim> discretization;
		std::array<BoundaryCondition, dim> boundary_conditions;

		unsigned int indexFromPoint(const Point<dim>& p) const;
	};


//Implementation:
//-----------------------------------------------------------------

	template<int dim>
	DiscreteField<dim>::DiscreteField()
	{}

	template<int dim>
	DiscreteField<dim>::DiscreteField(const Point<dim>& lower, const Point<dim>& upper,
		const std::array<unsigned int, dim>& disc)
		:
		bottom_left(lower),
		top_right(upper),
		discretization(disc)
	{
		reinitialize();
	}

	template<int dim>
	DiscreteField<dim>::DiscreteField(const DiscreteField<dim>& df)
	{
		field = df.field;
		bottom_left = df.bottom_left;
		top_right = df.top_right;
		discretization = df.discretization;
		boundary_conditions = df.boundary_conditions;
	}

	template<int dim>
	void DiscreteField<dim>::initialize(const Point<dim>& lower, const Point<dim>& upper,
		const std::array<unsigned int, dim>& disc, 
		const std::array<BoundaryCondition, dim>& bcs)
	{
		setBounds(lower,upper);
		setDiscretization(disc);
		boundary_conditions = bcs;
	}

	template<int dim>
	void DiscreteField<dim>::initialize(std::ifstream& infile, const Point<dim> &lower,
		double scale, double geometry_scale)
	{
		if(dim != 2)
			throw std::invalid_argument("File initialization for field only implented for dim = 2 for now");

		std::array<unsigned int, dim> disc;
		double dx = 1;
		double dy = 1;

		  unsigned int ic = 0;
		  std::string line;

		  // FOR FIELD FILE:
		  while(std::getline(infile,line)){
		    if(ic < 4){
		      // read first four lines give Nx, Ny, dx, dy: 
		      std::istringstream stream(line);
		      std::string varName;
		      double value; 
		      stream >> varName >> value;  

		      // cases:
		      if(varName.compare("Nx") == 0){disc[0] = value; ic++;}
		      if(varName.compare("Ny") == 0){disc[1] = value; ic++;}
		      if(varName.compare("dx") == 0){dx = value; ic++;}
		      if(varName.compare("dy") == 0){dy = value; ic++;}
		    } // for first 4 lines
		    else{
		      break;
		    } // else
		  } // while reading file

		  // lower bound already set and scaled by geometry
		  // rescale resolution:
		  dx *= geometry_scale;
		  dy *= geometry_scale;

		  Point<dim> upper;

		  upper[0] = lower[0] + disc[0]*dx;
		  upper[1] = lower[1] + disc[1]*dy;

		  std::array<BoundaryCondition,dim> bcs = {BoundaryCondition::REFLECT,
		  	BoundaryCondition::REFLECT};

		  initialize(lower,upper,disc,bcs);

		  unsigned int i = 0;
		  while(std::getline(infile,line,','))
		  {
		    std::istringstream stream(line);
		    double value;
		    stream >> value;

		    field[i] =  scale*value;
		    ++i;
		  }

		  //@todo : ``fix'' this so that we actually get rows into x 
		  // and columns into y
	}

	template<int dim>
	void DiscreteField<dim>::reinitialize()
	{
		unsigned int size = 1;
		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
			size *= discretization[dim_itr];

		field.assign(size, 0);
	}


	// accessors:
	template<int dim>
	unsigned int DiscreteField<dim>::getDiscretization(unsigned int dimension) const
	{
		if(dimension >= dim)
			throw std::invalid_argument("Cannot get discretization for non-existent dimension");

		return discretization[dimension];
	}

	template<int dim>
	BoundaryCondition DiscreteField<dim>::getBoundaryCondition(unsigned int dimension) const
	{
		if(dimension >= dim)
			throw std::invalid_argument("Cannot get boundarycondition for non-existent dimension");

		return boundary_conditions[dimension];
	}

	template<int dim> 
	Point<dim> DiscreteField<dim>::getLowerBound() const
	{
		return bottom_left;
	}

	template<int dim>
	Point<dim> DiscreteField<dim>::getUpperBound() const
	{
		return top_right;
	}

	template<int dim>
	double DiscreteField<dim>::getCellWidth(unsigned int dimension) const
	{
		if(dimension >= dim)
			throw std::invalid_argument("Cannot get cell width for non-existent dimension");

		return (top_right[dimension] - bottom_left[dimension])/discretization[dimension];
	}

	template<int dim>
	double DiscreteField<dim>::getInverseCellWidth(unsigned int dimension) const
	{
		return 1.0/getCellWidth(dimension);
	}

	template<int dim>
	double DiscreteField<dim>::getInverseCellWidthSquared(unsigned int dimension) const
	{
		return getInverseCellWidth(dimension)*getInverseCellWidth(dimension);
	}

	// dimension dependent accessors:
	template<int dim>
	double DiscreteField<dim>::at(unsigned int i) const
	{
		if(dim != 1)
			throw std::invalid_argument("Field dimension must equal 1 to use this accessor");
		if(i >= field.size())
			throw std::invalid_argument("Index out of bounds");

		return field[i];
	}


	template<int dim>
	double DiscreteField<dim>::at(unsigned int i, unsigned int j) const
	{
		if(dim != 2)
			throw std::invalid_argument("Field dimension must equal 2 to use this accessor");
		if( (i >= discretization[0]) || (j	>= discretization[1]) )
			throw std::invalid_argument("Index out of bounds");

		return field[ i + j*discretization[0] ];
	}

	template<int dim>
	double DiscreteField<dim>::at(unsigned int i, unsigned int j, unsigned int k) const
	{
		if(dim != 3)
			throw std::invalid_argument("Field dimension must equal 3 to use this accessor");
		if( (i >= discretization[0]) || (j >= discretization[1])
			|| (k >= discretization[2]) )
			throw std::invalid_argument("Index out of bounds");

		return field[ i + j*discretization[0] + k*discretization[0]*discretization[1] ];
	}

	// with boundary condition checking
	template<int dim>
	double DiscreteField<dim>::at_bc(int i) const
	{
		if(dim != 1)
			throw std::invalid_argument("Field dimension must equal 1 to use this accessor");
		// if(i >= field.size())
		// 	throw std::invalid_argument("Index out of bounds");

		if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::WRAP))
			i = discretization[0]-1;
		else if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
			i = 0;
		else if( (i == discretization[0]) && (boundary_conditions[0] == BoundaryCondition::WRAP))
			i = 0;
		else if( (i == discretization[0]) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
			i = discretization[0]-1;

		return field[i];
	}


	template<int dim>
	double DiscreteField<dim>::at_bc(int i, int j) const
	{
		if(dim != 2)
			throw std::invalid_argument("Field dimension must equal 2 to use this accessor");
		// if( (i >= discretization[0]) || (j	>= discretization[1]) )
		// 	throw std::invalid_argument("Index out of bounds");

		if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::WRAP))
			i = discretization[0]-1;
		else if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
			i = 0;
		else if( (i == discretization[0]) && (boundary_conditions[0] == BoundaryCondition::WRAP))
			i = 0;
		else if( (i == discretization[0]) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
			i = discretization[0]-1;

		if( (j == -1) && (boundary_conditions[1] == BoundaryCondition::WRAP))
			j = discretization[1]-1;
		else if( (j == -1) && (boundary_conditions[1] == BoundaryCondition::REFLECT))
			j = 0;
		else if( (j == discretization[1]) && (boundary_conditions[1] == BoundaryCondition::WRAP))
			j = 0;
		else if( (j == discretization[1]) && (boundary_conditions[1] == BoundaryCondition::REFLECT))
			j = discretization[1]-1;

		return field[ i + j*discretization[0] ];
	}


	template<int dim>
	double DiscreteField<dim>::at_bc(int i, int j, int k) const
	{
		if(dim != 3)
			throw std::invalid_argument("Field dimension must equal 3 to use this accessor");
		// if( (i >= discretization[0]) || (j >= discretization[1])
		// 	|| (k >= discretization[2]) )
		// 	throw std::invalid_argument("Index out of bounds");


		if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::WRAP))
			i = discretization[0]-1;
		else if( (i == -1) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
			i = 0;
		else if( (i == discretization[0]) && (boundary_conditions[0] == BoundaryCondition::WRAP))
			i = 0;
		else if( (i == discretization[0]) && (boundary_conditions[0] == BoundaryCondition::REFLECT))
			i = discretization[0]-1;

		if( (j == -1) && (boundary_conditions[1] == BoundaryCondition::WRAP))
			j = discretization[1]-1;
		else if( (j == -1) && (boundary_conditions[1] == BoundaryCondition::REFLECT))
			j = 0;
		else if( (j == discretization[1]) && (boundary_conditions[1] == BoundaryCondition::WRAP))
			j = 0;
		else if( (j == discretization[1]) && (boundary_conditions[1] == BoundaryCondition::REFLECT))
			j = discretization[1]-1;

		if( (k == -1) && (boundary_conditions[2] == BoundaryCondition::WRAP))
			k = discretization[2]-1;
		else if( (k == -1) && (boundary_conditions[2] == BoundaryCondition::REFLECT))
			k = 0;
		else if( (k == discretization[2]) && (boundary_conditions[2] == BoundaryCondition::WRAP))
			k = 0;
		else if( (k == discretization[2]) && (boundary_conditions[2] == BoundaryCondition::REFLECT))
			k = discretization[2]-1;


		return field[ i + j*discretization[0] + k*discretization[0]*discretization[1] ];
	}


	// read and write access
	template<int dim>
	double& DiscreteField<dim>::at(unsigned int i)
	{
		if(dim != 1)
			throw std::invalid_argument("Field dimension must equal 1 to use this accessor");
		if(i >= field.size())
			throw std::invalid_argument("Index out of bounds");

		return field[i];
	}


	template<int dim>
	double& DiscreteField<dim>::at(unsigned int i, unsigned int j)
	{
		if(dim != 2)
			throw std::invalid_argument("Field dimension must equal 2 to use this accessor");
		if( (i >= discretization[0]) || (j	>= discretization[1]) )
			throw std::invalid_argument("Index out of bounds");

		return field[ i + j*discretization[0] ];
	}

	template<int dim>
	double& DiscreteField<dim>::at(unsigned int i, unsigned int j, unsigned int k)
	{
		if(dim != 3)
			throw std::invalid_argument("Field dimension must equal 3 to use this accessor");
		if( (i >= discretization[0]) || (j >= discretization[1])
			|| (k >= discretization[2]) )
			throw std::invalid_argument("Index out of bounds");

		return field[ i + j*discretization[0] + k*discretization[0]*discretization[1] ];
	}



	// access by point:
	template<int dim>
	double DiscreteField<dim>::at(const Point<dim>& p) const
	{
		return field[ indexFromPoint(p) ];
	}

	template<int dim>
	double& DiscreteField<dim>::at(const Point<dim>& p)
	{
		return field[ indexFromPoint(p) ];
	}

	// mutators: 
	template<int dim>
	void DiscreteField<dim>::setBounds(const Point<dim>& lower, const Point<dim>& upper)
	{
		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
			if(lower[dim_itr] > upper[dim_itr])
				throw std::invalid_argument("Upper point must be greater than lower in constructing field bounds");
	
		bottom_left = lower;
		top_right = upper;
	}

	template<int dim>
	void DiscreteField<dim>::setDiscretization(const std::array<unsigned int, dim>& disc)
	{
		discretization = disc;
		reinitialize();
	}

	template<int dim>
	void DiscreteField<dim>::setBoundaryConditions(const std::array<BoundaryCondition, dim>& bcs)
	{
		boundary_conditions = bcs;
	}


	template<int dim>
	void DiscreteField<dim>::projectFunction(const DiscreteFunction<dim>& function)
	{
		if(dim ==2)
		{
			for(unsigned int j = 0; j < discretization[1]; j++)
				for(unsigned int i = 0; i < discretization[0]; i++)
					this->at(i,j) = function.value( getCellCenterPoint(i,j) );
		}
		else if(dim == 3)
		{
			for(unsigned int k = 0; k < discretization[2]; k++)
				for(unsigned int j = 0; j < discretization[1]; j++)
					for(unsigned int i = 0; i < discretization[0]; i++)
						this->at(i,j,k) = function.value( getCellCenterPoint(i,j,k) );
		}
		else
		{
			throw std::invalid_argument("function not implented for field dim != 2 or 3");
		}
	}

	template<int dim>
	Point<dim> DiscreteField<dim>::getCellCenterPoint(unsigned int i, unsigned int j) const
	{
		if(dim != 2)
			throw std::invalid_argument("Cannot get cell center point at(i,j) for dim != 2");

		const double x = bottom_left[0] + (i+ 0.5)*getCellWidth(0);
		const double y = bottom_left[1] + (j + 0.5)*getCellWidth(1); 

		return Point<dim>(x,y);

	}

	template<int dim>
	Point<dim> DiscreteField<dim>::getCellCenterPoint(unsigned int i, unsigned int j,
		unsigned int k) const
	{
		if(dim != 3)
			throw std::invalid_argument("Cannot get cell center point at(i,j,k) for dim != 3");

		const double x = bottom_left[0] + (i+ 0.5)*getCellWidth(0);
		const double y = bottom_left[1] + (j + 0.5)*getCellWidth(1); 
		const double z = bottom_left[2] + (k + 0.5)*getCellWidth(2);

		return Point<dim>(x,y,z);	
	}

	template<int dim>
	void DiscreteField<dim>::print(std::ostream& out) const
	{
		if(dim == 1)
		{
			for(unsigned int i = 0; i < field.size(); i++)
				out << field[i] << " ";
			out << std::endl;
		}
		else if(dim == 2)
		{
			for(int j = 0; j < discretization[1]; j++)
			{
				for(int i = 0; i < discretization[0]; i++)
				{
					out << at(i,j) << " ";
				}
				out << std::endl;
			}
		}
		else if(dim == 3)
		{
			for(int k = 0; k < discretization[2]; k++)
			{
				for(int j = 0; j < discretization[1]; j++)
				{
					for(int i = 0; i < discretization[0]; i++)
					{
						out << at(i,j,k) << " ";
					}
					out << std::endl;
				}
				out << std::endl;
			} // for
		} // if dim == 3
		else
		{
			throw std::runtime_error("Field print not currently implented for dimension higher than 3");
		}

	}


	template<int dim>
	void DiscreteField<dim>::printInfo(std::ostream& out) const
	{
		out << "DISCRETE FIELD INFO: \n" << std::endl
			<< "Bottom left: " << bottom_left << std::endl
			<< "Top right: " << top_right << std::endl;
			
		out << "\nDISCRETIZATION: " << std::endl;
		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
			out << discretization[dim_itr] << " ";
		out << std::endl;

		out << "\nBOUNDARY CONDITIONS: " << std::endl;
		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
			out << getBoundaryConditionString(boundary_conditions[dim_itr]) << " ";
		out << std::endl;		
	}

	template<int dim>
	unsigned int DiscreteField<dim>::indexFromPoint(const Point<dim>& p) const
	{
		// base cases:
		unsigned int i = 0;
		unsigned int discProd = 1;
		
		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
		{
			// temporary fix ***
			// double issue_compiler_warning;
			// @todo: think of a better fix, like a buffer or global tolerance or something

			double reset_value = p[dim_itr] - bottom_left[dim_itr]; // should be no less than 0

			if(reset_value < 0 && reset_value > -1e-8)
			{
				std::cout << "fixing wrong value of " << reset_value << std::endl;
				reset_value = 0;
			}

			// maybe same thing other way...

			i += discProd*floor( reset_value*getInverseCellWidth(dim_itr) );
			discProd *=  discretization[dim_itr];
		}

		if( i >= discProd )
		{
			std::ostringstream message;
			message << "Index " << i << " from point " << p 
				<< " is out of range. Max index is " << discProd-1;

			throw std::runtime_error(message.str());
		}
	

		return i;
	}

}



#endif

