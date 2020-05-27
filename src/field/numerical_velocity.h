#ifndef MICROBESIMULATOR_NUMERICAL_VELOCITY_H
#define MICROBESIMULATOR_NUMERICAL_VELOCITY_H

#include <deal.II/base/point.h>
using dealii::Point;
using dealii::Tensor;

#include <array>

#include "./discrete_field.h"

namespace MicrobeSimulator{

	template<int dim>
	class NumericalVelocity{
	public:
		NumericalVelocity();
		NumericalVelocity(std::ifstream &vx_file, std::ifstream& vy_file);
		NumericalVelocity(std::ifstream &vx_file, std::ifstream& vy_file,
			std::ifstream& vz_file);
		NumericalVelocity(const NumericalVelocity& nv);

		void initialize(std::ifstream &vx_file, std::ifstream& vy_file,
			const Point<dim>& bottom_left, 
			double velocity_scale, double geometry_scale);
		void initialize(std::ifstream &vx_file, std::ifstream& vy_file,
			std::ifstream& vz_file, 
			const Point<dim>& bottom_left, 
			double velocity_scale, double geometry_scale);

		Tensor<1,dim> value (const Point<dim> &p) const;

		void printInfo(std::ostream& out) const;
	private:
		std::array<DiscreteField<dim>, dim> velocity; 

		void readInFile(unsigned int dimension, std::ifstream infile);
	};


//IMPLEMENTATIONS
//-------------------------------------------------------------------------------------
	template<int dim>
	NumericalVelocity<dim>::NumericalVelocity()
	{}

	template<int dim>
	NumericalVelocity<dim>::NumericalVelocity(std::ifstream &vx_file, std::ifstream& vy_file)
	{
		initialize(vx_file,vy_file);
	}

	template<int dim>
	NumericalVelocity<dim>::NumericalVelocity(std::ifstream &vx_file, std::ifstream& vy_file,
		std::ifstream& vz_file)
	{
		initialize(vx_file,vy_file,vz_file);
	}

	template<int dim>	
	NumericalVelocity<dim>::NumericalVelocity(const NumericalVelocity& nv)
	{
		velocity = nv.velocity;
	}

	template<int dim>
	Tensor<1,dim> NumericalVelocity<dim>::value (const Point<dim> &p) const
	{
		Tensor<1,dim> result;

		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
			result[dim_itr] = velocity[dim_itr].at(p);

		return result;
	}


	template<int dim>
	void NumericalVelocity<dim>::initialize(std::ifstream &vx_file, std::ifstream& vy_file,
		const Point<dim>& bottom_left, double velocity_scale, double geometry_scale)
	{
		if(dim != 2)
			throw std::runtime_error("Cannot initialize velocity, dimension must be 2 to use this initializer");

		velocity[0].initialize(vx_file,bottom_left,velocity_scale, geometry_scale);
		velocity[1].initialize(vy_file,bottom_left,velocity_scale, geometry_scale);
	}

	template<int dim>
	void NumericalVelocity<dim>::initialize(std::ifstream &vx_file, std::ifstream& vy_file,
		 std::ifstream& vz_file, const Point<dim>& bottom_left, 
		 double velocity_scale,  double geometry_scale)
	{
		if(dim != 3)
			throw std::runtime_error("Cannot initialize velocity, dimension must be 3 to use this initializer");
			
		velocity[0].initialize(vx_file,bottom_left,velocity_scale, geometry_scale);
		velocity[1].initialize(vy_file,bottom_left,velocity_scale, geometry_scale);
		velocity[2].initialize(vz_file,bottom_left,velocity_scale, geometry_scale);
	}

	template<int dim>
	void NumericalVelocity<dim>::readInFile(unsigned int dimension, std::ifstream infile)
	{}


	template<int dim>
	void NumericalVelocity<dim>::printInfo(std::ostream& out) const
	{
		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
		{
			out << "\n..............................................\n";
			out << "NUMERICAL VELOCITY: component " << dim_itr << std::endl;
			velocity[dim_itr].printInfo(out); 
		}
	}

}

#endif

