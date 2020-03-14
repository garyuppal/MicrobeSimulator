#ifndef MICROBE_SIMULATOR_CELL_ITERATOR_MAP_H
#define MICROBE_SIMULATOR_CELL_ITERATOR_MAP_H

// #include <deal.II/base/quadrature_lib.h>
// #include <deal.II/base/logstream.h>
// #include <deal.II/base/function.h>
// #include <deal.II/base/utilities.h>

// #include <deal.II/grid/tria.h>
// #include <deal.II/grid/grid_generator.h>
// #include <deal.II/grid/tria_accessor.h>
// #include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
// #include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
// #include <deal.II/dofs/dof_renumbering.h>
// #include <deal.II/dofs/dof_accessor.h>
// #include <deal.II/dofs/dof_tools.h>

// #include <deal.II/fe/fe_q.h>
// #include <deal.II/fe/fe_system.h>
// #include <deal.II/fe/fe_values.h>

// #include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <fstream>
#include <memory>
#include <array>
#include <fstream>
#include <string>
#include <sstream>

// #include "../utility/argparser.h"
#include "../geometry/geometry.h"
// #include "../advection/advection.h"
// #include "../bacteria/bacteria.h"


namespace MicrobeSimulator{ namespace FETools{
	// using namespace dealii;

	template<int dim>
	class PointCellMap{
	public:
		PointCellMap();

		void initialize(const Geometry<dim>& geometry, 
						const dealii::DoFHandler<dim>& dofh,
						double resolution = 0.2); // done

		void initialize(const Geometry<dim>& geometry, 
						const dealii::DoFHandler<dim>& dofh,
						std::array<double, dim> resolutions); // done

		void initialize(const Geometry<dim>& geometry,
						const dealii::DoFHandler<dim>& dofh,
						std::array<unsigned int, dim> disc); // done

		std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, 
    		dealii::Point<dim> >  get_cell_point_pair(unsigned int i) const; // done

		std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, 
    		dealii::Point<dim> > get_cell_point_pair(const Point<dim>& p) const; // done


		// @todo: these are in common with discrete_field -- should probably implement
		// as a base class from which this and discrete field derive
    		// or a templated class...
		unsigned int getTotalSize() const; // done
		unsigned int getMapSize() const; //done // for checking... should be equal
		double getCellWidth(unsigned int dimension) const; // done
		double getInverseCellWidth(unsigned int dimension) const; // done
		Point<dim> getCellCenterPoint(unsigned int i, unsigned int j) const; // done
		Point<dim> getCellCenterPoint(unsigned int i, unsigned int j,unsigned int k) const; // done
		void printInfo(std::ostream& out) const;  // done


	private:
	    std::vector<std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, 
	    		dealii::Point<dim> > >  	iterator_map;

		const Geometry<dim>* 				geometry;
		std::array<unsigned int, dim> 		discretization;

		unsigned int indexFromPoint(const Point<dim>& p) const; // done
		void setupBase(const std::array<unsigned int, dim> disc); // done
	};



// IMPLEMENTATION:
//-------------------------------------------------------------------------------------------
	template<int dim>
	PointCellMap<dim>::PointCellMap()
	{}


	template<int dim>
	void 
	PointCellMap<dim>::setupBase(const std::array<unsigned int, dim> disc)
	{
		discretization = disc;
		unsigned int size = getTotalSize();
		iterator_map.reserve(size);
	}


	template<int dim>
	void
	PointCellMap<dim>::initialize(const Geometry<dim>& geo, 
						const dealii::DoFHandler<dim>& dofh,
						double resolution)
	{
		// use approx resolution to get discretization:
		std::array<unsigned int, dim> disc;

		for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
			disc[dim_itr] = std::ceil(geo.getWidth(dim_itr)/resolution);

		initialize(geo, dofh, disc);
	}


	template<int dim>
	void 
	PointCellMap<dim>::initialize(const Geometry<dim>& geo, 
					const dealii::DoFHandler<dim>& dofh,
					std::array<double, dim> resolutions)
	{
		// use array of resolutions to get different discretizations:
		std::array<unsigned int, dim> disc;

		for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
			disc[dim_itr] = std::ceil( geo.getWidth(dim_itr) / resolutions[dim_itr] );

		initialize(geo, dofh, disc);
	}

	template<int dim>
	void 
	PointCellMap<dim>::initialize(const Geometry<dim>& geo,
									const dealii::DoFHandler<dim>& dofh,
									std::array<unsigned int, dim> disc)
	{
		geometry = &geo;
		setupBase(disc);

		if(dim == 2)
		{
			for(unsigned int j = 0; j < discretization[1]; j++)
				for(unsigned int i = 0; i < discretization[0]; i++)
				{
					Point<dim> temp_point = getCellCenterPoint(i,j);

					if(geometry->isInDomain(temp_point))
					{
						std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, 
							dealii::Point<dim> >
						cell_point =
						dealii::GridTools::find_active_cell_around_point (
									dealii::StaticMappingQ1<dim>::mapping, 
						            dofh, 
						            temp_point);

						iterator_map.push_back(cell_point);
					}
					else // could use unordered map instead because of void spaces ***
					{
						// int switch_to_map; // @todo -- instead of vector, use unsorted_map

						std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, 
							dealii::Point<dim> >
						cell_point(dofh.end(), Point<dim>() );

						iterator_map.push_back(cell_point);
					}
		      	}
      	} // if dim == 2
      	else if(dim == 3)
      	{
      		for(unsigned int k = 0; k < discretization[2]; k++)
				for(unsigned int j = 0; j < discretization[1]; j++)
					for(unsigned int i = 0; i < discretization[0]; i++)
					{
						Point<dim> temp_point = getCellCenterPoint(i,j,k);

						if(geometry->isInDomain(temp_point))
						{		
							std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, 
								dealii::Point<dim> >
							cell_point =
							dealii::GridTools::find_active_cell_around_point (
										dealii::StaticMappingQ1<dim>::mapping, 
							            dofh, 
							            temp_point);

							iterator_map.push_back(cell_point);
						}
						else
						{
							std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, 
								dealii::Point<dim> >
							cell_point(dofh.end(), Point<dim>() );

							iterator_map.push_back(cell_point);
						}
			      	}
      	}
      	else
      	{
      		throw std::invalid_argument("cell map only for 2 or 3 dim");
      	}

	} // iniitalize from geometry


	template<int dim>
	unsigned int 
	PointCellMap<dim>::getTotalSize() const
	{
		unsigned int size = 1;
		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
			size *= discretization[dim_itr];

		return size;
	}


	template<int dim>
	unsigned int
	PointCellMap<dim>::getMapSize() const
	{
		return iterator_map.size();
	}


	template<int dim>
	std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, 
	    		dealii::Point<dim> >
	PointCellMap<dim>::get_cell_point_pair(unsigned int i) const
	{
		return iterator_map[i];
	}


	template<int dim>
	std::pair<typename dealii::DoFHandler<dim>::active_cell_iterator, 
	    		dealii::Point<dim> >
	PointCellMap<dim>::get_cell_point_pair(const Point<dim>& p) const
	{
		return iterator_map[indexFromPoint(p)];
	}


	template<int dim>
	double 
	PointCellMap<dim>::getCellWidth(unsigned int dimension) const
	{
		if(dimension >= dim)
			throw std::invalid_argument("Cannot get cell width for non-existent dimension");

		return (geometry->getTopRightPoint()[dimension] 
			- geometry->getBottomLeftPoint()[dimension])/discretization[dimension];
	}


	template<int dim>
	double 
	PointCellMap<dim>::getInverseCellWidth(unsigned int dimension) const
	{
		return 1.0/getCellWidth(dimension);
	}


	template<int dim>
	Point<dim> 
	PointCellMap<dim>::getCellCenterPoint(unsigned int i, unsigned int j) const
	{
		if(dim != 2)
			throw std::invalid_argument("Cannot get cell center point at(i,j) for dim != 2");

		const double x = geometry->getBottomLeftPoint()[0] + (i+ 0.5)*getCellWidth(0);
		const double y = geometry->getBottomLeftPoint()[1] + (j + 0.5)*getCellWidth(1); 

		return Point<dim>(x,y);
	}

	template<int dim>
	Point<dim> 
	PointCellMap<dim>::getCellCenterPoint(unsigned int i, unsigned int j,
		unsigned int k) const
	{
		if(dim != 3)
			throw std::invalid_argument("Cannot get cell center point at(i,j,k) for dim != 3");

		const double x = geometry->getBottomLeftPoint()[0] + (i+ 0.5)*getCellWidth(0);
		const double y = geometry->getBottomLeftPoint()[1] + (j + 0.5)*getCellWidth(1); 
		const double z = geometry->getBottomLeftPoint()[2] + (k + 0.5)*getCellWidth(2);

		return Point<dim>(x,y,z);	
	}


	template<int dim>
	unsigned int 
	PointCellMap<dim>::indexFromPoint(const Point<dim>& p) const
	{
		// base cases:
		unsigned int i = 0;
		unsigned int discProd = 1;

		// calculate buffer as max resolution:

		double min_buffer = 1e10;
		double max_buffer = 0;
		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
		{
			const double res_along_dim = getCellWidth(dim_itr);

			max_buffer = std::max(max_buffer,res_along_dim);
			min_buffer = std::min(min_buffer, res_along_dim);
		}
		const double buffer = 0.5*(min_buffer + max_buffer);

		// std::cout << "buffer is: " << buffer << std::endl;
		
		Point<dim> buffered_point; 

		geometry->addPointBuffer(buffer,p,buffered_point); 

		if(!geometry->isInDomain(buffered_point))
		{
			std::ostringstream message;
			message << "buffered_point: " << buffered_point << " out of domain" 
				<< " original point at: " << p <<  std::endl;
			throw std::runtime_error(message.str());
		}
			// shift point relative to bottom left and add buffer

		// std::cout << "\nbuffer = " << buffer << std::endl
		// 	<< "intial point: " << p << std::endl
		// 	<< "\tbuffered point: " << buffered_point << std::endl;

		for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
		{
			const double buffered_value = buffered_point[dim_itr] - geometry->getBottomLeftPoint()[dim_itr];

			i += discProd*floor( buffered_value*getInverseCellWidth(dim_itr) );
			discProd *=  discretization[dim_itr];
		}

		if( i >= discProd )
		{
			printInfo(std::cout); 

			for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
			{
				const double buffered_value = buffered_point[dim_itr] - geometry->getBottomLeftPoint()[dim_itr];

				unsigned int j = discProd*floor( buffered_value*getInverseCellWidth(dim_itr) );
				std::cout << "dim " << dim_itr << "disc: " << j << std::endl;
			}

			std::ostringstream message;
			message << "Index " << i << " from point " << p 
				<< " with buffered point: " << buffered_point 
				<< " is out of range. Max index is " << discProd-1;

			throw std::runtime_error(message.str());
		}
	

		return i;
	}


	template<int dim>
	void 
	PointCellMap<dim>::printInfo(std::ostream& out) const
	{
	    out << "\n\n-----------------------------------------------------" << std::endl
	    	 << "\t\t ITERATOR MAP INFO:"
	    	 << "\n-----------------------------------------------------" << std::endl;

	   	out << "MAP SIZE: " << iterator_map.size() << std::endl << std::endl;

		out << "DIMENSION: " << dim << std::endl
			<< "\t BottomLeft: " << geometry->getBottomLeftPoint() << std::endl
			<< "\t TopRight: " << geometry->getTopRightPoint() << std::endl << std::endl;

		out << "DISCRETIZATION: " << std::endl;
			for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
				out << discretization[dim_itr] << std::endl;

		out << "\n-----------------------------------------------------\n\n" << std::endl;
	}

}} // close namespaces
#endif