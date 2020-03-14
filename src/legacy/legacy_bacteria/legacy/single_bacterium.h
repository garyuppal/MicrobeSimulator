#ifndef SINGLE_BACTERIUM_H
#define SINGLE_BACTERIUM_H


#include <deal.II/base/point.h>
using dealii::Point;

#include "../geometry/geometry.h"
#include "../advection/velocity_interface.h"
#include "../fitness/fitness_base.h"

#include <iostream>
#include <string>
#include <sstream>
#include <array>

namespace MicrobeSimulator{

  template<int dim, int numchem>
  class SingleBacterium{
  private:
    Point<dim> location;

    std::array<double, numchem> secretionRates;

  public:
    // constructors:
    SingleBacterium();
    SingleBacterium(const Point<dim>& p, 
        const std::array<double, numchem>& rates);
    SingleBacterium(const SingleBacterium& b);

    // assignment
    SingleBacterium<dim, numchem>& operator=(const SingleBacterium<dim, numchem>& rhs);

    // accessors:
    Point<dim> getLocation() const; 
    std::array<double, numchem> getSecretionRates() const;
    double getSecretionRate(unsigned int i) const;

    // mutators:
    void setLocation(const Point<dim>& p);
    void setSecretionRates(const std::array<double, numchem>& rates);
    void setSecretionRate(unsigned int i, double rate); 
 
    // functions:
    void randomStep(double timeStep, double diffusionConstant, const Geometry<dim>& geometry,
        const VelocityInterface<dim>& velocity); 
    // void randomStep(double timeStep, double diffusionConstant, 
    // 	const Geometry<dim>* const geometry, const AdvectionField<dim>* const advection = NULL );

    // void randomStep(double timeStep, double diffusionConstant, 
    //     const Geometry<dim>* const geometry, const DoFHandler<dim>& stokes_dof,
    //     const BlockVector<double> stokes_solution);
    // // *** probably better to just pass a function pointer for advection and geometry boundary conditions ***

    // void randomStep(double timeStep, double diffusionConstant, 
    //     const Geometry<dim>* const geometry, const DoFHandler<dim>& stokes_dof,
    //     const BlockVector<double> stokes_solution, const PointCellMap<dim>& pcm);

    double getFitness(const FitnessBase<dim,numchem>& fitness_function);

    void mutate(double deltaSecretion, double original_secretion_rate, bool binary_mutation); 
        // what about for multiple public goods ?

    void printBacterium(std::ostream& out) const; 

  }; // class Bacterium


// IMPLEMENTATION:
// ---------------------------------------------------------------------------------------------------
    template<int dim, int numchem>
    SingleBacterium<dim, numchem>::SingleBacterium()
        : 
        location()
    {}


    template<int dim, int numchem>
    SingleBacterium<dim, numchem>::SingleBacterium(
        const Point<dim>& p,  
        const std::array<double, numchem>& rates)
        : 
        location(p),
        secretionRates(rates)
    {}


    template<int dim, int numchem>
    SingleBacterium<dim, numchem>::SingleBacterium(
        const SingleBacterium& b)
    {
        location = b.location;
        secretionRates = b.secretionRates;
    }


    template<int dim, int numchem>
    SingleBacterium<dim, numchem>&
    SingleBacterium<dim, numchem>::operator=(
        const SingleBacterium<dim, numchem>& rhs)
    {
       // check for self copy:
        if(this == &rhs)
          return *this;

        // copy:
        location = rhs.location;
        secretionRates = rhs.secretionRates;

        return *this;
    }


    // accessors:
    template<int dim, int numchem>
    Point<dim> 
    SingleBacterium<dim, numchem>::getLocation() const
    {
        return location;
    }


    template<int dim, int numchem>
    std::array<double, numchem> 
    SingleBacterium<dim, numchem>::getSecretionRates() const
    {
        return secretionRates; 
    }

    
    template<int dim, int numchem>
    double 
    SingleBacterium<dim, numchem>::getSecretionRate(unsigned int i) const 
    {
        return secretionRates.at(i); 
    }

    // mutators:
    template<int dim, int numchem>
    void SingleBacterium<dim, numchem>::setLocation(const Point<dim>& p)
    {
        location = p;
    }


    template<int dim, int numchem>
    void 
    SingleBacterium<dim, numchem>::setSecretionRates(
        const std::array<double, numchem>& rates)
    {
        secretionRates = rates;
    }

    template<int dim, int numchem>
    void 
    SingleBacterium<dim, numchem>::setSecretionRate(unsigned int i,
        double rate)
    {
        secretionRates[i] = rate;
    }


    // functions:
    template<int dim, int numchem>
    void 
    SingleBacterium<dim, numchem>::randomStep(double timeStep, double diffusionConstant, 
        const Geometry<dim>& geometry,
        const VelocityInterface<dim>& velocity)
    {
        Point<dim> old_location = location;

        const double theta = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));
        const double phi = dealii::numbers::PI*((double) rand() / (RAND_MAX));

        Point<dim> randomPoint = (dim == 2) ? Point<dim>(std::cos(theta),std::sin(theta))
                                    : Point<dim>(std::cos(phi), std::sin(phi)*std::cos(theta),
                                        std::sin(phi)*std::sin(theta));

        location += std::sqrt(2*dim*timeStep*diffusionConstant)*randomPoint
                + timeStep*velocity.value(location);

        geometry.checkBoundaries(old_location, location); 
    }
// legacy:
/*
    template<int dim, int numchem>
    void 
    SingleBacterium<dim, numchem>::randomStep(double timeStep, double diffusionConstant, 
        const Geometry<dim>* const geometry, const AdvectionField<dim>* const advection)
    {
        Point<dim> old_location = location;
        Point<dim> randomPoint;

        if(dim == 2)
        {
            const double theta = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));
            randomPoint(0) = std::cos(theta);
            randomPoint(1) = std::sin(theta);
        }
        else if(dim == 3)
        {
            const double theta = dealii::numbers::PI*((double) rand() / (RAND_MAX));
            const double phi = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));

            randomPoint(0) = std::cos(theta);
            randomPoint(1) = std::sin(theta)*std::cos(phi);
            randomPoint(2) = std::sin(theta)*std::sin(phi);
        }
        else
        {
            throw std::runtime_error("Random step not implemented for dim != 2,3");
        }

        location += std::sqrt(2*dim*timeStep*diffusionConstant)*randomPoint;

        if(advection != NULL)
            location += advection->value(old_location)*timeStep;

            // std::cout << " need to implement advection" << std::endl;
            // location += velocity(location)*dt

        geometry->checkBoundaries(old_location, location);
    }


    template<int dim, int numchem>
    void 
    SingleBacterium<dim, numchem>::randomStep(double timeStep, double diffusionConstant, 
        const Geometry<dim>* const geometry, const DoFHandler<dim>& stokes_dof,
        const BlockVector<double> stokes_solution)
    {
        Point<dim> old_location = location;
        Point<dim> randomPoint;

        Point<dim> advection_value; 

        if(dim == 2)
        {
            const double theta = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));
            randomPoint(0) = std::cos(theta);
            randomPoint(1) = std::sin(theta);
        }
        else if(dim == 3)
        {
            const double theta = dealii::numbers::PI*((double) rand() / (RAND_MAX));
            const double phi = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));

            randomPoint(0) = std::cos(theta);
            randomPoint(1) = std::sin(theta)*std::cos(phi);
            randomPoint(2) = std::sin(theta)*std::sin(phi);
        }
        else
        {
            throw std::runtime_error("Random step not implemented for dim != 2,3");
        }

        location += std::sqrt(2*dim*timeStep*diffusionConstant)*randomPoint;

        // compute velocity:
        Vector<double> velocity;
        velocity.reinit(dim + 1); // +1 for pressure

        dealii::VectorTools::point_value(stokes_dof,
                                stokes_solution,
                                old_location,
                                velocity);

        for(unsigned int dim_itr = 0; dim_itr<dim; dim_itr++)
            advection_value[dim_itr] = velocity[dim_itr];

        location += advection_value*timeStep;

        geometry->checkBoundaries(old_location, location);
    }


    template<int dim, int numchem>
    void 
    SingleBacterium<dim, numchem>::randomStep(double timeStep, double diffusionConstant, 
        const Geometry<dim>* const geometry, const DoFHandler<dim>& stokes_dof,
        const BlockVector<double> stokes_solution, const PointCellMap<dim>& pcm)
    {
        Point<dim> old_location = location;
        Point<dim> randomPoint;

        Point<dim> advection_value; 

        if(dim == 2)
        {
            const double theta = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));
            randomPoint(0) = std::cos(theta);
            randomPoint(1) = std::sin(theta);
        }
        else if(dim == 3)
        {
            const double theta = dealii::numbers::PI*((double) rand() / (RAND_MAX));
            const double phi = 2*dealii::numbers::PI*((double) rand() / (RAND_MAX));

            randomPoint(0) = std::cos(theta);
            randomPoint(1) = std::sin(theta)*std::cos(phi);
            randomPoint(2) = std::sin(theta)*std::sin(phi);
        }
        else
        {
            throw std::runtime_error("Random step not implemented for dim != 2,3");
        }

        location += std::sqrt(2*dim*timeStep*diffusionConstant)*randomPoint;

        // compute velocity:
        Vector<double> velocity;
        velocity.reinit(dim + 1); // +1 for pressure

        vector_point_value(pcm, *geometry, stokes_dof, stokes_solution, 
            old_location, velocity);
        // dealii::VectorTools::point_value(stokes_dof,
        //                         stokes_solution,
        //                         old_location,
        //                         velocity);

        for(unsigned int dim_itr = 0; dim_itr<dim; dim_itr++)
            advection_value[dim_itr] = velocity[dim_itr];

        location += advection_value*timeStep;

        geometry->checkBoundaries(old_location, location);
    }
*/


    template<int dim, int numchem>
    void 
    SingleBacterium<dim, numchem>::mutate(double deltaSecretion, 
        double original_secretion_rate, bool binary_mutation)
    {
        // for public goods -- up to numchem - 1
        // pick one to mutate

        int chem_id = rand() % (numchem-1); 

        if(binary_mutation == true)
        {
            secretionRates[chem_id] = (secretionRates[chem_id] == 0 ?
                                        original_secretion_rate :
                                        0);
        }
        else
        {
            const double rate = secretionRates[chem_id] 
                + deltaSecretion*( 2.*((double) rand() / (RAND_MAX)) - 1.);

            secretionRates[chem_id] = ( rate < 0 ?
                                        0 :
                                        rate );  // absorb at 0
        }
    }

    template<int dim, int numchem>
    double 
    SingleBacterium<dim, numchem>::getFitness(
        const FitnessBase<dim,numchem>& fitness_function)
    {
        return fitness_function.value(location, secretionRates);
    }


    template<int dim, int numchem>
    void 
    SingleBacterium<dim, numchem>::printBacterium(std::ostream& out) const
    {
        out << location;
        for(unsigned int i = 0; i < numchem; i++)
            out << " " << secretionRates[i];
        out << std::endl;
    }

}

#endif  // bacterium.h
