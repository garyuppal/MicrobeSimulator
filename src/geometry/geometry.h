#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cstddef>
#include <array>

#include "sphere.h"
#include "hyper_rectangle.h"
#include "line.h"

#include "../utility/parameter_handler.h"


namespace MicrobeSimulator{

	/** \brief Enum type for boundary conditions */
	enum class BoundaryCondition : int
	{
		WRAP = 0, REFLECT = 1, OPEN = 2 
	}; 

	/** \brief Convert boundary condition type to string */
	std::string getBoundaryConditionString(BoundaryCondition bc)
	{
		std::string result = "";
		if(bc == BoundaryCondition::WRAP)
			result = "Periodic";
		else if(bc == BoundaryCondition::REFLECT)
			result = "Neumann";
		else if(bc == BoundaryCondition::OPEN)
			result = "Open";
		else
			result = "ERROR";

		return result;
	}

	/** \brief Convert a valid string to enum boundary condition type */
	BoundaryCondition stringToBoundaryCondition(std::string value)
	{
		if( boost::iequals(value, "WRAP") )
			return BoundaryCondition::WRAP;
		else if( boost::iequals(value, "REFLECT") )
			return BoundaryCondition::REFLECT;
		else if( boost::iequals(value, "OPEN") )
			return BoundaryCondition::OPEN;
		else
		{
			std::string msg = "Invalid string! Could not convert <"
							+ value
							+ "> to boundary condition type.";
			throw std::invalid_argument(msg.c_str());
		}
	}

// GEOMETRY CLASS
// ------------------------------------------------------------------------------------
/** \brief Geometry class for mesh constuction support and
* handling bacteria collisions with boundaries
*/
template<int dim>
class Geometry{
public:
	Geometry();

	// accessors:
	Point<dim> getBottomLeftPoint() const;
	Point<dim> getTopRightPoint() const;

	double getWidth(unsigned int direction) const;
	std::array<BoundaryCondition, dim> getBoundaryConditions() const;

	// sphere accessors:
	std::vector<Sphere<dim> > getSpheres() const;
	unsigned int getNumberSpheres() const;
	Sphere<dim> getSphereAt(unsigned int i) const;

	// rectangle accessors:
	std::vector<HyperRectangle<dim> > getRectangles() const;
	unsigned int getNumberRectangles() const;
	HyperRectangle<dim> getRectangleAt(unsigned int i) const;

	// line accessors:
	std::vector<Line> getLines() const;
	unsigned int getNumberLines() const;

	// mutators: (used by builder to constuct geometry)
	void setBottomLeftPoint(const Point<dim>& lower); 
	void setTopRightPoint(const Point<dim>& upper); 

	void setBoundaryConditions(const std::array<BoundaryCondition, dim>& bcs);
	void setBoundaryCondition(unsigned int i, BoundaryCondition bc);
	
	void addSphere(const Sphere<dim>& sp);
	void addRectangle(const HyperRectangle<dim>& rect);
	void addLine(const Line& line);

	// HANDLING BOUNDARIES:
	void checkBoundaries(const Point<dim>& oldPoint, 
	                   Point<dim>& newPoint,
	                   const double buffer = 0.005) const; 

	bool isInDomain(const Point<dim>& location, double buffer=0) const; 

	void addPointBuffer(const double buffer,
	                   const Point<dim>& test_point,
	                   Point<dim>& buffered_point) const; 

	void printInfo(std::ostream& out) const;
	void outputGeometry(std::string output_directory = ".") const;

	std::vector<Point<dim> > getQuerryPoints(double resolution = 0.2) const; 

private:

	// bounding box:
	Point<dim> bottom_left;
	Point<dim> top_right;

	std::array<BoundaryCondition, dim> boundary_conditions; 

	// interior obstacles:
	std::vector<Sphere<dim> > spheres;
	std::vector<HyperRectangle<dim> > rectangles;

	// bounding lines:
	std::vector<Line> lines;
}; // class Geometry{}



// IMPLEMENTATION
// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------

/** \brief Constuctor */
template<int dim>
Geometry<dim>::Geometry()
{}


// ACCESSORS:
// ----------------------------------------------
/** \brief Return bottom left corner of bounding domain */
template<int dim>
Point<dim> 
Geometry<dim>::getBottomLeftPoint() const
{
	return bottom_left;
}

/** \brief Return top right corner of bounding domain */
template<int dim>
Point<dim> 
Geometry<dim>::getTopRightPoint() const
{
	return top_right;
}

/** \brief Return width of bounding geometry along specified dimension. */
template<int dim>
double 
Geometry<dim>::getWidth(unsigned int direction) const
{
	if(direction >= dim)
		throw std::invalid_argument("Desired dimension to get width does not exist");
	return top_right[direction] - bottom_left[direction];
}

/** \brief Return array of boundary conditions */
template<int dim>
std::array<BoundaryCondition, dim> 
Geometry<dim>::getBoundaryConditions() const
{
	return boundary_conditions;
}

// sphere accessors:

/** \brief Return vector copy of spheres in geometry */
template<int dim>
std::vector<Sphere<dim> > 
Geometry<dim>::getSpheres() const
{
	return spheres;
}

/** \brief Returns number of spheres in geometry */
template<int dim>
unsigned int 
Geometry<dim>::getNumberSpheres() const
{
	return spheres.size();
}

/** \brief Returns ith sphere in geometry */
/** Spheres are stored in no particular order. Typically access to spheres requires
* access to all spheres, so this method would usually be used in a loop over
* all spheres in the domain.
*/
template<int dim>
Sphere<dim> 
Geometry<dim>::getSphereAt(unsigned int i) const
{
	return spheres[i];
}


// rectangle accessors:

/** \brief Returns vector of all interior rectangles in geometry */
template<int dim>    
std::vector<HyperRectangle<dim> > 
Geometry<dim>::getRectangles() const
{
	return rectangles;
}

/** \brief Return number of rectangles in geometry.*/
template<int dim>  
unsigned int 
Geometry<dim>::getNumberRectangles() const
{
	return rectangles.size();
}

/** \brief Return copy of ith rectangle in geometry */
/** Rectangles are not stored in any particular order. 
*/
template<int dim>
HyperRectangle<dim> 
Geometry<dim>::getRectangleAt(unsigned int i) const
{
	return rectangles[i];
}

// line accessors:
/** \brief Return vector of lines in geometry */
template<int dim>
std::vector<Line> 
Geometry<dim>::getLines() const
{
	return lines;
}

/** \brief Return number of boundary lines in geometry */
template<int dim>
unsigned int 
Geometry<dim>::getNumberLines() const
{
	return lines.size();
}

// MUTATORS:
// ----------------------------------------------

/** \brief Set the bottom left corner of geometry */
template<int dim>
void 
Geometry<dim>::setBottomLeftPoint(const Point<dim>& lower)
{
	bottom_left = lower;
}

/** \brief Set the top right point of geometry */
template<int dim>
void 
Geometry<dim>::setTopRightPoint(const Point<dim>& upper)
{
	top_right = upper;
}

/** \brief Set boundary conditions to specified values */
template<int dim>
void 
Geometry<dim>::setBoundaryConditions(const std::array<BoundaryCondition, dim>& bcs)
{
	boundary_conditions = bcs;
}

/** \brief Set ith boundary condition to specified value */
template<int dim>
void
Geometry<dim>::setBoundaryCondition(unsigned int i, BoundaryCondition bc)
{
	boundary_conditions[i] = bc;
}

/** \brief Add sphere to geometry */
template<int dim>
void 
Geometry<dim>::addSphere(const Sphere<dim>& sp)
{
	spheres.emplace_back(sp);
}

/** \brief Add interior rectangle to geometry */
template<int dim>
void 
Geometry<dim>::addRectangle(const HyperRectangle<dim>& rect)
{
	rectangles.emplace_back(rect);
}

/** \brief Add line to geometry */
template<int dim>
void
Geometry<dim>::addLine(const Line& line)
{
	lines.emplace_back(line);
}


// FUNCTIONS:
// ------------------------------------------------------------------------

/** \brief Check if moving from old point to new point is valid and reflect off 
* obstacles as needed.
*/
/** Optional 3rd buffer parameter gives an extra ``cushion'' to boundaries. 
* Points within buffer distance to reflecting type boundaries (including interior boundaries)
* will be reflected. Wrapping boundaries are treated by moving the point to the 
* appropriate opposite end of the domain. Open boundaries are left alone. 
* Typically the bacteria class will handle open boundaries by removing microbes that
* move outside of the domain. 
*/
template<int dim>
void 
Geometry<dim>::checkBoundaries(const Point<dim>& oldPoint, Point<dim>& newPoint,
								const double buffer) const
{
	const double tolerance = 1e-8;

	// check interior and exterior spheres:
	unsigned int number_spheres = spheres.size();
	for(unsigned int sphere_id = 0; sphere_id < number_spheres; ++sphere_id)
		if( (spheres[sphere_id].distance_from_border(newPoint) - buffer) < tolerance)
		{
			spheres[sphere_id].reflectPoint(oldPoint, newPoint, buffer); 
			break;
		}

	// check interior rectangle,  buffer optional
	unsigned int number_rectangles = rectangles.size();
	for(unsigned int rect_id = 0; rect_id < number_rectangles; ++rect_id)
		if( rectangles[rect_id].distance_from_border(newPoint, buffer) < tolerance ) 
		{		
			rectangles[rect_id].reflectPoint(oldPoint, newPoint, buffer); 
			break;        
		}

	// check interior lines:
	unsigned int n_lines = lines.size();
	for(unsigned int i = 0; i < n_lines; ++i)
		if( !lines[i].is_in_bounds(newPoint, buffer) )
		{
			lines[i].reflectPoint(oldPoint, newPoint, buffer);
			break;
		}

	// check bounding box:
	for(unsigned int i = 0; i < dim; i++)
	{
		// bottom_left gives lower boundaries, top_right gives upper
		if( newPoint[i] < (bottom_left[i]+buffer) )
		{
			if( ( newPoint[i] < bottom_left[i] ) && 
				(boundary_conditions[i] == BoundaryCondition::WRAP) )
				 newPoint[i] = newPoint[i] + (top_right[i] - bottom_left[i]);
			else //if(boundary_conditions[i] == BoundaryCondition::REFLECT)
				 newPoint[i] = 2.0*(bottom_left[i]+buffer) - newPoint[i]; 
			 	// -- use reflect for open as well on left
		}
		else if( newPoint[i] > (top_right[i]-buffer) )
		{
			if( (newPoint[i] > top_right[i]) &&
				(boundary_conditions[i] == BoundaryCondition::WRAP) )
				 newPoint[i] = newPoint[i] - (top_right[i] - bottom_left[i]);
			else if(boundary_conditions[i] == BoundaryCondition::REFLECT) 
				 newPoint[i]= 2.0*(top_right[i]-buffer) - newPoint[i];
			// else --- is open
		}
	} // for dim
} // check_boundaries()


/** \brief Checks if point is within the domain and not inside any
* interior obstacles.
*/
template<int dim>
bool 
Geometry<dim>::isInDomain(const Point<dim>& location, double buffer) const
{
	// check bounding box:
	for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
		if( location[dim_itr] < (bottom_left[dim_itr] + buffer)
			|| location[dim_itr] > (top_right[dim_itr] - buffer) )
	   		return false;

	// check interior/exterior spheres:
	unsigned int number_spheres = spheres.size();
	for(unsigned int sphere_id = 0; sphere_id < number_spheres; ++sphere_id)
		if( spheres[sphere_id].isInSphere(location, buffer))
	   		return false;

	// check interior rectangles:
	unsigned int number_rectangles = rectangles.size();
	for(unsigned int rect = 0; rect < number_rectangles; ++rect)
	 	if( rectangles[rect].distance_from_border(location, buffer) < 1e-8)
	   		return false;

	// check interior lines:
	unsigned int n_lines = lines.size();
	for(unsigned int i = 0; i < n_lines; ++i)
		if( !lines[i].is_in_bounds(location, buffer) )
			return false;

	return true;
}

/** \brief Pushes test point away from boundary by buffer if needed. */
template<int dim>
void 
Geometry<dim>::addPointBuffer(const double buffer, const Point<dim>& test_point,
               Point<dim>& buffered_point) const
{
	buffered_point = test_point;

	// add buffer to spheres-- assuming already outside of actual sphere:
	unsigned int number_spheres = spheres.size();
	for(unsigned int sphere_id = 0; sphere_id < number_spheres; ++sphere_id)
		if(spheres[sphere_id].isInSphere(test_point, buffer))
		{
			Tensor<1,dim> normal = spheres[sphere_id].getSphereNormalVector(test_point);
			buffered_point = buffered_point + buffer*normal;
			break;
		}

	// add buffer to rectangles:
	unsigned int number_rectangles = rectangles.size();
	for(unsigned int rect_id = 0; rect_id < number_rectangles; ++rect_id)
		if( rectangles[rect_id].distance_from_border(buffered_point) < buffer )
		{
			/** @todo double check this */
			// add buffer to appropriate direction ...
			Tensor<1, dim> normal_vector = rectangles[rect_id].getNormalVector(buffered_point);
			buffered_point += (buffer*normal_vector);
			break;
		} 

	// add buffer to lines:
	unsigned int n_lines = lines.size();
	for(unsigned int i = 0; i < n_lines; ++i)
		if( lines[i].distance_from_line(buffered_point) < buffer )
		{
			Tensor<1, dim> normal_vector = lines[i].getNormalVector(buffered_point);
			buffered_point += (buffer*normal_vector);
			break;
		} 
}

/** \brief Returns vector of points in domain at specified resolution.
* Also gives points around any obstacles are a higher resolution.
* @todo generalize to 3D
*/
template<int dim>
std::vector<Point<dim> > 
Geometry<dim>::getQuerryPoints(double resolution) const
{
	if(dim != 2)
		throw std::invalid_argument("getQuerryPoints() not implemented for dim != 2");
	
	const unsigned int number_spheres = spheres.size();
	const unsigned int circlePoints = 32; // can maybe pass this in too

	unsigned int gridPoints = 1;
	for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
		gridPoints *= ceil( (this->getWidth(dim_itr))/resolution );

	std::vector<Point<dim> > querry_points;
	querry_points.reserve(gridPoints + circlePoints*number_spheres);

	for(unsigned int i = 0; i < spheres.size(); i++)
	{
		const double radius = spheres[i].getRadius();
		const Point<dim> center = spheres[i].getCenter();

		double theta = 0; 
		for(unsigned int j = 0; j < circlePoints; j++)
		{
			const double x = radius*std::cos(theta) + center[0]; 
			// @todo *** may need to add eps to be in domain
			const double y = radius*std::sin(theta) + center[1];

			const Point<2> p(x,y);
			if( this->isInDomain(p) )
			{
				querry_points.push_back(p);
			} // if point in domain

			theta += dealii::numbers::PI/16.0;
		} // for circle points 
	} // for spheres


	Point<dim> p = bottom_left;
	for(unsigned int i=0; i<gridPoints; i++)
	{
		if( this->isInDomain(p) )
			querry_points.push_back(p);

		// update p:
		p[0] += resolution;

		if(p[0] > top_right[0])
		{
			p[0] = bottom_left[0];
			p[1] += resolution;
		}
	} 

	return querry_points;
} // getQuerryPoints()

/** \brief Output informaiton about this geometry object */
template<int dim>
void 
Geometry<dim>::printInfo(std::ostream& out) const
{
	out << std::endl << std::endl << Utility::medium_line << std::endl
		<< "\t\tGEOMETRY INFO:"
		<< std::endl << Utility::medium_line << std::endl;

	out << "DIMENSION: " << dim << std::endl
		<< "\t BottomLeft: " << bottom_left << std::endl
		<< "\t TopRight: " << top_right << std::endl;

	out << "\nBOUNDARY CONDITIONS: " << std::endl;
	for(unsigned int i=0; i < dim; i++)
		out << "\t " << i << "th boundary: " 
			<< getBoundaryConditionString(boundary_conditions[i]) << std::endl; // enum to ostream?

	out << "\nSPHERES: " << spheres.size() << std::endl;
	for(unsigned int i = 0; i < spheres.size(); i++)
		out << "\t Sphere Center: " << spheres[i].getCenter() 
			<< " \tRadius: " << spheres[i].getRadius() << std::endl
			<< " \tType: " 
				<< getObstacleTypeString(spheres[i].getObstacleType()) << std::endl;

	out << "\nRECTANGLES: " << rectangles.size() << std::endl;
	for(unsigned int i = 0; i < rectangles.size(); i++)
		out << "\t Bottom left: " << rectangles[i].getBottomLeft() 
			<< "\n\t Top right: " << rectangles[i].getTopRight() << std::endl;

	out << "\nLINES: " << lines.size() << std::endl;
	for(unsigned int i = 0; i < lines.size(); ++i)
		lines[i].print(out);

	out << Utility::medium_line << std::endl
			<< std::endl << std::endl;
}

/** \brief Output geometry information to file to assist post-processing */
template<int dim>
void 
Geometry<dim>::outputGeometry(std::string output_directory) const
{
	// general info:
	std::ofstream geo_out(output_directory + "/geometryInfo.dat");
	this->printInfo(geo_out);

	// boundary:
	std::ofstream boundary_out(output_directory + "/boundary.dat");
	boundary_out << bottom_left << std::endl
		<< top_right << std::endl;

	// spheres:
	std::ofstream spheres_out(output_directory + "/spheres.dat");
	for(unsigned int i = 0; i < spheres.size(); ++i)
		spheres_out << spheres[i].getCenter() << " " 
			<< spheres[i].getRadius() << std::endl;

	// rectangles:
	std::ofstream rect_out(output_directory + "/rectangles.dat");
	for(unsigned int i = 0; i < rectangles.size(); ++i)
		rect_out << rectangles[i].getBottomLeft() << " " 
			<< rectangles[i].getTopRight() << std::endl;

	// lines:
	std::ofstream lines_out(output_directory + "/lines.dat");
	for(unsigned int i = 0; i < lines.size(); ++i)
		lines[i].print(lines_out); 
} 


// LEGACY:

// FILE INIT:

// template<int dim>
// void Geometry<dim>::initialize(std::string geometryFile, std::string meshFile)
// {
//    std::cout << "...Initializing from Geometry File: " << geometryFile << std::endl;

//    mesh_file = meshFile;

//    std::ifstream infile(geometryFile);

//    if(!infile)
//      throw std::invalid_argument("ERROR: GEOMETRY FILE DOES NOT EXIST");

//    std::string line;
//    std::string delimiter = " ";

//    bool usingDefaultDiscretization = true;
//    // scale = 1;

//    /// FILTER PARAMETERS:
//    unsigned int number_channels;
//    double wall_thickness;
//    double channel_thickness;
//    double filter_left;
//    double filter_right;
//    double filter_center;

//    // Input Format: variable value \n
//    while(std::getline(infile,line)){
//    //  std::cout << line << std::endl;
//    //  std::cout << std::endl;

//      // "tokenize" line:
//      size_t pos = 0;
//      std::string token;
//      while((pos = line.find(delimiter)) != std::string::npos){
//        token = line.substr(0,pos);
//   //    std::cout << "token is " << token << std::endl;
//        if(token.compare("Domain") == 0){
//          std::istringstream numRead(line);
//          std::string category;
//          unsigned int numLines; 
//          numRead >> category >> numLines; // get number of lines to read for boundary
//          for(unsigned int i = 0; i < numLines; i++){
//            std::getline(infile,line);
//            std::istringstream stream(line);
//            std::string varName;
//            double value; 
//            stream >> varName;
//            if(varName.compare("bottom_left") == 0)
//            {
//              Point<dim> temp;
//              for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
//              {
//                stream >> value;
//                temp[dim_itr] = value;  
//              }  
//              bottom_left = temp;
//            }
//            if(varName.compare("top_right") == 0)
//            {
//              Point<dim> temp;
//              for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
//              {
//                stream >> value;
//                temp[dim_itr] = value;  
//              }  
//              top_right = temp;
//            }
//          } // for domain boundary
//          std::getline(infile,line); // move to next line
//        } // read in boundary lines
//        else if(token.compare("Boundaries") == 0)
//        {
//          std::istringstream lineRead(line);
//          std::string category;
//          lineRead >> category;
//          int value;
//          for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
//          {
//            lineRead >> value; 
//            boundary_conditions[dim_itr] = (BoundaryCondition)value;  
//          }
//          // move to next line
//          std::getline(infile,line); 
//        }
//        else if(token.compare("Discretization") == 0)
//        {
//          usingDefaultDiscretization = false;
//          std::istringstream lineRead(line);
//          std::string category;
//          lineRead >> category;
//          int value;
//          for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
//          {
//            lineRead >> value; 
//            discretization[dim_itr] = (unsigned int)value;  
//          }
//          // move to next line
//          std::getline(infile,line); 
//        }
//        else if(token.compare("Scale") == 0)
//        {
//          std::istringstream lineRead(line);
//          std::string category;
//          lineRead >> category;
//          double value;
//          lineRead >> value;
//          // scale = value;
//          // move to next line
//          std::getline(infile,line); 
//        }
//        else if(token.compare("Spheres") == 0)
//        {
//          std::istringstream numRead(line);
//          std::string category;
//          unsigned int numCircles; 
//          numRead >> category >> numCircles; // get number of lines to read for boundary
//          for(unsigned int i = 0; i < numCircles; i++){
//            std::getline(infile,line);
//            std::istringstream stream(line);
//            double value, radius;
//            Point<dim> center;
//            for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
//            {
//              stream >> value;
//              center[dim_itr] = value;
//            }
//            stream >> radius;
//            spheres.push_back(Sphere<dim>(center,radius));
//          } // for boundary lines
//          // move to next line
//          std::getline(infile,line); 
//        } // read in circle lines
//        else if(token.compare("Rectangles") == 0){
//          std::istringstream numRead(line);
//          std::string category;
//          unsigned int numRectangles; 
//          numRead >> category >> numRectangles; // get number of lines to read for boundary
//          for(unsigned int i = 0; i < numRectangles; i++){
//            std::getline(infile,line);
//            std::istringstream stream(line);
//            Point<dim> lower, upper;
//            double value;
//            for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
//            {
//              stream >> value;
//              lower[dim_itr] = value;
//            }
//           for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
//            {
//              stream >> value;
//              upper[dim_itr] = value;
//            }
//            rectangles.push_back(HyperRectangle<dim>(lower,upper));
//          } // for boundary lines
//          std::getline(infile,line); // move to next line
//        } // read in rectangles
//        else if(token.compare("Mesh") == 0){
//          std::istringstream stream(line);
//          std::string category;
//          int mType;
//          stream >> category >> mType;
//          mesh_type = (MeshType)mType;

//          // move to next line
//          std::getline(infile,line); 
//        } // read in mesh file
//        else if(token.compare("number_channels") ==0){
//          std::cout << "reading number channels... " << std::endl;
//          std::istringstream inStream(line);
//          std::string category;
//          inStream >> category;
//          inStream >> number_channels;
//          std::getline(infile,line);
//        }
//        else if(token.compare("channel_thickness") ==0){
//          std::cout << "reading channel_thickness... " << std::endl;
//          std::istringstream inStream(line);
//          std::string category;
//          inStream >> category;
//          inStream >> channel_thickness;
//          std::getline(infile,line);
//        }
//        else if(token.compare("wall_thickness") ==0){
//          std::cout << "reading wall wall_thickness... " << std::endl;
//          std::istringstream inStream(line);
//          std::string category;
//          inStream >> category;
//          inStream >> wall_thickness;
//          std::getline(infile,line);
//        }
//        else if(token.compare("left_length") ==0){
//          std::istringstream inStream(line);
//          std::string category;
//          inStream >> category;
//          inStream >> filter_left;
//          std::getline(infile,line);
//        }
//        else if(token.compare("center_length") ==0){
//          std::istringstream inStream(line);
//          std::string category;
//          inStream >> category;
//          inStream >> filter_center;
//          std::getline(infile,line);
//        }
//        else if(token.compare("right_length") ==0){
//          std::istringstream inStream(line);
//          std::string category;
//          inStream >> category;
//          inStream >> filter_right;
//          std::getline(infile,line);
//        }
//        else{
//          line.erase(0,pos + delimiter.length()); // otherwise might be infinite loop
//        } // otherwise keep parsing

//      } // while tokenizing line
//   //   std:: cout << line << "\n\n";

//    }  // while reading lines


//    if(mesh_type == MeshType::FILTER)
//    {
//      create_filter_geometry(number_channels,
//                            channel_thickness,
//                            wall_thickness,
//                            filter_left,
//                            filter_center,
//                            filter_right);
//    }


//    if(usingDefaultDiscretization)
//    {
//      for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
//        discretization[dim_itr] = (unsigned int) round(5*
//          (top_right[dim_itr] - bottom_left[dim_itr]) );
//    }

//    // if(scale != 1)
//    //   rescale();
// } // initialize() -- probably want to break/clean this up...

} // CLOSE NAMESPACE
#endif 