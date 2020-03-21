#ifndef GEOMETRY_H
#define GEOMETRY_H

// TO DO:
// have sphere implement its own reflection
// clean up names and intialization
// do we need all the below headers??? -- remove uncessary ones
// remove scale parameter

#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cstddef>
#include <array>

#include "sphere.h"
#include "hyper_rectangle.h"

#include "../utility/enum_types.h" // clean up names and definitions...
#include "../utility/parameter_handler.h"

// maybe remove this:
#include "./geo_types.h" 


namespace MicrobeSimulator{

template<int dim>
class Geometry{
public:

	// CONSTRUCTORS:
	Geometry();
	Geometry(const Point<dim>& lower, 
	         const Point<dim>& upper,
	         const std::array<BoundaryCondition, dim>& bcs);  

	// Geometry(const Geometry &geo); // copy constructor

	static void declare_parameters(ParameterHandler& prm);

	void init(const ParameterHandler& prm);

	// LEGACY INITIALIZATION:
	void initialize(std::string geoFile, std::string meshFile = "");
	void initialize(GeoTypes::Filter filter);
	void initialize(GeoTypes::Mixer mixer); 
	void initialize(GeoTypes::Pipe pipe); 
	void initialize(GeoTypes::Splitter splitter); 

	// create_filter_geometry(unsigned int number_channels,
	//                             double channel_thickness,
	//                             double wall_thickness,
	//                             double filter_left,
	//                             double filter_center,
	//                             double filter_right);


	// OPERATORS:
	Geometry<dim>& operator=(const Geometry& rhs);

	template<int dimension>
	friend std::ostream& operator<<(std::ostream& out, const Geometry<dimension>& geo);

	// accessors:
	Point<dim> getBottomLeftPoint() const;
	Point<dim> getTopRightPoint() const;
	double getScaleFactor() const;
	std::array<unsigned int, dim> getDiscretization() const;
	double getWidth(unsigned int direction) const;
	std::array<BoundaryCondition, dim> getBoundaryConditions() const;
	MeshType getMeshType() const;
	std::string getMeshFile() const;

	/// sphere accessors:
	std::vector<Sphere<dim> > getSpheres() const;
	unsigned int getNumberSpheres() const;
	Sphere<dim> getSphereAt(unsigned int i) const;

	/// rectangle accessors:
	std::vector<HyperRectangle<dim> > getRectangles() const;
	unsigned int getNumberRectangles() const;
	HyperRectangle<dim> getRectangleAt(unsigned int i) const;

	// mutators:
	void setBottomLeftPoint(const Point<dim>& lower);
	void setTopRightPoint(const Point<dim>& upper);
	void setScaleFactor(double s);
	void setBoundaryConditions(const std::array<BoundaryCondition, dim>& bcs);
	void setMeshType(MeshType mtype);
	void addSphere(const Sphere<dim>& sp);
	void addRectangle(const HyperRectangle<dim>& rect);

	void rescale();

	// interface:
	// functions:
	void checkBoundaries(const Point<dim>& oldPoint, 
	                   Point<dim>& newPoint,
	                   const double buffer = 0.005) const; // only function bacteria uses***

	bool isInDomain(const Point<dim>& location) const;

	void addPointBuffer(const double buffer,
	                   const Point<dim>& test_point,
	                   Point<dim>& buffered_point) const; 

	void printInfo(std::ostream& out) const;
	void outputGeometry(std::string output_directory = ".") const;

	// for debugging:
	std::vector<Point<dim> > getQuerryPoints(double resolution = 0.2) const; 
private:
	// bounding box:
	Point<dim> bottom_left;
	Point<dim> top_right;

	double scale; // (legacy)

	std::array<BoundaryCondition, dim> boundary_conditions; 

	// for discrete field: (legacy)
	std::array<unsigned int, dim> discretization;

	// bounding circle or cylinder
	// double boundaryRadius; // perhaps change to just have ``interior obstacle...''

	// obstacles:
	std::vector<Sphere<dim> > spheres;
	std::vector<HyperRectangle<dim> > rectangles;
	// cylinders
	// ellipsoids ...
	// ... each can be bounding or internal....

	MeshType mesh_type;
	std::string mesh_file;

	// reflecting off spheres: 
	/// these functions should actually be implemented in the sphere class...
	void reflectSphere(const unsigned int sphere_id, 
	              const Point<dim>& oldPoint, 
	              Point<dim>& newPoint,
	               const double buffer = 0.) const;

	bool isInSphere(const unsigned int sphere_id, 
	            const Point<dim>& location,
	            const double buffer = 0.) const;

	Point<dim> getLineSphereIntersection(const unsigned int sphere_id, 
	                                const Point<dim>& oldPoint, 
	                                const Point<dim>& newPoint, 
	                                const double buffer = 0.) const;

	Tensor<1, dim> getSphereNormalVector(const unsigned int sphere_id, 
	                                const Point<dim>& intersection_point
	                                ) const;

	static void declare_filter_parameters(ParameterHandler& prm);
	static void declare_mixer_parameters(ParameterHandler& prm);
	static void declare_splitter_parameters(ParameterHandler& prm);

	void create_filter_geometry(unsigned int number_channels,
	                   double channel_thickness,
	                   double wall_thickness,
	                   double filter_left,
	                   double filter_center,
	                   double filter_right);

	void create_mixer_geometry(double left_length,
	                   double right_length,
	                   double height,
	                   double radius);

	void create_splitter_geometry(double left_length,
	                   double right_length,
	                   double height,
	                  double radius);
}; // class Geometry{}

// IMPLEMENTATION
// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------


// CONSTRUCTORS:
template<int dim>
Geometry<dim>::Geometry()
	:
	scale(1)
{}


template<int dim>
Geometry<dim>::Geometry(const Point<dim>& lower, 
 const Point<dim>& upper, const std::array<BoundaryCondition, dim>& bcs)
	:
	bottom_left(lower),
	top_right(upper),
	scale(1),
	boundary_conditions(bcs)
{
	for(unsigned int i = 0; i < dim; i++)
		if(bottom_left[i] > top_right[i])
			throw std::invalid_argument("Geometry error: bottom left point cannot" 
				"be greater than top right.");
}

template<>
void
Geometry<2>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.declare_entry("Geometry type",
		          "Box",
		          Patterns::Selection("Box|Filter|Mixer|Splitter|File"),
		          "Geometry type. Options are Box, Filter, Mixer "
		          "Splitter, of File. For File, name of file must also be"
		          "provided in the \"Geometry file\" parameter. ");
		prm.declare_entry("Geometry file",
		          "",
		          Patterns::Anything(),
		          "Name of geometry file to initialize geometry from. "
		          "Mesh type of geometry file should match mesh type "
		          "given for grid construction.");
		prm.declare_entry("Bottom left","{0,0}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Top right","{5,5}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Boundary conditions",
		          "{WRAP,WRAP}",
		          Patterns::List(Patterns::Selection("WRAP|REFLECT|OPEN")),
		          "Boundary conditions.");
		declare_filter_parameters(prm);
		declare_mixer_parameters(prm);
		declare_splitter_parameters(prm);
	prm.leave_subsection();
}

template<>
void
Geometry<3>::declare_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Geometry");
		prm.declare_entry("Geometry type",
		          "Box",
		          Patterns::Selection("Box|Filter|Mixer|Splitter|File"),
		          "Geometry type. Options are Box, Filter, Mixer "
		          "Splitter, of File. For File, name of file must also be"
		          "provided in the \"Geometry file\" parameter. ");
		prm.declare_entry("Geometry file",
		          "",
		          Patterns::Anything(),
		          "Name of geometry file to initialize geometry from. "
		          "Mesh type of geometry file should match mesh type "
		          "given for grid construction.");
		prm.declare_entry("Bottom left","{0,0,0}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Top right","{5,5,5}",Patterns::List(Patterns::Double()));
		prm.declare_entry("Boundary conditions",
		          "{WRAP,WRAP,WRAP}",
		          Patterns::List(Patterns::Selection("WRAP|REFLECT|OPEN")),
		          "Boundary conditions.");
		declare_filter_parameters(prm);
		declare_mixer_parameters(prm);
		declare_splitter_parameters(prm);
	prm.leave_subsection();
}

template<int dim>
void 
Geometry<dim>::declare_filter_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Filter");
		prm.declare_entry("Number channels","1",Patterns::Double());
		prm.declare_entry("Channel thickness","1",Patterns::Double());
		prm.declare_entry("Wall thickness","1",Patterns::Double());
		prm.declare_entry("Left length","1",Patterns::Double());
		prm.declare_entry("Center length","1",Patterns::Double());
		prm.declare_entry("Right length","1",Patterns::Double());
	prm.leave_subsection();
}

template<int dim>
void 
Geometry<dim>::declare_mixer_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Mixer");
		prm.declare_entry("Left length","1",Patterns::Double());
		prm.declare_entry("Right length","1",Patterns::Double());
		prm.declare_entry("Height","1",Patterns::Double());
		prm.declare_entry("Radius","1",Patterns::Double());
	prm.leave_subsection();
}

template<int dim>
void 
Geometry<dim>::declare_splitter_parameters(ParameterHandler& prm)
{
	prm.enter_subsection("Splitter");
		prm.declare_entry("Left length","1",Patterns::Double());
		prm.declare_entry("Right length","1",Patterns::Double());
		prm.declare_entry("Height","1",Patterns::Double());
		prm.declare_entry("Radius","1",Patterns::Double());
	prm.leave_subsection();
}

template<int dim>
void 
Geometry<dim>::init(const ParameterHandler& prm)
{
	std::string section = "Geometry";
	std::string geo_type_str = prm.get_string(section, "Geometry type");
	MeshType geo_type = stringToMeshType(geo_type_str); 

	if(geo_type == MeshType::BOX_MESH)
	{
	 // dimension should be checked by parser
	 std::vector<double> bl = prm.get_double_vector(section, "Bottom left");
	 std::vector<double> tr = prm.get_double_vector(section, "Top right");
	 std::vector<std::string> bnd = prm.get_string_vector(section,"Boundary conditions");

	 mesh_type = MeshType::BOX_MESH;
	 for(unsigned int i = 0; i < dim; ++i)
	 {
	   bottom_left[i] = bl[i];
	   top_right[i] = tr[i];

	  boundary_conditions[i] = stringToBoundaryCondition(bnd[i]);
	 }
	}
	else if(geo_type == MeshType::FILTER)
	{
	 std::string subsection = section + ".Filter";

	 create_filter_geometry(prm.get_unsigned(subsection, "Number channels"),
	                         prm.get_double(subsection, "Channel thickness"),
	                         prm.get_double(subsection, "Wall thickness"),
	                         prm.get_double(subsection, "Left length"),
	                         prm.get_double(subsection, "Center length"),
	                         prm.get_double(subsection, "Right length") );

	 mesh_type = MeshType::FILTER;

	 boundary_conditions[0] = BoundaryCondition::OPEN;
	 boundary_conditions[1] = BoundaryCondition::REFLECT;
	}
	else if(geo_type == MeshType::MIXER)
	{
	 std::string subsection = section + ".Mixer";

	 create_mixer_geometry(prm.get_double(subsection, "Left length"),
	                 prm.get_double(subsection, "Right length"),
	                 prm.get_double(subsection, "Height"),
	                 prm.get_double(subsection, "Radius") );

	 mesh_type = MeshType::MIXER;

	 boundary_conditions[0] = BoundaryCondition::OPEN;
	 boundary_conditions[1] = BoundaryCondition::REFLECT;
	}
	else if(geo_type == MeshType::SPLITTER)
	{
		std::string subsection = section + ".Splitter";

	 create_splitter_geometry(prm.get_double(subsection, "Left length"),
	                   prm.get_double(subsection, "Right length"),
	                   prm.get_double(subsection, "Height"),
	                   prm.get_double(subsection, "Radius") );

	 mesh_type = MeshType::SPLITTER;

	 boundary_conditions[0] = BoundaryCondition::OPEN;
	 boundary_conditions[1] = BoundaryCondition::REFLECT;
	}
	else if(geo_type == MeshType::FILE_MESH)
	{
		// making use of legacy function for now:
		initialize(prm.get_string(section, "Geometry file"), 
					prm.get_string("Mesh", "Mesh file") );
	}
	else
	{
	 std::cout << "GEOMETRY TYPE NOT IMPLEMENTED" << std::endl;
	}
}

// INITIALIZATON:

template<int dim>
void
Geometry<dim>::initialize(GeoTypes::Filter filter)
{
 create_filter_geometry(filter.number_channels,
                         filter.channel_thickness,
                         filter.wall_thickness,
                         filter.left_length,
                         filter.center_length,
                         filter.right_length);

 mesh_type = MeshType::FILTER;

 boundary_conditions[0] = BoundaryCondition::OPEN;
 boundary_conditions[1] = BoundaryCondition::REFLECT;
}

template<int dim>
void
Geometry<dim>::initialize(GeoTypes::Mixer mixer)
{
 create_mixer_geometry(mixer.left_length,
                       mixer.right_length,
                       mixer.height,
                       mixer.radius);

 mesh_type = MeshType::MIXER;

 boundary_conditions[0] = BoundaryCondition::OPEN;
 boundary_conditions[1] = BoundaryCondition::REFLECT;
}

template<int dim>
void 
Geometry<dim>::initialize(GeoTypes::Splitter splitter)
{
 create_splitter_geometry(splitter.left,
                         splitter.right,
                         splitter.height,
                         splitter.radius);

 mesh_type = MeshType::SPLITTER;

 boundary_conditions[0] = BoundaryCondition::OPEN;
 boundary_conditions[1] = BoundaryCondition::REFLECT;
}

template<int dim>
void 
Geometry<dim>::initialize(GeoTypes::Pipe pipe)
{
if(dim != 2)
 throw std::runtime_error("geometry pipe initialization not implemented for 2d");

bottom_left = Point<2>(pipe.xmin,pipe.ymin);
top_right = Point<2>(pipe.xmax,pipe.ymax);

mesh_type = MeshType::BOX_MESH;

boundary_conditions[0] = BoundaryCondition::OPEN;
boundary_conditions[1] = BoundaryCondition::REFLECT; 
}

// ASSIGNMENT OPERATORS:
template<int dim>
Geometry<dim>& Geometry<dim>::operator=(const Geometry& rhs)
{
// check for self copy:
if(this == &rhs)
 return *this;

// copy:
bottom_left = rhs.bottom_left;
top_right = rhs.top_right;
scale = rhs.scale;
boundary_conditions = rhs.boundary_conditions; // should work without loop? 
spheres = rhs.spheres;
rectangles = rhs.rectangles;

return *this;
}


template<int dim>
std::ostream& operator<<(std::ostream& out, const Geometry<dim>& geo)
{
out << std::endl << "Dimension: " << dim << std::endl
 << "BottomLeft: " << geo.bottom_left << std::endl
 << "TopRight: " << geo.top_right << std::endl;

out << "Boundary Conditions: " << std::endl;
for(unsigned int i=0; i < dim; i++)
 out << i << "th boundary: " << geo.boundary_conditions[i] << std::endl; // enum to ostream?

// HERE, overload << for spheres instead...
out << "Number of Spheres: " << geo.spheres.size() << std::endl;
for(unsigned int i = 0; i < geo.spheres.size(); i++)
 out << "\t Sphere Center: " << geo.spheres[i].getCenter() 
   << " Radius: " << geo.spheres[i].getRadius() << std::endl;

return out;
}



// ACCESSORS:
template<int dim>
Point<dim> Geometry<dim>::getBottomLeftPoint() const
{
return bottom_left;
}

template<int dim>
Point<dim> Geometry<dim>::getTopRightPoint() const
{
return top_right;
}

template<int dim>
double Geometry<dim>::getScaleFactor() const
{
return scale;
}

template<int dim>
double Geometry<dim>::getWidth(unsigned int direction) const
{
if(direction >= dim)
 throw std::invalid_argument("Desired dimension to get width does not exist");
return top_right[direction] - bottom_left[direction];
}


template<int dim>
std::array<BoundaryCondition, dim> Geometry<dim>::getBoundaryConditions() const
{
return boundary_conditions;
}


template<int dim>
std::array<unsigned int, dim> Geometry<dim>::getDiscretization() const
{
return discretization;
}

template<int dim>
MeshType Geometry<dim>::getMeshType() const
{
return mesh_type;
}

template<int dim>
std::string Geometry<dim>::getMeshFile() const
{
return mesh_file;
}

template<int dim>
std::vector<Sphere<dim> > Geometry<dim>::getSpheres() const
{
return spheres;
}

template<int dim>
unsigned int Geometry<dim>::getNumberSpheres() const
{
return spheres.size();
}

template<int dim>
Sphere<dim> Geometry<dim>::getSphereAt(unsigned int i) const
{
return spheres[i];
}


/// rectangle accessors:

template<int dim>    
std::vector<HyperRectangle<dim> > 
Geometry<dim>::getRectangles() const
{
return rectangles;
}

template<int dim>  
unsigned int 
Geometry<dim>::getNumberRectangles() const
{
return rectangles.size();
}

template<int dim>
HyperRectangle<dim> Geometry<dim>::getRectangleAt(unsigned int i) const
{
return rectangles[i];
}

// MUTATORS:
template<int dim>
void Geometry<dim>::setBottomLeftPoint(const Point<dim>& lower)
{
bottom_left = lower;
}

template<int dim>
void Geometry<dim>::setTopRightPoint(const Point<dim>& upper)
{
top_right = upper;
}

template<int dim>
void Geometry<dim>::setScaleFactor(double s)
{
scale = s;
}

template<int dim>
void Geometry<dim>::setBoundaryConditions(const std::array<BoundaryCondition, dim>& bcs)
{
boundary_conditions = bcs;
}

template<int dim>
void Geometry<dim>::setMeshType(MeshType mtype)
{
mesh_type = mtype;
}

template<int dim>
void Geometry<dim>::addSphere(const Sphere<dim>& sp)
{
spheres.push_back(sp);
}

template<int dim>
void Geometry<dim>::addRectangle(const HyperRectangle<dim>& rect)
{
rectangles.push_back(rect);
}

template<int dim>
void Geometry<dim>::initialize(std::string geometryFile, std::string meshFile)
{
   std::cout << "...Initializing from Geometry File: " << geometryFile << std::endl;

   mesh_file = meshFile;

   std::ifstream infile(geometryFile);

   if(!infile)
     throw std::invalid_argument("ERROR: GEOMETRY FILE DOES NOT EXIST");

   std::string line;
   std::string delimiter = " ";

   bool usingDefaultDiscretization = true;
   scale = 1;

   /// FILTER PARAMETERS:
   unsigned int number_channels;
   double wall_thickness;
   double channel_thickness;
   double filter_left;
   double filter_right;
   double filter_center;

   // Input Format: variable value \n
   while(std::getline(infile,line)){
   //  std::cout << line << std::endl;
   //  std::cout << std::endl;

     // "tokenize" line:
     size_t pos = 0;
     std::string token;
     while((pos = line.find(delimiter)) != std::string::npos){
       token = line.substr(0,pos);
  //    std::cout << "token is " << token << std::endl;
       if(token.compare("Domain") == 0){
         std::istringstream numRead(line);
         std::string category;
         unsigned int numLines; 
         numRead >> category >> numLines; // get number of lines to read for boundary
         for(unsigned int i = 0; i < numLines; i++){
           std::getline(infile,line);
           std::istringstream stream(line);
           std::string varName;
           double value; 
           stream >> varName;
           if(varName.compare("bottom_left") == 0)
           {
             Point<dim> temp;
             for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
             {
               stream >> value;
               temp[dim_itr] = value;  
             }  
             bottom_left = temp;
           }
           if(varName.compare("top_right") == 0)
           {
             Point<dim> temp;
             for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
             {
               stream >> value;
               temp[dim_itr] = value;  
             }  
             top_right = temp;
           }
         } // for domain boundary
         std::getline(infile,line); // move to next line
       } // read in boundary lines
       else if(token.compare("Boundaries") == 0)
       {
         std::istringstream lineRead(line);
         std::string category;
         lineRead >> category;
         int value;
         for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
         {
           lineRead >> value; 
           boundary_conditions[dim_itr] = (BoundaryCondition)value;  
         }
         // move to next line
         std::getline(infile,line); 
       }
       else if(token.compare("Discretization") == 0)
       {
         usingDefaultDiscretization = false;
         std::istringstream lineRead(line);
         std::string category;
         lineRead >> category;
         int value;
         for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
         {
           lineRead >> value; 
           discretization[dim_itr] = (unsigned int)value;  
         }
         // move to next line
         std::getline(infile,line); 
       }
       else if(token.compare("Scale") == 0)
       {
         std::istringstream lineRead(line);
         std::string category;
         lineRead >> category;
         double value;
         lineRead >> value;
         scale = value;
         // move to next line
         std::getline(infile,line); 
       }
       else if(token.compare("Spheres") == 0)
       {
         std::istringstream numRead(line);
         std::string category;
         unsigned int numCircles; 
         numRead >> category >> numCircles; // get number of lines to read for boundary
         for(unsigned int i = 0; i < numCircles; i++){
           std::getline(infile,line);
           std::istringstream stream(line);
           double value, radius;
           Point<dim> center;
           for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
           {
             stream >> value;
             center[dim_itr] = value;
           }
           stream >> radius;
           spheres.push_back(Sphere<dim>(center,radius));
         } // for boundary lines
         // move to next line
         std::getline(infile,line); 
       } // read in circle lines
       else if(token.compare("Rectangles") == 0){
         std::istringstream numRead(line);
         std::string category;
         unsigned int numRectangles; 
         numRead >> category >> numRectangles; // get number of lines to read for boundary
         for(unsigned int i = 0; i < numRectangles; i++){
           std::getline(infile,line);
           std::istringstream stream(line);
           Point<dim> lower, upper;
           double value;
           for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
           {
             stream >> value;
             lower[dim_itr] = value;
           }
          for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
           {
             stream >> value;
             upper[dim_itr] = value;
           }
           rectangles.push_back(HyperRectangle<dim>(lower,upper));
         } // for boundary lines
         std::getline(infile,line); // move to next line
       } // read in rectangles
       else if(token.compare("Mesh") == 0){
         std::istringstream stream(line);
         std::string category;
         int mType;
         stream >> category >> mType;
         mesh_type = (MeshType)mType;

         // move to next line
         std::getline(infile,line); 
       } // read in mesh file
       else if(token.compare("number_channels") ==0){
         std::cout << "reading number channels... " << std::endl;
         std::istringstream inStream(line);
         std::string category;
         inStream >> category;
         inStream >> number_channels;
         std::getline(infile,line);
       }
       else if(token.compare("channel_thickness") ==0){
         std::cout << "reading channel_thickness... " << std::endl;
         std::istringstream inStream(line);
         std::string category;
         inStream >> category;
         inStream >> channel_thickness;
         std::getline(infile,line);
       }
       else if(token.compare("wall_thickness") ==0){
         std::cout << "reading wall wall_thickness... " << std::endl;
         std::istringstream inStream(line);
         std::string category;
         inStream >> category;
         inStream >> wall_thickness;
         std::getline(infile,line);
       }
       else if(token.compare("left_length") ==0){
         std::istringstream inStream(line);
         std::string category;
         inStream >> category;
         inStream >> filter_left;
         std::getline(infile,line);
       }
       else if(token.compare("center_length") ==0){
         std::istringstream inStream(line);
         std::string category;
         inStream >> category;
         inStream >> filter_center;
         std::getline(infile,line);
       }
       else if(token.compare("right_length") ==0){
         std::istringstream inStream(line);
         std::string category;
         inStream >> category;
         inStream >> filter_right;
         std::getline(infile,line);
       }
       else{
         line.erase(0,pos + delimiter.length()); // otherwise might be infinite loop
       } // otherwise keep parsing

     } // while tokenizing line
  //   std:: cout << line << "\n\n";

   }  // while reading lines


   if(mesh_type == MeshType::FILTER)
   {
     create_filter_geometry(number_channels,
                           channel_thickness,
                           wall_thickness,
                           filter_left,
                           filter_center,
                           filter_right);
   }


   if(usingDefaultDiscretization)
   {
     for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
       discretization[dim_itr] = (unsigned int) round(5*
         (top_right[dim_itr] - bottom_left[dim_itr]) );
   }

   if(scale != 1)
     rescale();
} // initialize() -- probably want to break/clean this up...


template<int dim>
void 
Geometry<dim>::create_filter_geometry(unsigned int number_channels,
                                 double channel_thickness,
                                 double wall_thickness,
                                 double filter_left,
                                 double filter_center,
                                 double filter_right)
{
std::cout << "\n\t CREATING FILTER GEOMETRY: " << std::endl
 << "number_channels = " << number_channels << std::endl
 << "channel_thickness = " << channel_thickness << std::endl
 << "wall_thickness = " << wall_thickness << std::endl
 << "filter_left = " << filter_left << std::endl
 << "filter_center = " << filter_center << std::endl
 << "filter_right = " << filter_right << std::endl;

if(dim != 2)
 throw std::runtime_error("filter geometry currently only implemented for dim == 2");

// check parameters...
const unsigned int num_rect = number_channels - 1;

const double width = filter_left + filter_center + filter_right;
const double height = number_channels*channel_thickness + num_rect*wall_thickness;

/// actually, change below, get corners from filter parameters: -- maybe display warning
/*
if( ((top_right[0] - bottom_left[0]) - width) > 1e-8)
 std::cout << "WARNING: Filter parameters should give same width as bounding box corners: " << std::endl
     << "\t corner width: " << top_right[0] - bottom_left[0] << " != mesh width: " 
     << width << std::endl
     << "reassigning corners to match filter parameters...." << std::endl;
if( ((top_right[1] - bottom_left[1]) - height) > 1e-8)
 std::cout << "WARNING: Filter parameters should give same height as bounding box corners: " << std::endl
     << "\t corner height: " << top_right[1] - bottom_left[1] << " != mesh height: " 
     << height << std::endl
     << "reassigning corners to match filter parameters...." << std::endl;
if( bottom_left[0] != 0 || bottom_left[1] != 0)
 std::cout << "WARNING: Reassgning corners to place bottom left at origin..." << std::endl;
*/

// reassign anyways:
for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
 bottom_left[dim_itr] = 0;
top_right[0] = width;
top_right[1] = height;

// add rectangles:
const double x_left = filter_left;
const double x_right = filter_left + filter_center;

double y_bottom = channel_thickness;
double y_top = channel_thickness + wall_thickness;

for(unsigned int i = 0; i < num_rect; ++i)
{
 rectangles.push_back(HyperRectangle<2>(Point<2>(x_left,y_bottom),
                                        Point<2>(x_right,y_top)  ) );
 y_bottom += (channel_thickness + wall_thickness);
 y_top += (channel_thickness + wall_thickness);
}

}


template<int dim>
void 
Geometry<dim>::create_mixer_geometry(double left_length,
                               double right_length,
                               double height,
                               double radius)
{
if(dim != 2)
 throw std::invalid_argument("Mixer only implemented for dim == 2");

std::cout << "\n\t CREATE MIXER GEOMETRY:" << std::endl
 << "left_length = " << left_length << std::endl
 << "right_length = " << right_length << std::endl
 << "height = " << height << std::endl
 << "radius = " << radius << std::endl; 

if(height < 2.*radius)
 throw std::invalid_argument("Mixer height must be greater than sphere diameter");

const double width = left_length + 2.*radius + right_length;

for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
 bottom_left[dim_itr] = 0;

top_right[0] = width;
top_right[1] = height;

// add spheres:
const double center_x = left_length + radius;
spheres.push_back(Sphere<2>(Point<2>(center_x, 0.) ,radius));
spheres.push_back(Sphere<2>(Point<2>(center_x, height) ,radius));
}


template<int dim>
void 
Geometry<dim>::create_splitter_geometry(double left_length,
                     double right_length,
                     double height,
                     double radius)
{
// if(dim != 2)
//   throw std::invalid_argument("Mixer only implemented for dim == 2");

std::cout << "\n\t CREATE SPLITTER GEOMETRY:" << std::endl
 << "left_length = " << left_length << std::endl
 << "right_length = " << right_length << std::endl
 << "height = " << height << std::endl
 << "radius = " << radius << std::endl; 

const double width = left_length + 2.*radius + right_length;

for(unsigned int dim_itr = 0; dim_itr < dim; ++dim_itr)
 bottom_left[dim_itr] = 0;

top_right[0] = width;
top_right[1] = height;

// add spheres:
const double center_x = left_length + radius;
const double center_y = 0.5*height;
spheres.push_back(Sphere<2>(Point<2>(center_x, center_y) ,radius));
}


template<int dim>
void Geometry<dim>::rescale()
{
// rescale bounding box:
bottom_left = scale*bottom_left;
top_right = scale*top_right;

// rescale spheres:
for(unsigned int i = 0; i < spheres.size(); i++)
{
 spheres[i].setCenter(scale*spheres[i].getCenter());
 spheres[i].setRadius(scale*spheres[i].getRadius());
}

// rescale hyper_rectangles:
for(unsigned int i = 0; i < rectangles.size(); i++)
{
 rectangles[i].setBottomLeft(rectangles[i].getBottomLeft());
 rectangles[i].setTopRight(rectangles[i].getTopRight());
}
}


// FUNCTIONS:
template<int dim>
void 
Geometry<dim>::checkBoundaries(const Point<dim>& oldPoint, Point<dim>& newPoint,
								const double buffer) const
{
	const double tolerance = 1e-8;

	// check interior spheres:
	unsigned int number_spheres = spheres.size();
	for(unsigned int sphere_id = 0; sphere_id < number_spheres; ++sphere_id)
		// if( (spheres[sphere_id].distance_from_border(newPoint) - buffer) < tolerance)
		// {
		// 	spheres[sphere_id].reflectPoint(oldPoint, newPoint, buffer); 
		// 	break;
		// }
		if(isInSphere(sphere_id,newPoint, buffer))
		{
			reflectSphere(sphere_id,oldPoint,newPoint, buffer); /** @todo move to sphere class */
			break; // assuming don't hit mutliple circles -- otherwise need to do something else
		}

	// check interior rectangles:
	unsigned int number_rectangles = rectangles.size();
	for(unsigned int rect_id = 0; rect_id < number_rectangles; ++rect_id)
		if( (rectangles[rect_id].distance_from_border(newPoint) - buffer) < tolerance ) // add buffer here too!!
		{
			rectangles[rect_id].reflectPoint(oldPoint, newPoint, buffer); // buffer optional
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

template<int dim>
void 
Geometry<dim>::reflectSphere(const unsigned int sphere_id, 
                         const Point<dim>& oldPoint, 
                         Point<dim>& newPoint,
                         const double buffer) const
{
const Point<dim> intersection = getLineSphereIntersection(sphere_id, 
                                                         oldPoint, 
                                                         newPoint,
                                                         buffer);

const Tensor<1,dim> normal = getSphereNormalVector(sphere_id, 
                                                 intersection); 

const Tensor<1,dim> incident = newPoint - oldPoint;
const Tensor<1,dim> transmitted = newPoint - intersection;

Tensor<1,dim> reflected_point;
reflected_point = incident - 2.0*( incident*normal )*normal;

// rescale:
reflected_point *= (transmitted.norm())/(reflected_point.norm());

// recenter (shift vector origin)
newPoint = intersection + reflected_point;
}


template<int dim> 
Point<dim> Geometry<dim>::getLineSphereIntersection(const unsigned int sphere_id, 
                                                 const Point<dim>& oldPoint, 
                                                 const Point<dim>& newPoint,
                                                 const double buffer) const
{
const double radius = spheres[sphere_id].getRadius() + buffer;
const Point<dim> center = spheres[sphere_id].getCenter();
// line origin:
const Point<dim> origin = oldPoint;
// direction of line:  ( line = origin + d*direction)
Tensor<1,dim> direction = newPoint - oldPoint;
direction /= direction.norm(); // unit vector

// Joachimsthal's Equation:
// d = -b +/- sqrt[ b^2 - c] ... a == 1
const double b = direction*( origin - center );
const double c = (origin - center)*(origin - center) - radius*radius;

const double discriminant = b*b - c;

if(discriminant < 0)
 throw std::runtime_error("Error: Line does not intersect sphere");

if(discriminant == 0)
 return origin + (-b)*direction;

const Point<dim> first_intersection = origin + (-b + std::sqrt(discriminant))*direction;
const Point<dim> second_intersection = origin + (-b - std::sqrt(discriminant))*direction;

// pick point closest to old point:
if( oldPoint.distance(first_intersection) < oldPoint.distance(second_intersection) )
 return first_intersection;
else
 return second_intersection;
}


/** these (v) should be members of the sphere class instead 
*/

   template<int dim>
   Tensor<1, dim> Geometry<dim>::getSphereNormalVector(const unsigned int sphere_id, 
                                               const Point<dim>& intersection_point) const
   {
     Tensor<1,dim> normal = intersection_point - spheres[sphere_id].getCenter();

     // rescale (normalize):
     normal /= normal.norm(); 

     return normal;
   }


   template<int dim>
   bool 
   Geometry<dim>::isInSphere(unsigned int sphere_id, 
                             const Point<dim>& location,
                             const double buffer) const
   {
     unsigned int number_spheres = spheres.size();

     if(sphere_id >= number_spheres)
       throw std::invalid_argument("Trying to access non-existing sphere.");

     if( location.distance(spheres[sphere_id].getCenter()) < 
           (spheres[sphere_id].getRadius() + buffer) )
       return true;

     return false;
   }

/** these (^) should be members of the sphere class instead 
*/

template<int dim>
bool Geometry<dim>::isInDomain(const Point<dim>& location) const
{
	// should have that point is in box, but still check:
	for(unsigned int dim_itr = 0; dim_itr < dim; dim_itr++)
	 if(location[dim_itr] < bottom_left[dim_itr] || location[dim_itr] > top_right[dim_itr])
	   return false;

	// check obstacles...
	unsigned int number_spheres = spheres.size();
	for(unsigned int sphere_id = 0; sphere_id < number_spheres; ++sphere_id)
	 if(isInSphere(sphere_id, location))
	   return false;

	unsigned int number_rectangles = rectangles.size();
	for(unsigned int rect = 0; rect < number_rectangles; ++rect)
	 if(rectangles[rect].distance_from_border(location) < 1e-8)
	   return false;

	return true;
}


template<int dim>
void 
Geometry<dim>::addPointBuffer(const double buffer, const Point<dim>& test_point,
               Point<dim>& buffered_point) const
{
buffered_point = test_point;

// add buffer to spheres-- assuming already outside of actual sphere:
unsigned int number_spheres = spheres.size();
for(unsigned int sphere_id = 0; sphere_id < number_spheres; ++sphere_id)
 if(isInSphere(sphere_id, test_point, buffer))
 {
   Tensor<1,dim> normal = getSphereNormalVector(sphere_id,
                                               test_point);
  
   buffered_point = buffered_point + buffer*normal;

   break;
 }

// add buffer to rectangles:
 unsigned int number_rectangles = rectangles.size();
 for(unsigned int rect_id = 0; rect_id < number_rectangles; ++rect_id)
   if( rectangles[rect_id].distance_from_border(buffered_point) < buffer )
   {
     // std::cout << "original point: " << buffered_point << std::endl;

     // add buffer to appropriate direction ...
     Tensor<1, dim> normal_vector = rectangles[rect_id].getNormalVector(buffered_point);
     buffered_point += (buffer*normal_vector);

     // std::cout << "\t buffered point: " << buffered_point << std::endl;
     break;
   } 
}


template<int dim>
std::vector<Point<dim> > Geometry<dim>::getQuerryPoints(double resolution) const
{
	if(dim != 2)
		throw std::invalid_argument("getQuerryPoints() not implemented for dim != 2");
	
	const unsigned int number_spheres = spheres.size();

	// if(number_spheres > 0 && dim != 2)
	//   throw std::invalid_argument("getQuerryPoints() not implemented for spheres in 3d");

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
		} // for circle points @todo : can add more points if querrying 3d spheres
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


template<int dim>
void Geometry<dim>::printInfo(std::ostream& out) const
{

out << std::endl << std::endl << Utility::medium_line << std::endl;
out << "\t\tGEOMETRY INFO:";
out << std::endl << Utility::medium_line << std::endl;

out << "DIMENSION: " << dim << std::endl
 << "\t BottomLeft: " << bottom_left << std::endl
 << "\t TopRight: " << top_right << std::endl;

out << "\nSCALE: " << scale << std::endl;

out << "\nBOUNDARY CONDITIONS: " << std::endl;
for(unsigned int i=0; i < dim; i++)
 out << "\t " << i << "th boundary: " 
 << getBoundaryConditionString(boundary_conditions[i]) << std::endl; // enum to ostream?

//@todo: enums to strings ... maybe c class/function that takes the enum and returns a string?

out << "\nSPHERES: " << spheres.size() << std::endl;
for(unsigned int i = 0; i < spheres.size(); i++)
 out << "\t Sphere Center: " << spheres[i].getCenter() 
   << " \tRadius: " << spheres[i].getRadius() << std::endl;

out << "\nRECTANGLES: " << rectangles.size() << std::endl;
for(unsigned int i = 0; i < rectangles.size(); i++)
 out << "\t Bottom left: " << rectangles[i].getBottomLeft() 
   << "\n\t Top right: " << rectangles[i].getTopRight() << std::endl;

out << "\nMESH TYPE: " << getMeshTypeString(mesh_type) << std::endl;
if(mesh_type == MeshType::FILE_MESH)
 out << "\t Mesh File: " << mesh_file << std::endl;

	out << Utility::medium_line << std::endl
		<< std::endl << std::endl;
}


template<int dim>
void 
Geometry<dim>::outputGeometry(std::string output_directory) const
{
// boundary:
std::ofstream boundary_out(output_directory + "/boundary.dat");
boundary_out << bottom_left << std::endl
<< top_right << std::endl;

// spheres:
std::ofstream spheres_out(output_directory + "/spheres.dat");
for(unsigned int i = 0; i < spheres.size(); ++i)
 spheres_out << spheres[i].getCenter() << " " << spheres[i].getRadius() << std::endl;

// rectangles:
std::ofstream rect_out(output_directory + "/rectangles.dat");
for(unsigned int i = 0; i < rectangles.size(); ++i)
 rect_out << rectangles[i].getBottomLeft() << " " 
     << rectangles[i].getTopRight() << std::endl;
}



}


#endif 
