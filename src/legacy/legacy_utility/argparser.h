#ifndef MICROBE_SIMULATOR_ARGPARSER2_H
#define MICROBE_SIMULATOR_ARGPARSER2_H

// for points:
#include <deal.II/base/point.h>

#include <ctime>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cstddef>
// #include <list>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

// to explore directory:
#include <dirent.h>
#include <sys/stat.h> 
#include "./enum_types.h"
#include "../geometry/geo_types.h"

/** @todo
* change ''sphere refinement'' to ''obstacle refinement''
* move some methods to parser tools?
*/

namespace MicrobeSimulator{

	// DATE-TIME FUNCTION:
	const std::string currentDateTime() 
	{
	    time_t     now = time(0);
	    struct tm  tstruct;
	    char       buf[80];
	    tstruct = *localtime(&now);
	    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	    // for more information about date/time format
	    strftime(buf, sizeof(buf), "%a-%m-%d-%y_%H-%M-%S", &tstruct);

	    return buf;
	} // for outputting data --- used in ArgParser

	bool stringToBool(std::string value)
	{
		if( (value.compare("TRUE") == 0) || (value.compare("1") == 0) )
			return true;
		else if( (value.compare("FALSE") == 0) || (value.compare("0") == 0) )
			return false;
		else
			throw std::invalid_argument("Invalid string! Could not convert string to boolean");
	}

	std::string boolToString(bool value)
	{
		if(value == true)
			return "True";
		else if(value == false)
			return "False";
		else
			throw std::invalid_argument("Invalid value! Could not convert bool to string");
	}

	// for parameters we're looping over:
	struct LoopedParameter{
		double start_value;
		double step_size;
		double end_value;
		std::string name;
	};


	class ArgParser{
	public:
		ArgParser(int argc, char** argv);

		// ACCESSORS:
		// geometry:
		std::string getGeometryFile() const;

		// advection:
		VelocityType getVelocityType() const;
		std::string getVelocityFile_X() const;
		std::string getVelocityFile_Y() const;
		double getMaximumVelocity() const;
    	double getVortexRadius() const;
    	double getVortexRotation() const;
    	unsigned int getStokesRefinement() const;

    	// common:
    	SimulatorType getSimulationType() const;
    	unsigned int getDimension() const;
    	double getTimeStep() const;
    	double getRunTime() const;
    	double getSavePeriod() const;
    	int getJobID() const;
    	std::string getOutputDirectory() const;
    	RunMode getRunMode() const;
    	unsigned int getNumberRunCycles() const;

    	// mesh:
    	unsigned int getGlobalRefinement() const;
    	unsigned int getSphereRefinement() const;
    	unsigned int getBoundaryRefinement() const;
    	std::string getMeshFile() const;
    	MeshType getMeshType() const;

    	// chemicals:
		unsigned int getNumberChemicals() const;
    	double getGoodDiffusionConstant() const;
    	double getWasteDiffusionConstant() const;
    	double getGoodDecayConstant() const;
    	double getWasteDecayConstant() const;
    	bool isSavingChemicals() const;

    	// bacteria:
    	unsigned int getNumberBacteria() const;
		unsigned int getNumberGroups() const;
		double getBacteriaDiffusionConstant() const;
		double getWasteSecretionRate() const;
		double getGoodSecretionRate() const;
		double getMutationRate() const;
		unsigned int getInitialNumberCheaters() const;
			// deterministic mutation:
		unsigned int getNumberMutate() const;
		double getDeterministicMutationTime() const;
		std::vector<dealii::Point<2> > getInitialLocations() const;
		bool isInitializeBacteriaFromField() const;
			// intialization subdomain:
		double getLeftSubdomainLength() const;

		// for reintroducing new groups:
		bool isAddingNewGroups() const;
		double getReintroductionPeriod() const;

		
		// fitness:
		double getAlphaGood() const;
		double getAlphaWaste() const;
		double getSecretionCost() const;
		double getGoodSaturation() const;
		double getWasteSaturation() const;

		// for specialization:
		double getInefficiencyPenalty() const;
		bool isUsingANDFitness() const;

		// DEBUGGING:
		bool isDebuggingCode() const;
		bool isPrintGrid() const;
		bool isPrintVelocity() const;
		bool isCheckMass() const;
		bool isPointSource() const;
		bool isInitialGaussian() const;
		bool isReproduceBacteria() const;
		double getFlowDelay() const;
		double getReproductionDelay() const;
		double getMutationDelay() const;

		// numerical experimentation:
		double getTimeStepFactor() const;
		double getViscosityBeta() const;

		// SOURCES:
		SourceImplementation getSourceImplementation() const;
		double getSourceResolution() const;

		// GEOMETRY:
		BoundaryCondition getXBoundaryBC() const;
		BoundaryCondition getYBoundaryBC() const;
		GeoTypes::Filter getFilter() const;
		GeoTypes::Mixer getMixer() const;
		GeoTypes::Pipe getPipe() const;
		GeoTypes::Splitter getSplitter() const;
		
		// output:
		void print(std::ostream& out) const;
		void outputParameters() const;

		// FOR INTERMITTENT CHEATING SIMULATIONS:
		unsigned int getNumberRegularBacteria() const;
		unsigned int getNumberICBacteria() const;
		unsigned int getNumberASBacteria() const;
		double getSpeciesMutationRate() const;
		double getICSwitchingRate() const;
		double getICOnPeriod() const;
		double getICOffPeriod() const;
	
		// FOR AGING:
		double getCellDensity() const;
		double getCoarseGrainResolution() const;

	private:
		void parse(int argc, char** argv);
		void checkParameters();
    	void parseDirectory(const char* directory);
    	// void assignParameters(const char* parameterFile); // legacy

    	void parseParameterFile(const char* parameterFile);
		void tokenizeLine(std::string line);

		void parseInitialLocations(std::string line);

		void assignParameter(std::string parameter, std::string value);
		void createLoopedParameter(std::string parameter, 
				std::string start, std::string step, std::string stop);
		void realizeLoopedParameters();
		std::vector<double> getGridParameterSet(std::vector<std::vector<double> > all_values, 
											unsigned int job_number,
											unsigned int total_grid_size) const;
		void assignGeometryMeshFiles();

		// geometry:
		std::string geometry_file;

		// looping over many geometry files:
		unsigned int number_geometry_files;
		std::string geometry_directory; // not used if files == 1 or lower

    	// advection:
    	VelocityType velocity_type;
    	std::string x_velocity_file;
    	std::string y_velocity_file;

    	double maximum_velocity;
    	double vortex_radius;
    	double vortex_rotation;
    	unsigned int stokes_refinement;

    	// common:
    	SimulatorType simulator_type;
    	unsigned int dimension;
    	double time_step;
    	double run_time;
    	double save_period;
    	int job_ID;
    	std::string output_directory;
    	RunMode run_mode;
    	unsigned int number_run_cycles;

    	// mesh:
    	unsigned int global_refinement;
    	unsigned int sphere_refinement;
    	unsigned int boundary_refinement;
    	std::string mesh_file;
    	 // for super simulator -- to swtich between simulation types
		MeshType mesh_type;

    	// chemicals:
    	unsigned int number_chemicals;
    	double good_diffusion_constant;
    	double waste_diffusion_constant;
    	double good_decay_constant;
    	double waste_decay_constant;
			// numerical experimentation:
		double time_factor_ck;
		double viscosity_beta;

    	// bacteria:
    	unsigned int number_bacteria;
    	unsigned int number_groups;
    	double bacteria_diffusion_constant;
    	double waste_secretion_rate;
    	double good_secretion_rate;
    	double mutation_rate;
    	unsigned int initial_cheaters;
    		// deterministic mutation:
    	unsigned int number_mutate;
    	double deterministic_mutation_time;
    	std::vector<dealii::Point<2> > initial_locations;
    	bool field_intialization;
    	double left_subdomain_length;
    	bool adding_new_groups;
    	double reintroduction_period;

    	// fitness:
		double alpha_good;
		double alpha_waste;
		double secretion_cost;
		double good_saturation;
		double waste_saturation;

		double inefficiency_penalty;
		bool use_and_fitness;

		// DEBUGGING:
		bool debug_code;
		bool print_grid;
		bool print_velocity;
		bool check_mass;
		bool point_source;
		bool initial_gaussian;
		bool reproduce_bacteria;
		bool save_chemicals;

		double flow_delay;
		double reproduction_delay;
		double mutation_delay;

		// compare source methods:
		SourceImplementation source_implementation;
		double source_resolution; 

		// FOR INTERMITTENT CHEATING:
		unsigned int number_reg_bacteria;
		unsigned int number_ic_bacteria;
		unsigned int number_as_bacteria;
		double species_mutation_rate;
		double ic_switching_rate;
		double ic_on_period;
		double ic_off_period;

		// For looping parameters:
		std::vector<LoopedParameter> looped_parameters;

		// FOR GEOMETRY:

		// boundaries:
		BoundaryCondition x_boundary;
		BoundaryCondition y_boundary;
		
		// filter:
		// GeoTypes::Filter filter; 
		double filter_height;
		double filter_number_walls; // can loop over -- this is an int!!
		double filter_wall_thickness;
		double filter_left_length;
		double filter_center_length;
		double filter_right_length;

		// mixer:
		// GeoTypes::Mixer mixer;
		double mixer_height;
		double mixer_left_length;
		double mixer_right_length;
		double mixer_radius;

		// variable 2d pipe ...
		double xmin;
		double xmax;
		double ymin;
		double ymax;

		// splitter:
		double splitter_height;
		double splitter_left;
		double splitter_right;
		double splitter_radius;

		// FOR AGING SIMULATIONS:
		// double aging_alpha;
		double cell_density;
		double coarse_grain_resolution;
		// also add hill constant k...

	};


// IMPLEMENTATION:
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------

// CONSTRUCTOR:
	ArgParser::ArgParser(int argc, char** argv)
		:
		number_geometry_files(0),
		velocity_type(VelocityType::NO_FLOW), maximum_velocity(0), 
		vortex_radius(0), vortex_rotation(0), stokes_refinement(3),
		job_ID(0), output_directory("./"), number_run_cycles(1),
		global_refinement(0), sphere_refinement(0), boundary_refinement(0), mesh_file(""),
		mutation_rate(0), number_mutate(0), initial_locations(), field_intialization(false),
		left_subdomain_length(-1),
		debug_code(false), print_grid(false), print_velocity(false), check_mass(false),
		point_source(false), initial_gaussian(false), reproduce_bacteria(false)
	{
		parse(argc,argv);
	}



// ACCESSORS:
//------------------------------------------------------------------------------------------
	// geometry:
	std::string ArgParser::getGeometryFile() const {return geometry_file;}

	// advection:
	VelocityType ArgParser::getVelocityType() const {return velocity_type;}
	std::string ArgParser::getVelocityFile_X() const {return x_velocity_file;}
	std::string ArgParser::getVelocityFile_Y() const {return y_velocity_file;}
	double ArgParser::getMaximumVelocity() const {return maximum_velocity;}
	double ArgParser::getVortexRadius() const {return vortex_radius;}
	double ArgParser::getVortexRotation() const {return vortex_rotation;}
	unsigned int ArgParser::getStokesRefinement() const {return stokes_refinement;}

	// common:
	SimulatorType ArgParser::getSimulationType() const {return simulator_type;}
	unsigned int ArgParser::getDimension() const {return dimension;}
	double ArgParser::getTimeStep() const {return time_step;}
	double ArgParser::getRunTime() const {return run_time;}
	double ArgParser::getSavePeriod() const {return save_period;}
	int ArgParser::getJobID() const {return job_ID;}
	std::string ArgParser::getOutputDirectory() const {return output_directory;}
	RunMode ArgParser::getRunMode() const {return run_mode;}
	unsigned int ArgParser::getNumberRunCycles() const {return number_run_cycles;}

	// mesh:
	unsigned int ArgParser::getGlobalRefinement() const {return global_refinement;}
	unsigned int ArgParser::getSphereRefinement() const {return sphere_refinement;}
	unsigned int ArgParser::getBoundaryRefinement() const {return boundary_refinement;}
	std::string ArgParser::getMeshFile() const {return mesh_file;}
	MeshType ArgParser::getMeshType() const {return mesh_type;}

	// chemicals:
	unsigned int ArgParser::getNumberChemicals() const {return number_chemicals;}
	double ArgParser::getGoodDiffusionConstant() const {return good_diffusion_constant;}
	double ArgParser::getWasteDiffusionConstant() const {return waste_diffusion_constant;}
	double ArgParser::getGoodDecayConstant() const {return good_decay_constant;}
	double ArgParser::getWasteDecayConstant() const {return waste_decay_constant;}
	bool ArgParser::isSavingChemicals() const {return save_chemicals;}

		// numerical experimentation:
	double ArgParser::getTimeStepFactor() const {return time_factor_ck;}
	double ArgParser::getViscosityBeta() const {return viscosity_beta;}

	// bacteria:
	unsigned int ArgParser::getNumberBacteria() const {return number_bacteria;}
	unsigned int ArgParser::getNumberGroups() const {return number_groups;}
	double ArgParser::getBacteriaDiffusionConstant() const {return bacteria_diffusion_constant;}
	double ArgParser::getWasteSecretionRate() const {return waste_secretion_rate;}
	double ArgParser::getGoodSecretionRate() const {return good_secretion_rate;}
	double ArgParser::getMutationRate() const {return mutation_rate;}
	unsigned int ArgParser::getInitialNumberCheaters() const {return initial_cheaters;}
		// deterministic mutation:
	unsigned int ArgParser::getNumberMutate() const {return number_mutate;}
	double ArgParser::getDeterministicMutationTime() const {return deterministic_mutation_time;}
	std::vector<dealii::Point<2> > ArgParser::getInitialLocations() const {return initial_locations;}
	bool ArgParser::isInitializeBacteriaFromField() const {return field_intialization;}
		// subdomain:
	double ArgParser::getLeftSubdomainLength() const {return left_subdomain_length;}

	// for adding in new groups:
	bool ArgParser::isAddingNewGroups() const {return adding_new_groups;}
	double ArgParser::getReintroductionPeriod() const {return reintroduction_period;}

	// fitness:
	double ArgParser::getAlphaGood() const {return alpha_good;}
	double ArgParser::getAlphaWaste() const {return alpha_waste;}
	double ArgParser::getSecretionCost() const {return secretion_cost;}
	double ArgParser::getGoodSaturation() const {return good_saturation;}
	double ArgParser::getWasteSaturation() const {return waste_saturation;}

	double ArgParser::getInefficiencyPenalty() const {return inefficiency_penalty;}
	bool ArgParser::isUsingANDFitness() const {return use_and_fitness;}

	//DEBUGGING:
	bool ArgParser::isDebuggingCode() const {return debug_code;}
	bool ArgParser::isPrintGrid() const {return print_grid;}
	bool ArgParser::isPrintVelocity() const {return print_velocity;}
	bool ArgParser::isCheckMass() const {return check_mass;}
	bool ArgParser::isPointSource() const {return point_source;}
	bool ArgParser::isInitialGaussian() const {return initial_gaussian;}
	bool ArgParser::isReproduceBacteria() const {return reproduce_bacteria;}
	double ArgParser::getFlowDelay() const {return flow_delay;}
	double ArgParser::getReproductionDelay() const {return reproduction_delay;}
	double ArgParser::getMutationDelay() const {return mutation_delay;}

	// SOURCES:
	SourceImplementation ArgParser::getSourceImplementation() const {return source_implementation;}
	double ArgParser::getSourceResolution() const {return source_resolution;}

	// FOR INTERMITTENT CHEATING SIMULATIONS:
	unsigned int ArgParser::getNumberRegularBacteria() const { return number_reg_bacteria;}
	unsigned int ArgParser::getNumberICBacteria() const {return number_ic_bacteria;}
	unsigned int ArgParser::getNumberASBacteria() const {return number_as_bacteria;}
	double ArgParser::getSpeciesMutationRate() const {return species_mutation_rate;}
	double ArgParser::getICSwitchingRate() const {return ic_switching_rate;}
	double ArgParser::getICOnPeriod() const {return ic_on_period;}
	double ArgParser::getICOffPeriod() const {return ic_off_period;}

	// FOR AGING:
	double ArgParser::getCellDensity() const {return cell_density;}
	double ArgParser::getCoarseGrainResolution() const {return coarse_grain_resolution;}

	// GEOMETRY:
	BoundaryCondition ArgParser::getXBoundaryBC() const {return x_boundary;}
	BoundaryCondition ArgParser::getYBoundaryBC() const {return y_boundary;}

	GeoTypes::Filter ArgParser::getFilter() const
	{
		GeoTypes::Filter filter;
		filter.number_channels = (unsigned int)filter_number_walls + 1;
		filter.wall_thickness = filter_wall_thickness;
		filter.left_length = filter_left_length;
		filter.center_length = filter_center_length;
		filter.right_length = filter_right_length;

		filter.channel_thickness = (filter_height - (filter_wall_thickness*filter_number_walls) ) 
			/ filter.number_channels;

		return filter;
	}

	GeoTypes::Mixer ArgParser::getMixer() const
	{
		GeoTypes::Mixer mixer;

		mixer.left_length = mixer_left_length;
		mixer.right_length = mixer_right_length;
		mixer.height = mixer_height;
		mixer.radius = mixer_radius;

		return mixer;
	}
		
	GeoTypes::Pipe ArgParser::getPipe() const
	{
		GeoTypes::Pipe pipe;

		pipe.xmin = xmin;
		pipe.xmax = xmax;
		pipe.ymin = ymin;
		pipe.ymax = ymax;

		return pipe;
	}

	GeoTypes::Splitter ArgParser::getSplitter() const
	{
		GeoTypes::Splitter splitter;

		splitter.height = splitter_height;
		splitter.left = splitter_left;
		splitter.right = splitter_right;
		splitter.radius = splitter_radius;

		return splitter;
	}

// PARSERS:
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
	void ArgParser::parse(int argc, char** argv)
	{

		std::cout << "\n\n...Parsing parameters\n" << std::endl;
		bool odSet = false;
		std::string dirTemp;

		if(argc == 1)
		{
			std::cout << "...No input files given. Using *UNITIIALIZED* parameters. \n"; // won't work...
		}
		else if((argc % 2) == 1) // program + flags + values == odd
		{
			for(int i = 1; i < argc; i += 2)
			{
			  std::string flag = argv[i]; 

			  if(flag.compare("-d") == 0) { parseDirectory(argv[i+1]); }
			  else if(flag.compare("-fp") == 0) { parseParameterFile(argv[i+1]); }
			  else if(flag.compare("-fg") == 0) { geometry_file = argv[i+1]; }
  			  else if(flag.compare("-fm") == 0) { mesh_file = argv[i+1]; }
			  else if(flag.compare("-fvx") == 0) { x_velocity_file = argv[i+1]; }
			  else if(flag.compare("-fvy") == 0) { y_velocity_file = argv[i+1]; }
			  else if(flag.compare("-o") == 0) { odSet = true; dirTemp = argv[i+1]; }
			  else if(flag.compare("-id") == 0 ) { job_ID = atoi(argv[i+1]); }
			  // over write parameters:
			  else if(flag.compare("-diff1") == 0) { good_diffusion_constant = atof(argv[i+1]); }
			  else if(flag.compare("-diff2") == 0) { waste_diffusion_constant = atof(argv[i+1]); }
			  else if(flag.compare("-vscale") == 0) { maximum_velocity = atof(argv[i+1]); }
			  else if(flag.compare("-mrate") == 0) { mutation_rate = atof(argv[i+1]); }
			  else if(flag.compare("-gref") == 0) { global_refinement = (unsigned int)atof(argv[i+1]); }
			  else if(flag.compare("-tmax") == 0) { run_time = atof(argv[i+1]); }
			  else
			  {
			    std::ostringstream message;
			    message << std::endl
			      << "***UNKNOWN FLAG*** \n"
			      << "Usage: ./program [-flags] (file or directory) \n"
			      << "\t -d \"directory\" \n"
			      << "\t -fp \"parameter file\" \n"
			      << "\t -fg \"geometry file\" \n"
			      << "\t -fvx \"x-velocity file\" \n"
			      << "\t -fvy \"y-velocity file\" \n" 
			      << "\t -o \"output directory\" \n"
			      << "\t -id \"job id\" \n" << std::endl;
			    throw std::invalid_argument(message.str());  // SHOULD THROW AN EXECPTION HERE
			  }

			} // for all inputs
		} // if even inputs (-flag value)
		else 
		{
			std::ostringstream message;
			message << std::endl
			  << "***INVALID NUMBER OF INPUT PARAMETERS*** \n"
			  << "Usage: ./program [-flags] (file or directory) \n"
			  << "\t -d \"directory\" \n"
			  << "\t -fp \"parameter file\" \n"
			  << "\t -fg \"geometry file\" \n"
			  << "\t -fvx \"x-velocity file\" \n"
			  << "\t -fvy \"y-velocity file\" \n" 
			  << "\t -o \"output directory\" \n"
			  << "\t -id \"job id\" \n" << std::endl;
			throw std::invalid_argument(message.str()); 
		} // if correct number of inputs

		 // CREATE OUTPUT DIRECTORY:
		if(odSet)
			output_directory =  "./data_" + dirTemp + "_jobID_" 
				+  std::to_string(job_ID) + "_" + currentDateTime();
		else
			output_directory =  "./data_jobID_" +  std::to_string(job_ID) + "_" + currentDateTime();

		int check;
		check = mkdir(output_directory.c_str(), 0777);

		if(check == 0)
			std::cout << "...Data output directory created: " << output_directory << std::endl;
		else
			throw std::runtime_error("\n*** FAILED TO CREATE OUTPUT DIRECTORY ***\n"); 

		// note order important here -- realizedLoopedParameters() 
			// outputs text file with grid parameters into output_directory
		// CREATE LOOPED PARAMETER GRID AND REALIZATION IF NEEDED:
		realizeLoopedParameters();

		// assign (possibly overwrite) geometry file if using many files:
		if(number_geometry_files > 1)
			assignGeometryMeshFiles();

		// CHECK PARAMETERS:
		checkParameters();
	}


	void ArgParser::parseDirectory(const char* directory)
	{
	  std::string direc(directory);
	  std::string::iterator it = direc.end();
	  --it;
	  if( (*it) != '/')
	    direc += "/"; // add the slash

	  std::cout << "...Configuration directory: " << direc << std::endl;

	  DIR *dir; // pointer to open directory
	  struct dirent *entry; // stuff in directory
	  // struct stat info; // information about each entry

	  //1 open
	  dir = opendir(directory);

	  if(!dir)
  	    throw std::runtime_error("\n***INPUT PARAMETER DIRECTORY NOT FOUND***\n");

	  // 2 read:
	  while( (entry = readdir(dir)) != NULL )
	  {
	    if(entry->d_name[0] != '.')
	    {
	      std::string fileName = entry->d_name;

	      std::size_t found = fileName.find_last_of(".");
	      std::string extention = fileName.substr(found+1);

	      if(extention.compare("msh") == 0){ mesh_file = direc + fileName; }
	      if(extention.compare("dat") == 0)
	      {
	        std::size_t fnd = fileName.find_last_of("_");
	        std::string type = fileName.substr(fnd+1);

	        if(type.compare("velx.dat") == 0){ x_velocity_file = direc + fileName; }
	        if(type.compare("vely.dat") == 0){ y_velocity_file = direc + fileName; }
	        if(type.compare("parameters.dat") == 0 || type.compare("para.dat") == 0)
	        {
	          std::string pfile = direc + fileName;
	          parseParameterFile(pfile.c_str());
	        }
	        if(type.compare("geo.dat") == 0)
	        {
	          geometry_file = direc + fileName;
	        }

	      }

	    } // if not current or parent directory -- should generalize to is not directory... 
	  } // while files in directory

	  // 3 close:
	  closedir(dir);

	} // parseDirectory()



	void ArgParser::parseParameterFile(const char* parameterFile)
	{
		std::cout << "...Reading parameter file: " << parameterFile << std::endl;

		std::ifstream infile(parameterFile);
		std::string line;

		// read in each line:
		while(std::getline(infile,line)){
			// tokenize each line:
			tokenizeLine(line); // assigns parameters and looped parameters
		}
	}


	void ArgParser::tokenizeLine(std::string line)
	{
		std::vector<std::string> tokens;
		std::string trimmed_line = line; 
		boost::trim(trimmed_line);

		// ignore comments:
		if(trimmed_line.size() == 0)
			return;
		else if(trimmed_line.size() >= 2)
			if(trimmed_line[0]=='/' && trimmed_line[1]=='/')
				return;

		// tokenize line:
		boost::split(tokens,trimmed_line,boost::is_any_of(" ,:="),boost::token_compress_on);

		const unsigned int number_tokens = tokens.size();

		// handle initial locations separately:
		if(number_tokens != 0)
		{
			if(tokens[0].compare("initial_locations") == 0)
			{
				parseInitialLocations(trimmed_line);
				return;
			}
			else if(tokens[0].compare("geometry_files") == 0)
			{
				number_geometry_files = std::stoul(tokens[1]);
				geometry_directory = tokens[2];
				return;
			} // if using many geometries
		}

		// assign parameter or create looped parameter:
		if(number_tokens == 2)
		{
			assignParameter(tokens[0], tokens[1]);
		}
		else if(number_tokens == 4)
		{
			createLoopedParameter(tokens[0], tokens[1], tokens[2], tokens[3]);
		}
		else
		{
			std::ostringstream message;
			message << "Invalid parameter entry: \n\t<" << line << ">\n ";
			throw std::invalid_argument(message.str());
		}
	}

	
	void ArgParser::parseInitialLocations(std::string line)
	{
		/** @todo 
		* NOTE THIS IS CURRENTLY IMPLEMENTED FOR 2D POINTS!!!
		* NEED TO GENERALIZE TO 3D STILL
		*/

		// std::cout << "... parsing intial bacteria locations\n\n " << std::endl;

		// tokenize line:
		std::vector<std::string> tokens;
		boost::split(tokens,line,boost::is_any_of(" ,:={}()"),boost::token_compress_on);

		// gives empty token at end.... check for string length
		// first token corresponds to ``intial location'' parameter
		for(unsigned int i = 2; i < tokens.size(); i += 2)
			if(tokens[i].size() != 0)
				initial_locations.push_back(dealii::Point<2>( std::stod(tokens[i-1]), 
							std::stod(tokens[i]) ) );
	}


	void ArgParser::assignGeometryMeshFiles()
	{
		if(job_ID < 1)
			throw std::invalid_argument("job id must be at least 1 for using looped geometry files");

		const unsigned int geo_ID = ((job_ID - 1) % number_geometry_files) + 1;

		// look through geometry directory:

		DIR *dir; // pointer to open directory
		struct dirent *entry; // stuff in directory

		// open:
		dir = opendir(geometry_directory.c_str());

		if(!dir)
  	    	throw std::runtime_error("\n***GEOMETRY DIRECTORY NOT FOUND***\n");

  	    std::cout << "Reading geometry data from " << geometry_directory << std::endl;
  	    // read:
  	    while( (entry = readdir(dir)) != NULL )
  	    {
  	    	if(entry->d_name[0] != '.')
  	    	{
  	    		std::string fileName = entry->d_name;
  	    		std::string tag = "ID_" + std::to_string(geo_ID) + "_";

  	    		std::size_t found = fileName.find(tag);

  	    		if (found!=std::string::npos)
  	    		{
	  	    		std::string extention = fileName.substr(found+tag.length());
	  	    		// std::cout << "file prefix: " << extention << std::endl;

	  	    		// get mesh and geometry files from extension:
	  	    		std::size_t fnd = extention.find_last_of(".");
	  	    		std::string file_exten = extention.substr(fnd+1);

	  	    		if(file_exten.compare("msh") == 0){mesh_file = geometry_directory 
	  	    				+ tag + extention;}
  	    			else if(file_exten.compare("dat") == 0){geometry_file = geometry_directory
  	    					+ tag + extention;}
	  	    	}

  	    	}
  	    }
	}


	void ArgParser::assignParameter(std::string parameter, std::string value)
	{
		//CASES:

		// common:
		if(parameter.compare("simulator_type") == 0){simulator_type = stringToSimulatorType(value);}
		else if(parameter.compare("dimension") == 0){dimension = std::stoul(value);}
		else if(parameter.compare("run_mode") == 0){run_mode = stringToRunMode(value);} // not used

	    else if(parameter.compare("run_time") == 0){run_time = std::stod(value);}
	    else if(parameter.compare("time_step") == 0){time_step = std::stod(value);} // not really used ...
	    else if(parameter.compare("save_period") == 0){save_period = std::stod(value);}       
	    else if(parameter.compare("number_run_cycles") == 0){number_run_cycles = std::stoul(value);}

	    //mesh:
	    else if(parameter.compare("global_refinement") == 0){global_refinement = std::stoul(value);}
	    else if(parameter.compare("sphere_refinement") == 0){sphere_refinement = std::stoul(value);}
	    else if(parameter.compare("boundary_refinement") == 0){boundary_refinement = std::stoul(value);}
	    else if(parameter.compare("mesh_type") == 0){mesh_type = stringToMeshType(value);}

	    //bacteria:
	    else if(parameter.compare("number_bacteria") == 0){number_bacteria = std::stoul(value);} 
	    else if(parameter.compare("groups") == 0){number_groups = std::stoul(value);} 
	    else if(parameter.compare("bacteria_diffusion") == 0)
	    	{bacteria_diffusion_constant = std::stod(value);}
	    else if(parameter.compare("good_secretion") == 0){good_secretion_rate = std::stod(value);}
	    else if(parameter.compare("waste_secretion") == 0){waste_secretion_rate = std::stod(value);}
	    else if(parameter.compare("mutation_rate") == 0){mutation_rate = std::stod(value);}
	    else if(parameter.compare("initial_cheaters") == 0){initial_cheaters = std::stoul(value);}
	    	// deterministic mutation:
	    else if(parameter.compare("number_mutate") == 0){number_mutate = std::stoul(value);}
	    else if(parameter.compare("deterministic_mutation_time") == 0)
	    	{deterministic_mutation_time = std::stod(value);}
	    else if(parameter.compare("field_intialization") == 0){field_intialization = stringToBool(value);}
	    else if(parameter.compare("left_subdomain_length") == 0){left_subdomain_length = std::stod(value);}
	  	
		  	// for adding new groups:
	    else if(parameter.compare("adding_new_groups") == 0){adding_new_groups = stringToBool(value);}
	    else if(parameter.compare("reintroduction_period") == 0){reintroduction_period = std::stod(value);}
	    
	    //fitness:
		else if(parameter.compare("alpha_good") == 0){alpha_good = std::stod(value);}
	    else if(parameter.compare("alpha_waste") == 0){alpha_waste = std::stod(value);}
	    else if(parameter.compare("good_saturation") == 0){good_saturation = std::stod(value);}
	    else if(parameter.compare("waste_saturation") == 0){waste_saturation = std::stod(value);}
	    else if(parameter.compare("secretion_cost") == 0){secretion_cost = std::stod(value);}
	    else if(parameter.compare("inefficiency_penalty") == 0){inefficiency_penalty = std::stod(value);}
	    else if(parameter.compare("use_and_fitness") == 0){use_and_fitness = stringToBool(value);}

	    // chemicals:
	    else if(parameter.compare("number_chemicals") == 0){number_chemicals = std::stoul(value);}
	    else if(parameter.compare("good_diffusion") == 0){good_diffusion_constant = std::stod(value);}
	    else if(parameter.compare("waste_diffusion") == 0){waste_diffusion_constant = std::stod(value);}
	    else if(parameter.compare("good_decay") == 0){good_decay_constant = std::stod(value);}
	    else if(parameter.compare("waste_decay") == 0){waste_decay_constant = std::stod(value);}
	    else if(parameter.compare("save_chemicals") == 0){save_chemicals = stringToBool(value);}

    		// numerical experimentation:
		else if(parameter.compare("time_factor_ck") == 0){time_factor_ck = std::stod(value);}
		else if(parameter.compare("viscosity_beta") == 0){viscosity_beta = std::stod(value);}


	    // advection:
	    else if(parameter.compare("velocity_type") == 0 ){velocity_type = stringToVelocityType(value);} 
		else if(parameter.compare("maximum_velocity") == 0){maximum_velocity = std::stod(value);}
		else if(parameter.compare("vortex_radius") == 0){vortex_radius = std::stod(value);}
		else if(parameter.compare("vortex_rotation") == 0){vortex_rotation = std::stod(value);}
		else if(parameter.compare("stokes_refinement") == 0){stokes_refinement = std::stoul(value);}
		// sources:
		else if(parameter.compare("source_implementation") == 0)
			{source_implementation = stringToSourceImplementation(value);}
		else if(parameter.compare("source_resolution") == 0){source_resolution = std::stod(value);}

		//DEBUGGING:
		else if(parameter.compare("debug_mode") == 0){debug_code = stringToBool(value);}
		else if(parameter.compare("print_grid") == 0){print_grid = stringToBool(value);}
		else if(parameter.compare("print_velocity") == 0){print_velocity = stringToBool(value);}
		else if(parameter.compare("check_mass") == 0){check_mass = stringToBool(value);}
		else if(parameter.compare("point_source") == 0){point_source = stringToBool(value);}
		else if(parameter.compare("initial_gaussian") == 0){initial_gaussian = stringToBool(value);}
		else if(parameter.compare("reproduce_bacteria") == 0){reproduce_bacteria = stringToBool(value);}
		else if(parameter.compare("flow_delay") == 0){flow_delay = std::stod(value);}
		else if(parameter.compare("reproduction_delay") == 0){reproduction_delay = std::stod(value);}
		else if(parameter.compare("mutation_delay") == 0){mutation_delay = std::stod(value);}

		// FOR INTERMITTENT CHEATING:
		else if(parameter.compare("number_reg_bacteria") == 0){number_reg_bacteria = std::stoul(value);}
		else if(parameter.compare("number_ic_bacteria") == 0){number_ic_bacteria = std::stoul(value);}
		else if(parameter.compare("number_as_bacteria") == 0){number_as_bacteria = std::stoul(value);}
		else if(parameter.compare("species_mutation_rate") == 0){species_mutation_rate = std::stod(value);}
		else if(parameter.compare("ic_switching_rate") == 0){ic_switching_rate = std::stod(value);}
		else if(parameter.compare("ic_on_period") == 0){ic_on_period = std::stod(value);}
		else if(parameter.compare("ic_off_period") == 0){ic_off_period = std::stod(value);}
		
		// GEOMETRY:

		// boundaries:
		else if(parameter.compare("x_boundary") == 0){x_boundary = stringToBoundaryCondition(value);}
		else if(parameter.compare("y_boundary") == 0){y_boundary = stringToBoundaryCondition(value);}

		// filter:
		else if(parameter.compare("filter_height") == 0){filter_height = std::stod(value);}
		else if(parameter.compare("filter_number_walls") == 0){filter_number_walls = std::stod(value);}
		else if(parameter.compare("filter_wall_thickness") == 0){filter_wall_thickness = std::stod(value);}
		else if(parameter.compare("filter_left_length") == 0){filter_left_length = std::stod(value);}
		else if(parameter.compare("filter_center_length") == 0){filter_center_length = std::stod(value);}
		else if(parameter.compare("filter_right_length") == 0){filter_right_length = std::stod(value);}

		//  mixer:
		else if(parameter.compare("mixer_height") == 0){mixer_height = std::stod(value);}
		else if(parameter.compare("mixer_left_length") == 0){mixer_left_length = std::stod(value);}
		else if(parameter.compare("mixer_right_length") == 0){mixer_right_length = std::stod(value);}
		else if(parameter.compare("mixer_radius") == 0){mixer_radius = std::stod(value);}

		// pipe:
		else if(parameter.compare("xmin") == 0){xmin = std::stod(value);}
		else if(parameter.compare("xmax") == 0){xmax = std::stod(value);}
		else if(parameter.compare("ymin") == 0){ymin = std::stod(value);}
		else if(parameter.compare("ymax") == 0){ymax = std::stod(value);}

		// splitter:
		else if(parameter.compare("splitter_height") == 0){splitter_height = std::stod(value);}
		else if(parameter.compare("splitter_left") == 0){splitter_left = std::stod(value);}
		else if(parameter.compare("splitter_right") == 0){splitter_right = std::stod(value);}
		else if(parameter.compare("splitter_radius") == 0){splitter_radius = std::stod(value);}
		
		// FOR AGING SIMULATIONS:
		else if(parameter.compare("coarse_grain_resolution")==0){coarse_grain_resolution = std::stod(value);}
		else if(parameter.compare("cell_density")==0){cell_density = std::stod(value);}

		// ELSE: THROW EXCEPTION:
		else{
			std::string message = "INVALID PARAMETER ARUGMENT: <" + parameter + ">\n";
			throw std::invalid_argument(message.c_str());
		}
	}


	void ArgParser::createLoopedParameter(std::string parameter, 
			std::string start, std::string step, std::string stop)
	{
		std::cout << "...adding looped parameter: " << parameter << std::endl;

		LoopedParameter para;
		para.start_value = std::stod(start);
		para.step_size = std::stod(step);
		para.end_value = std::stod(stop);
		para.name = parameter;

		looped_parameters.push_back(para);

		std::cout << "\t>> start: " << para.start_value << " | step: " << para.step_size
			<< " | stop: " << para.end_value << std::endl;
	}


	// create looped parameter grid and realization if needed:
	void ArgParser::realizeLoopedParameters()
	{
		const unsigned int number_looping = looped_parameters.size();
		if( number_looping == 0 )
			return;

		if( job_ID < 1 )
			throw std::invalid_argument("Job id must be greater than or equal to 1 to use looped parameters");

		std::ostringstream grid;
		grid << "-----------------------------------------------------" << std::endl
			<< "\t\t PARAMETER GRID:" << std::endl
			<< "-----------------------------------------------------\n" << std::endl;

		grid << "Job ID \t";

		std::vector<std::vector<double> > all_values;
		// create grid ...
		// job-id p1 p2 ... pn this(*)
		unsigned int total_grid_size = 1;

		for(unsigned int i = 0; i < number_looping; ++i)
		{
			grid << looped_parameters[i].name << "\t";
			std::vector<double> values;
			double current = looped_parameters[i].start_value;
			do
			{
				values.push_back(current);
				current += looped_parameters[i].step_size;
			} while (current <= looped_parameters[i].end_value);

			total_grid_size *= values.size();

			all_values.push_back(values);
		}
		grid << "this" << std::endl;

		unsigned int job_number = 1;
		// grid values:
		for(unsigned int i = 0; i < total_grid_size; ++i)
		{
			grid << job_number << "\t\t";

			std::vector<double> parameter_set = getGridParameterSet(all_values, 
					job_number, total_grid_size);

			for(unsigned int j = 0; j < parameter_set.size(); ++j)
				grid << parameter_set[j] << " \t\t";

			if( (((job_ID - 1) % total_grid_size) + 1) == job_number )
			{
				grid << "X";
				for(unsigned int k = 0; k < parameter_set.size(); ++k)
					assignParameter(looped_parameters[k].name, 
						std::to_string(parameter_set[k])); // works for now
			}
			grid << std::endl;

			++job_number;
		} // for all possible grid combinations

		grid << "\nTotal grid size is: " << total_grid_size << std::endl 
			<< "-----------------------------------------------------\n" << std::endl;

		// output grid info to screen and file
		std::cout << std::endl << std::endl;
		std::cout << grid.str(); 

		// output to file:
		std::string outFile = output_directory + "/parameterGrid.dat";
		std::ofstream out(outFile);
		out << grid.str();
	}


	std::vector<double> ArgParser::getGridParameterSet(std::vector<std::vector<double> >all_values, 
														unsigned int job_number,
														unsigned int total_grid_size) const
	{
		std::vector<double> these_values;

		const unsigned int this_job = ((job_number -1) % total_grid_size);

		// find index along each parameter:
		const unsigned int number_looped = all_values.size();
		std::vector<unsigned int> all_index(number_looped,0); // intial multi-index

		unsigned int divide_factor = 1;
		for(unsigned int i = 0; i < number_looped; ++i)
		{
			const unsigned int current_index = ((this_job/divide_factor) % all_values[i].size());

			all_index[i] += current_index;
			
			// increase divide factor for next parameter set:
			divide_factor *= all_values[i].size(); 

			these_values.push_back(all_values[i][all_index[i]]);
		}

		return these_values;
	}


	void ArgParser::checkParameters()
	{
	  if(bacteria_diffusion_constant < 0.0 || good_diffusion_constant < 0.0 || 
	  		waste_diffusion_constant < 0.0)
	    throw std::invalid_argument("ERROR: MUST CHOOSE POSITIVE DIFFUSION CONSTANTS.");
	  if(good_decay_constant < 0.0 || waste_decay_constant < 0.0)
	    throw std::invalid_argument("ERROR: MUST CHOOSE POSITIVE DECAY CONSTANTS");
	  if(good_secretion_rate < 0.0 || waste_secretion_rate < 0.0)
	    throw std::invalid_argument("ERROR: MUST CHOOSE POSITIVE SECRETION CONSTANTS");
	  if(alpha_good < 0.0 || alpha_waste < 0.0 || secretion_cost < 0.0)
	    throw std::invalid_argument("ERROR: MUST CHOOSE POSITIVE FITNESS CONSTANTS");
	  if(good_saturation < 0.0 || waste_saturation < 0.0)
	    throw std::invalid_argument("ERROR: MUST CHOOSE POSITIVE SATURATION CONSTANTS");
	  if(time_step < 0 || run_time < 0)
	    throw std::invalid_argument("ERROR: TIME CONSTANTS MUST BE POSITIVE"); // save rate is unsigned
	  // if(Nb < 0)
	  //   throw std::invalid_argument("ERROR: MUST CHOOSE POSITIVE NUMBER OF INITIAL BACTERIA"); // Nb is unsigned...
	} // checkParameters()


// OUTPUT:
// ----------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------

	void ArgParser::print(std::ostream& out) const
	{
		std::ostringstream programData;

		programData<< "\n\n-----------------------------------------------------" << std::endl
			<< "\t\tSYSTEM PARAMETERS:" << std::endl
			<< "-----------------------------------------------------" << std::endl
			<< "\nSIMULATOR: " << getSimulatorTypeString(simulator_type) << std::endl 
		 	<< "\nDATE_TIME: " << currentDateTime() << std::endl
		 	<< "\nINPUT FILES: " << std::endl
		 	<< "\t Geometry File: " << geometry_file << std::endl
		 	<< "\t X-Velocity FIle: " << x_velocity_file << std::endl
		 	<< "\t Y-Velocity File: " << y_velocity_file << std::endl
		 	<< "\nVELOCITY:" << std::endl
		 	<< "\t Velocity Type: " << getVelocityTypeString(velocity_type) << std::endl
		 	<< "\t Maximum velocity: " << maximum_velocity << std::endl
		 	<< "\t Vortex radius: " << vortex_radius << std::endl
		 	<< "\t Vortex rotation: " << vortex_rotation << std::endl
		 	<< "\t Stokes refinement: " << stokes_refinement << std::endl
		 	<< "\nCHEMICALS: " << std::endl
		 	<< "\t Number chemicals: " << number_chemicals << std::endl
		 	<< "\t Good diffusion: " << good_diffusion_constant << std::endl
		 	<< "\t Waste diffusion: " << waste_diffusion_constant << std::endl
		 	<< "\t Good decay rate: " << good_decay_constant << std::endl
		 	<< "\t Waste decay rate: " << waste_decay_constant << std::endl	
		 	<< "\t Saving chemicals: " << boolToString(save_chemicals) << std::endl
		 	<< "\nBACTERIA: " << std::endl
		 	<< "\t Number initial bacteria: " << number_bacteria << std::endl
		 	<< "\t Number initial groups: " << number_groups << std::endl
		 	<< "\t Bacteria diffusion: " << bacteria_diffusion_constant << std::endl
		 	<< "\t Good secretion rate: " << good_secretion_rate << std::endl
		 	<< "\t Waste secretion rate: " << waste_secretion_rate << std::endl
		 	<< "\t Initial number cheaters: " << initial_cheaters << std::endl
		 	<< "\t Mutation rate: " << mutation_rate << std::endl
		 	<< "\t Number deterministic mutate: " << number_mutate << std::endl
		 	<< "\t Deterministic mutation time: " << deterministic_mutation_time << std::endl
		 	<< "\t Intialize from field: " << boolToString(field_intialization) << std::endl
		 	<< "\t Left subdomain length: " << left_subdomain_length << std::endl;

		 if( !initial_locations.empty() && !field_intialization )
		 {
		 	programData << "\t Intial locations:" << std::endl;
		 	for(unsigned int i = 0; i < initial_locations.size(); ++i )
		 		programData << "\t\t" << initial_locations[i] << std::endl;
		 }

		 programData << "\nFITNESS: " << std::endl
		 	<< "\t Alpha good: " << alpha_good << std::endl
		 	<< "\t Alpha waste: " << alpha_waste << std::endl
		 	<< "\t Secretion cost: " << secretion_cost << std::endl
		 	<< "\t Good saturation: " << good_saturation << std::endl
		 	<< "\t Waste saturation: " << waste_saturation << std::endl
		 	<< "\t Inefficiency penalty: " << inefficiency_penalty << std::endl
		 	<< "\nMESH REFINEMENT:" << std::endl
		 	<< "\t Global refinement: " << global_refinement << std::endl
		 	<< "\t Sphere refinement: " << sphere_refinement << std::endl
		 	<< "\t Boundary refinement: " << boundary_refinement << std::endl
		 	<< "\t Mesh file: " << mesh_file << std::endl
		 	<< "\t Mesh type: " << getMeshTypeString(mesh_type) << std::endl;

		 if(mesh_type == MeshType::FILTER)
		 {
		 	GeoTypes::Filter filter = getFilter();
		 	filter.printInfo(programData);
		 }
		 else if(mesh_type == MeshType::MIXER)
		 {
		 	GeoTypes::Mixer mixer = getMixer();
		 	mixer.printInfo(programData);
		 }
		 else if(mesh_type == MeshType::BOX_MESH)
		 {
		 	GeoTypes::Pipe pipe = getPipe();
		 	pipe.printInfo(programData);
		 }
		 else if(mesh_type == MeshType::SPLITTER)
		 {
		 	GeoTypes::Splitter splitter = getSplitter();
		 	splitter.printInfo(programData);
		 }

		 if(simulator_type == SimulatorType::AGING)
		 	programData << "\nAGING: " << std::endl
				<< "\t Coarse grain resolution: " << coarse_grain_resolution << std::endl
				<< "\t Cell density: " << cell_density << std::endl;
		 	
		 programData << "\nRUNNING/SAVING: " << std::endl		//@todo print run mode and if debugging
		 	<< "\t Dimension: " << dimension << std::endl
		 	<< "\t Run mode: " << getRunModeString(run_mode) << std::endl
		 	<< "\t Time step: " << time_step << std::endl
		 	<< "\t Run time: " << run_time << std::endl
		 	<< "\t Save period: " << save_period << std::endl
		 	<< "\t Job ID: " << job_ID << std::endl
		 	<< "\t Output directory: " << output_directory << std::endl
		 	<< "\t Number cycles: " << number_run_cycles << std::endl
		 	<< "\n\n USING SOURCE IMPLEMENTATION: " 
		 		<< getSourceImplementationString(source_implementation) << std::endl
	 		<< "\t Source resolution: " << source_resolution << std::endl;
		 	
		if(debug_code)
		{
			programData << "\n***RUNNING DEBUG PARAMETERS***" << std::endl
				<< "\t Print grid: " << boolToString(print_grid) << std::endl
				<< "\t Print velocity: " << boolToString(print_velocity) << std::endl
				<< "\t Check mass: " << boolToString(check_mass) << std::endl
				<< "\t Point source: " << boolToString(point_source) << std::endl
				<< "\t Initial gaussian: " <<  boolToString(initial_gaussian) << std::endl
				<< "\t Reproduce bacteria: " << boolToString(reproduce_bacteria) << std::endl
				<< "\t Flow delay: " << flow_delay << std::endl
				<< "\t Reproduction delay: " << reproduction_delay << std::endl
				<< "\t Mutation delay: " << mutation_delay << std::endl
				<< "\t Time factor ck:" << time_factor_ck << std::endl
				<< "\t Viscosity beta: " << viscosity_beta << std::endl;
		}

		programData << "-----------------------------------------------------\n\n\n" << std::endl;


		out << programData.str(); 
	}

	void ArgParser::outputParameters() const
	{
		  std::string outFileName = output_directory + "/programData.dat";
		  std::ofstream outFile(outFileName);
		  print(outFile);
	}


}  // namespace MicrobeSimulator
#endif 