#ifndef CHEMOTAXIS_ARGPARSER_BASE_H	
#define CHEMOTAXIS_ARGPARSER_BASE_H	

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

#include "./utility.h"

namespace Chemotaxis{ 

class ArgParserBase{
public:
	ArgParserBase(int argc, char** argv);
	virtual ~ArgParserBase() {}

	// OUTPUT:
	void print(std::ostream& out) const;
	void outputParameters() const;

	// ACCESSORS: // included in derived classes

protected:
	void parse(int argc, char** argv); // done
	virtual void assignParameter(std::string parameter, std:string value);

private:
	// METHODS:
	void checkParameters(); 
	void parseDirectory(const char* directory); // done

	void parseParameterFile(const char* parameterFile); // done
	void tokenizeLine(std::string line); // done

	void parseInitialLocations(std::string line); // done

	void createLoopedParameter(std::string parameter,
		std::string start, std::string step, std::string stop);
	void realizeLoopedParameters();
	std::vector<double> getGridParameterSet(std::vector<std::vector<double> > all_values,
											unsigned int job_number,
											unsigned int total_grid_size) const;

	// For looping parameters:
	std::vector<LoopedParameter> looped_parameters;

	// PARAMETERS:

};

// IMPL
// ------------------------------------------------------------------------------------------------

// CONSTRUCTOR:
ArgParserBase::ArgParserBase(int argc, char** argv)
	// :
{
	parse(argc, argv);
}

// ACCESSORS:
// -----------------------------------------------------------------------------------------


// PARSER METHODS:
// -----------------------------------------------------------------------------------------
void ArgParserBase::parse(int argc, char** argv)
{
	std::cout << "\n\n...Parsing parameters\n" << std::endl;
	bool odSet = false;
	std::string dirTemp;

	if(argc == 1)
	{
		std::cout << "...No input files given. Using *UNITIIALIZED* parameters. \n"; // won't work... // could run tests...
	}
	else if((argc % 2) == 1) // program + flags + values == odd
	{
		for(int i = 1; i < argc; i += 2)
		{
		  std::string flag = argv[i]; 

		  if(flag.compare("-d") == 0) { parseDirectory(argv[i+1]); }
		  else if(flag.compare("-fp") == 0) { parseParameterFile(argv[i+1]); }
		  // else if(flag.compare("-fg") == 0) { geometry_file = argv[i+1]; }
			  // else if(flag.compare("-fm") == 0) { mesh_file = argv[i+1]; }
		  // else if(flag.compare("-fvx") == 0) { x_velocity_file = argv[i+1]; }
		  // else if(flag.compare("-fvy") == 0) { y_velocity_file = argv[i+1]; }
		  else if(flag.compare("-o") == 0) { odSet = true; dirTemp = argv[i+1]; }
		  else if(flag.compare("-id") == 0 ) { job_ID = atoi(argv[i+1]); }
		  // over write parameters:
		  // else if(flag.compare("-diff1") == 0) { good_diffusion_constant = atof(argv[i+1]); }
		  // else if(flag.compare("-diff2") == 0) { waste_diffusion_constant = atof(argv[i+1]); }
		  // else if(flag.compare("-vscale") == 0) { maximum_velocity = atof(argv[i+1]); }
		  // else if(flag.compare("-mrate") == 0) { mutation_rate = atof(argv[i+1]); }
		  // else if(flag.compare("-gref") == 0) { global_refinement = (unsigned int)atof(argv[i+1]); }
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
	// if(number_geometry_files > 1)
	// 	assignGeometryMeshFiles();

	// CHECK PARAMETERS:
	checkParameters(); ///@todo implement this...
}

void ArgParserBase::checkParameters()
{} /// @todo

void ArgParserBase::parseDirectory(const char* directory)
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

			// if(extention.compare("msh") == 0){ mesh_file = direc + fileName; }
			if(extention.compare("dat") == 0)
			{
				std::size_t fnd = fileName.find_last_of("_");
				std::string type = fileName.substr(fnd+1);

				// if(type.compare("velx.dat") == 0){ x_velocity_file = direc + fileName; }
				// if(type.compare("vely.dat") == 0){ y_velocity_file = direc + fileName; }
				if(type.compare("parameters.dat") == 0 || type.compare("para.dat") == 0)
				{
					std::string pfile = direc + fileName;
					parseParameterFile(pfile.c_str());
				}
			// if(type.compare("geo.dat") == 0)
			// {
			// 	geometry_file = direc + fileName;
			// }

		 	} // if .dat file

		} // if not current or parent directory -- should generalize to is not directory... 
	} // while files in directory

  // 3 close:
  closedir(dir);

} // parseDirectory()

void ArgParserBase::parseParameterFile(const char* parameterFile)
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


void ArgParserBase::tokenizeLine(std::string line)
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


void ArgParserBase::parseInitialLocations(std::string line)
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


void ArgParserBase::createLoopedParameter(std::string parameter, 
		std::string start, std::string step, std::string stop)
{
	std::cout << "...adding looped parameter: " << parameter << std::endl;

	Utility::LoopedParameter para;
	para.start_value = std::stod(start);
	para.step_size = std::stod(step);
	para.end_value = std::stod(stop);
	para.name = parameter;

	looped_parameters.push_back(para);

	std::cout << "\t>> start: " << para.start_value << " | step: " << para.step_size
		<< " | stop: " << para.end_value << std::endl;
}


// create looped parameter grid and realization if needed:
void ArgParserBase::realizeLoopedParameters()
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


std::vector<double> ArgParserBase::getGridParameterSet(std::vector<std::vector<double> >all_values, 
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

} // CLOSE NAMESPACE
#endif