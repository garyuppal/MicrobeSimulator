#ifndef MICROBESIMULATOR_COMMAND_LINE_PARSER_H	
#define MICROBESIMULATOR_COMMAND_LINE_PARSER_H	

#include <iostream>
#include <fstream>
// to explore directory:
#include <dirent.h>
#include <sys/stat.h> 

#include "./utility.h"

namespace MicrobeSimulator{

class CommandLineParameters{
public:
	CommandLineParameters(int argc, char** argv);

	void create_output_directory();
	void print(std::ostream& out) const;

	// ACCESSORS:
	std::string getParameterFile() const {return parameter_file;}
	unsigned int getJobID() const {return job_ID;}
	unsigned int getDimension() const {return dimension;}
	std::string getOutputTag() const {return output_tag;}
	std::string getOutputDirectory() const {return output_directory;}
private:
	std::string parameter_file;
	unsigned int job_ID;
	unsigned int dimension;
	std::string output_tag;
	std::string output_directory;
	bool create_directory;

	void parse(int argc, char** argv);
};

// IMPL
// -----------------------------------------------------------------------
CommandLineParameters::CommandLineParameters(int argc, char** argv)
	:
	parameter_file(""),
	job_ID(1),
	dimension(2),
	output_tag(""),
	output_directory("./"),
	create_directory(true)
{
	parse(argc,argv);

	if(create_directory)
		create_output_directory();
}

void CommandLineParameters::parse(int argc, char** argv)
{
	std::ostringstream message;
		    message << std::endl
		      << "***UNKNOWN FLAG*** \n"
		      << "Usage: ./program [-flags] (file or directory) \n"
		      << "\t -f \"parameter file\" \n"
		      << "\t -p \"parameter file\" \n"
		      << "\t -o \"output tag\" \n"
		      << "\t -id \"job id\" \n"
		      << "\t -d \"dimension\" \n" 
		      << "\t -dim \"dimension\" \n" 
		      << "\t -od \"0 or 1\" \n" << std::endl;

	std::cout << "\n\n...Parsing command line\n" << std::endl;

	if((argc % 2) == 1) // program + flags + values == odd
	{
		for(int i = 1; i < argc; i += 2)
		{
			std::string flag = argv[i]; 

			if(flag.compare("-f") == 0) { parameter_file = argv[i+1]; }
			else if(flag.compare("-p") == 0) { parameter_file = argv[i+1]; }
			else if(flag.compare("-o") == 0) { output_tag = argv[i+1]; output_tag += "_";}
			else if(flag.compare("-id") == 0 ) { job_ID = atoi(argv[i+1]); }
			else if(flag.compare("-d") == 0) { dimension = atoi(argv[i+1]); }
			else if(flag.compare("-dim") == 0) { dimension = atoi(argv[i+1]); }
			else if(flag.compare("-od") == 0) { create_directory = Utility::stringToBool(argv[i+1]); }
			else
			{  
				throw std::invalid_argument(message.str()); 
			}
		} // for all inputs
	} // if even inputs (-flag value)
	else 
	{
		throw std::invalid_argument(message.str()); 
	} // if correct number of inputs
} // parse()


void CommandLineParameters::create_output_directory()
{
	// create name from tag:
	output_directory =  "./data_" + output_tag + "jobID_" 
		+  std::to_string(job_ID) + "_" + Utility::currentDateTime();

	// make directory and check if successful:
	int check;
	check = mkdir(output_directory.c_str(), 0777);

	if(check == 0)
		std::cout << "...Data output directory created: " << output_directory << std::endl;
	else
		throw std::runtime_error("\n*** FAILED TO CREATE OUTPUT DIRECTORY ***\n"); 
} // create_output_directory()

void CommandLineParameters::print(std::ostream& out) const
{
	out << std::endl << Utility::medium_line << std::endl
		<< "\t COMMAND LINE ARGUMENTS" << std::endl
		<< Utility::medium_line << std::endl
		<< "Parameter File: " << parameter_file << std::endl
		<< "Job ID: " << job_ID << std::endl
		<< "Dimension: " << dimension << std::endl
		<< "Output Directory: " << output_directory << std::endl
		<< Utility::medium_line << std::endl
		<< std::endl << std::endl;
}

} // CLOSE NAMESPACE
#endif