#ifndef MICROBE_SIMULATOR_PARAMETER_HANDLER_H
#define MICROBE_SIMULATOR_PARAMETER_HANDLER_H

#include <deal.II/base/point.h>
using dealii::Point;

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include <boost/archive/basic_archive.hpp>
#include <boost/property_tree/ptree_serialization.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
namespace pt = boost::property_tree;

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "./parameter_patterns.h"
#include "./utility.h"

// TO DO:
// better error checking with parameter patterns
// realize looped parameters ( default parameter before specific assignment)
// help infomation print out
// cleaner print out of parameter tree
// print out of parameter grid // done?

namespace MicrobeSimulator{

// ------------------------------------------------------------------------------
// 	Looping parameter data
// ------------------------------------------------------------------------------

/** \brief Class to handle looped parameter parsing */
class MultiParameterData{
public:
	MultiParameterData(const std::string& psec,
						const std::string& pname,
						const std::string& para_list);

	std::string get_parameter(unsigned int i) const;
	std::vector<std::string> get_parameter_list() const;
	unsigned int get_number_parameters() const;
	std::string get_section() const;
	std::string get_name() const;

	std::unique_ptr<MultiParameterData> clone() const;

private:	
	// need to reassign value of parameters at end of parsing
	std::string section;
	std::string parameter_name; 

	// values ( valid only for doubles at the moment)
	std::string parameter_list;

	std::vector<std::string> fixed_parameter_list;
	/** @todo std::string pattern; // ideally we'd just supply the pattern
	* also deal with looping over selection type parameters */
	unsigned int number_parameters;

	std::vector<std::string> set_fixed_parameter_list() const;
};

/** \brief Constuctor */
MultiParameterData::MultiParameterData(const std::string& psec,
						const std::string& pname,
						const std::string& para_list)
	:
	section(psec),
	parameter_name(pname),
	parameter_list(para_list)
{
	fixed_parameter_list = set_fixed_parameter_list();
	number_parameters = get_parameter_list().size(); // to assign number parameters
}


/** \brief Get ith parameter in loop*/
/** Called by parameter handler to assign specific parameter.
* In case of list type parameters, this returns a list type string.
* The parameter handler class then can manipulate the list string as needed
*/
std::string MultiParameterData::get_parameter(unsigned int i) const
{
	return get_parameter_list()[i];
}

std::vector<std::string> MultiParameterData::get_parameter_list() const
{
	return fixed_parameter_list;
}


/** \brief Get looped list of parameters */
/** This is the main method to find full list of looped parameter.
* For a list of looped parameters, each looped list element is treated 
* separately. If only a subset of list elements are looped, the other
* is treated as a looping parameter of size 1. 
* Also need this method to provide the correct number of parameters.
* That's simply done by taking the size of the output list. For a 
* multiply-looped list, the output is a list of total size given
* the product of sizes of each sub-looped parameter.
*/
std::vector<std::string> MultiParameterData::set_fixed_parameter_list() const
{
	// for doubles for now; in general will depend on pattern
	std::vector<std::string> result, tokens;

	// first check if list (should be able to pass this info from patterns)
	unsigned int n_commas = std::count(parameter_list.begin(), parameter_list.end(), ',');

	// std::cout << "there are " << n_commas << " commas in \n\t<" << parameter_list
	// 	<< ">" << std::endl;

	if(n_commas == 0)
	{
		boost::split(tokens, parameter_list, boost::is_any_of(":"), boost::token_compress_on);
 
 		// assertion still holds here, since n_commas == 0
		assert(tokens.size() == 3); /** @todo instead, count number of :'s, create lists as needed!
										* may need to add a list property to this class */
		const double start = std::stod(tokens[0]);
		const double step = std::stod(tokens[1]);
		const double stop = std::stod(tokens[2]);

		for(double value = start; value <= stop; value += step)
			result.emplace_back(std::to_string(value));
	}
	else if(n_commas == 1)
	{
		// process each subset of list:, assuming each part is looped:
		// first parameter:
		// from first { to ,:
		unsigned int first = parameter_list.find_first_of("{");
		unsigned int middle = parameter_list.find_first_of(",");
		unsigned int last = parameter_list.find_first_of("}");
		std::string first_parameter_list = parameter_list.substr(first+1,middle-first-1);
		boost::trim(first_parameter_list);

		/** Recursively find parameter list.
		* @todo Need to account for case with only a subset looped, not all
		*/
		MultiParameterData fplist = MultiParameterData(section, parameter_name, first_parameter_list);
		std::vector<std::string> firstList = fplist.get_parameter_list();
		unsigned int n_first = fplist.get_number_parameters();
		// std::cout << "there are " << n_first << " first parameters" << std::endl;
		// second parameter:
		std::string second_parameter_list = parameter_list.substr(middle+1,last-middle-1);
		boost::trim(second_parameter_list);

		MultiParameterData splist = MultiParameterData(section, parameter_name, second_parameter_list);
		std::vector<std::string> secondList = splist.get_parameter_list();
		unsigned int n_second = splist.get_number_parameters();
		// std::cout << "there are " << n_second << " second parameters" << std::endl;

		// make all tuples and append to result:
		result.reserve(n_first*n_second);

		for(unsigned int i = 0; i < n_first; ++i)
			for(unsigned int j = 0; j < n_second; ++j)
			{
				std::string currentParameter = "{" 
									+ firstList[i]
									+ ","
									+ secondList[j]
									+ "}";
				result.emplace_back(currentParameter);
			}

		// for(unsigned int i = 0; i < result.size(); ++i)
		// 	std::cout << "\t<" << result[i] << ">" << std::endl;
		// std::cout << "------------------------------------------\n\n\n" << std::endl;
	}
	else
	{
		std::cout << "Multi parameter looping for lists greater than two parameters"
			" not yet impemented" << std::endl;
		assert(false);
	}

	return result; // a tuple type if input is as such (1){a,b,c} (2){d,e,f} ...
}


	/* parameter list is whole list from { to } with list seperated by , 's
	* ... want to store each , separated thing separately, and reassign as a list...
	*/

	// std::cout << "Tokens for " << section + "." + parameter_name << ": " << std::endl;
	// for(unsigned int i = 0; i < tokens.size(); ++i)
	// 	std::cout << "<" << tokens[i] << ">" << std::endl;
	// std::cout << "............." << std::endl << std::endl;


unsigned int MultiParameterData::get_number_parameters() const
{
	return number_parameters;
}

std::string MultiParameterData::get_section() const
{
	return section;
}

std::string MultiParameterData::get_name() const
{
	return parameter_name;
}

std::unique_ptr<MultiParameterData> MultiParameterData::clone() const
{
	return std::unique_ptr<MultiParameterData>(new 
		MultiParameterData(section,parameter_name,parameter_list));
}















// ------------------------------------------------------------------------------
// 	PARAMETER HANDLER
// ------------------------------------------------------------------------------


/** \brief Parameter handling class */
class ParameterHandler{
public:
	ParameterHandler(const std::string& pfile, unsigned int jid = 0);

	void declare_entry(const std::string& parameter_name,
							const std::string& default_value,
							const Patterns::PatternBase& data_type,
							const std::string& help_msg = "");

	void parse_parameter_file(); 

	// TRAVERSAL:
	void enter_subsection(const std::string& subsection);
	void leave_subsection();

	// ACCESSORS:

	// base:
	std::string get(const std::string& entry) const; 
	std::string get_string(const std::string& entry) const;  // same as get
	double get_double(const std::string& entry) const;
	int get_int(const std::string& entry) const;
	unsigned int get_unsigned(const std::string& entry) const;
	bool get_bool(const std::string& entry) const;

	// access by section:
	std::string get(const std::string& section, const std::string& entry) const; 
	std::string get_string(const std::string& section, const std::string& entry) const;
	double get_double(const std::string& section, const std::string& entry) const;
	int get_int(const std::string& section, const std::string& entry) const;
	unsigned int get_unsigned(const std::string& section, const std::string& entry) const;
	bool get_bool(const std::string& section, const std::string& entry) const;

	// list access:
	std::vector<std::string> get_list(const std::string& entry) const;
	std::vector<Point<2> > get_point_list(const std::string& entry) const; // generalize to n-d?
	std::vector<std::string> get_string_vector(const std::string& entry) const;
	std::vector<double> get_double_vector(const std::string& entry) const;

	// list access by section:
	std::vector<std::string> get_list(const std::string& section, const std::string& entry) const;
	std::vector<Point<2> > get_point_list(const std::string& section, const std::string& entry) const;
	std::vector<std::string> get_string_vector(const std::string& section, const std::string& entry) const;
	std::vector<double> get_double_vector(const std::string& section, const std::string& entry) const;

	// OUTPUT:
	void print(std::ostream& out) const;
	void print_simple(std::ostream& out) const;
	void printLoopedParameterGrid(std::ostream& out) const;
	// std::ostream& print_parameters(std::ostream& out, const OutputStyle style) const;

private:
	static const char path_separator = '.'; // for traversing ptree
	static const int invalid_multiparameter_index = -1;

	std::vector<std::string> current_path;

	std::string parameter_file;
	unsigned int job_ID; // possibly needed for looped parameters *** will need to know
		// which parameters are looped, and how many, in advance ...

	std::unique_ptr<pt::ptree> entries; // does this need to be allocated on the heap?

	std::vector<std::unique_ptr<const Patterns::PatternBase> > patterns;

	std::vector<std::unique_ptr<const MultiParameterData> > looped_parameters;

	std::string	get_current_path() const;
	std::string	get_current_full_path(const std::string &name) const;

	// looped tools:
	void assign_looped_parameters();
	unsigned int get_grid_size() const; 
	std::vector<unsigned int> get_local_job_indices() const;
	std::vector<unsigned int> get_local_job_indices(unsigned int jid) const;
	
	// file parse tools:
	void parse_file(std::ifstream& infile);
	void process_line(std::string line);
	std::string process_possible_looped_value(std::string entry, 
		std::string value);

		// ADD ENTRIES:
	void set(const std::string& entry_path, const std::string& value);
	void reset(const std::string& section,
				const std::string& name,
				const std::string& value);
	// print:
	void print_recursive(std::ostream& out, 
		 const pt::ptree &pt, int level) const;
	void print_recursive_simple(std::ostream& out, 
		 const pt::ptree &pt, int level) const;

};

// IMPL
// ----------------------------------------------------------------------------
ParameterHandler::ParameterHandler(const std::string& pfile, unsigned int jid)
	:
	parameter_file(pfile),
	job_ID(jid),
	entries(new pt::ptree())
{}

void ParameterHandler::declare_entry(const std::string& parameter_name,
						const std::string& default_value,
						const Patterns::PatternBase& pattern,
						const std::string& help_msg)
{
	entries->put(get_current_full_path(parameter_name) + path_separator + "value",
                default_value);
    entries->put(get_current_full_path(parameter_name) + path_separator + "default_value",
                default_value);
    entries->put(get_current_full_path(parameter_name) + path_separator + "documentation",
                help_msg);
 
   patterns.reserve(patterns.size() + 1);
   patterns.emplace_back(pattern.clone());
   entries->put(get_current_full_path(parameter_name) + path_separator + "pattern",
                static_cast<unsigned int>(patterns.size() - 1)); // record pattern index

   // // also store the description of
   // // the pattern. we do so because we
   // // may wish to export the whole
   // // thing as XML or any other format
   // // so that external tools can work
   // // on the parameter file; in that
   // // case, they will have to be able
   // // to re-create the patterns as far
   // // as possible
   // entries->put(get_current_full_path(entry) + path_separator +
   //                "pattern_description",
   //              patterns.back()->description());
 	
 	assert(pattern.match(default_value));

   // // as documented, do the default value checking at the very end
   // AssertThrow(pattern.match(default_value),
   //             ExcValueDoesNotMatchPattern(default_value,
   //                                         pattern.description()));
}

// TRAVERSAL:
void ParameterHandler::enter_subsection(const std::string& subsection)
{
	// std::cout << "entering section: <" << subsection << "> " << std::endl;

	current_path.emplace_back(subsection);
}


void ParameterHandler::leave_subsection()
{
	// check if already at top
	assert(current_path.size() != 0);

	current_path.pop_back();
}

void ParameterHandler::set(const std::string& entry_string,
			const std::string& new_value)
{
	   // resolve aliases before looking up the correct entry
   std::string path = get_current_full_path(entry_string);
   // if (entries->get_optional<std::string>(path + path_separator + "alias"))
   //   path = get_current_full_path(
   //     entries->get<std::string>(path + path_separator + "alias"));
 
   // get the node for the entry. if it doesn't exist, then we end up
   // in the else-branch below, which asserts that the entry is indeed
   // declared

   if (entries->get_optional<std::string>(
   		get_current_full_path(entry_string) + path_separator + "value"))
     {
       // verify that the new value satisfies the provided pattern
       // const unsigned int pattern_index =
       //   entries->get<unsigned int>(path + path_separator + "pattern"); 

     //***  // how do we store patterns?
       // AssertThrow(patterns[pattern_index]->match(new_value),
       //             ExcValueDoesNotMatchPattern(new_value,
       //                                         entries->get<std::string>(
       //                                           path + path_separator +
       //                                           "pattern_description")));
 
       // // then also execute the actions associated with this
       // // parameter (if any have been provided)
       // const boost::optional<std::string> action_indices_as_string =
       //   entries->get_optional<std::string>(path + path_separator + "actions");
       // if (action_indices_as_string)
       //   {
       //     std::vector<int> action_indices = Utilities::string_to_int(
       //       Utilities::split_string_list(action_indices_as_string.get()));
       //     for (const unsigned int index : action_indices)
       //       if (actions.size() >= index + 1)
       //         actions[index](new_value);
       //   }

     	std::string usingValue = process_possible_looped_value(entry_string, new_value);
 
       // finally write the new value into the database
       entries->put(path + path_separator + "value", usingValue);
     }
   else{
   		std::cout << "parameter <" 
   			<< path <<  "> was not declared" << std::endl;
     assert(false); // AssertThrow(false, ExcEntryUndeclared(entry_string));
   }
}

std::string ParameterHandler::process_possible_looped_value(std::string entry, 
		std::string value)
{
	// check if loooped (if there's colons)
	std::size_t found = value.find_first_of(":");
	if(found!=std::string::npos)
	{
		// create looped parameter:
		looped_parameters.emplace_back( std::unique_ptr<MultiParameterData>(
			new MultiParameterData(get_current_path(),entry,value)) );

	}

	// double error_finish_this_function; // need to know section ...
	// else
		return value;
}

void ParameterHandler::assign_looped_parameters()
{
	std::vector<unsigned int> local_indices = get_local_job_indices();

	for(unsigned int i = 0; i < looped_parameters.size(); ++i)
	{
		std::string value = looped_parameters[i]->get_parameter(local_indices[i]);
		std::string section = looped_parameters[i]->get_section();
		std::string name = looped_parameters[i]->get_name();

		reset(section, name, value);
	}
}

void ParameterHandler::reset(const std::string& section,
				const std::string& name,
				const std::string& value)
{
	// if gotten this far, value should exist
	// std::string entry_string = section 
	// 	+ path_separator + name + path_separator + "value";

	// std::cout << "resetting <" << entry_string << ">" << std::endl;
	entries->put(section 
		+ path_separator + name + path_separator + "value", value);
}

unsigned int ParameterHandler::get_grid_size() const
{
	unsigned int result = 1;
	for(unsigned int i = 0; i < looped_parameters.size(); ++i)
		result *= looped_parameters[i]->get_number_parameters();

	return result;
}

std::vector<unsigned int> ParameterHandler::get_local_job_indices() const
{
	return get_local_job_indices(job_ID);
	// assert(job_ID > 0); 

	// std::vector<unsigned int> local_indices;
	// local_indices.reserve(looped_parameters.size());

	// const unsigned int total_grid_size = get_grid_size();
	// const unsigned int this_job = ((job_ID -1) % total_grid_size);

	// // get index based on job id, go through order of looped parameters
	// unsigned int divide_factor = 1;
	// for(unsigned int i = 0; i < looped_parameters.size(); ++i)
	// {
	// 	const unsigned int current_size = looped_parameters[i]->get_number_parameters();
	// 	const unsigned int current_index = ((this_job/divide_factor) % current_size);

	// 	local_indices.emplace_back(current_index); 
	// 	divide_factor *= current_size;
	// }

	// return local_indices;
}

std::vector<unsigned int> ParameterHandler::get_local_job_indices(unsigned int jid) const
{
	assert(jid > 0); 

	std::vector<unsigned int> local_indices;
	local_indices.reserve(looped_parameters.size());

	const unsigned int total_grid_size = get_grid_size();
	const unsigned int this_job = ((jid -1) % total_grid_size);

	// get index based on job id, go through order of looped parameters
	unsigned int divide_factor = 1;
	for(unsigned int i = 0; i < looped_parameters.size(); ++i)
	{
		const unsigned int current_size = looped_parameters[i]->get_number_parameters();
		const unsigned int current_index = ((this_job/divide_factor) % current_size);

		local_indices.emplace_back(current_index); 
		divide_factor *= current_size;
	}

	return local_indices;
}

void ParameterHandler::printLoopedParameterGrid(std::ostream& out) const
{
	out << std::endl << std::endl
		<< Utility::medium_line
		<< std::endl
		<< "\t\t PARAMETER GRID"
		<< std::endl
		<< Utility::medium_line
		<< std::endl;
	out << "Job ID";
	for(unsigned int i = 0; i < looped_parameters.size(); ++i)
	{
		out << "\t" << looped_parameters[i]->get_name();
	}
	out << "\t\t THIS" << std::endl;

	const unsigned int total_grid_size = get_grid_size();
	for(unsigned int i = 1; i <= total_grid_size; ++i) // for all possible parameters
	{
		std::vector<unsigned int> job_indices = get_local_job_indices(i);
		out << "  " << i;
		for(unsigned int j = 0; j < job_indices.size(); ++j)
		{
			out << "\t" << looped_parameters[j]->get_parameter(job_indices[j]);
		}
		if( ((job_ID-1)%total_grid_size) == (i-1) )
		{
			out << "\t\t X";
		}
		out << std::endl;
	}

	out << Utility::medium_line
		<< std::endl << std::endl << std::endl;

}


// ACCESSORS:
// std::string ParameterHandler::get(const std::vector<std::string>& subpath, 
// 							const std::string& entry_string) const
// {
// 	if (boost::optional<std::string> value = entries->get_optional<std::string>(
// 			subpath + "." + entry_string + "." +
//         	// get_current_full_path(subpath, entry_string) +
//         	path_separator + "value"))
//     	return value.get();
//    else
//    {
//    		std::cout << "ERROR: entry not declared" << std::endl;
//        // Assert(false,
//        //        ExcEntryUndeclared(demangle(
//        //          get_current_full_path(entry_subsection_path, entry_string))));
//        return "";
//     }
// }

std::string ParameterHandler::get(const std::string& entry_string) const
{
   // assert that the entry is indeed
   // declared
   if (boost::optional<std::string> value = entries->get_optional<std::string>(
         get_current_full_path(entry_string) + path_separator + "value"))
     return value.get();
   else
     {
     	std::cout << "parameter <" << entry_string << "> is undeclared" << std::endl;
     	assert(false);
       // Assert(false, ExcEntryUndeclared(entry_string));
       return "";
     }
}

std::string ParameterHandler::get_string(const std::string& entry) const
{
	return get(entry);
}

double ParameterHandler::get_double(const std::string& entry) const
{
	///@todo
	// need to do pattern checking ...
	return std::stod(get(entry));
}

int ParameterHandler::get_int(const std::string& entry) const
{
	///@todo
	// need to do pattern checking ...
	return std::stoi(get(entry));
}

unsigned int ParameterHandler::get_unsigned(const std::string& entry) const
{
	///@todo
	// need to do pattern checking...
	return std::stoul(get(entry));
}

bool ParameterHandler::get_bool(const std::string& entry) const
{
	///@todo catch exception possibly thrown by string to bool
	return Utility::stringToBool(get(entry));
}

// SECITON ACCESS:
std::string ParameterHandler::get(const std::string& section, const std::string& entry) const
{
	   // assert that the entry is indeed
   // declared
   if (boost::optional<std::string> value = entries->get_optional<std::string>(
   		section + get_current_full_path(entry) + path_separator + "value"))
     return value.get();
   else
     {
     	std::cout << "parameter <" << entry << "> is undeclared" << std::endl;
     	assert(false);
       // Assert(false, ExcEntryUndeclared(entry_string));
       return "";
     }
}

std::string ParameterHandler::get_string(const std::string& section, const std::string& entry) const
{
	return get(section, entry);
}

double ParameterHandler::get_double(const std::string& section, const std::string& entry) const
{
	// const std::string value_string = get(section, entry);
	// const double value = std::stod(value_string);
	// std::cout << "retrieving <" << entry << ">: " << " value string: >"
	// 	<< value_string << "> converted to double: " << value << std::endl;
	return std::stod(get(section,entry));
}

int ParameterHandler::get_int(const std::string& section, const std::string& entry) const
{
	return std::stoi(get(section,entry));
}

unsigned int ParameterHandler::get_unsigned(const std::string& section, const std::string& entry) const
{
	return std::stoul(get(section, entry));
}

bool ParameterHandler::get_bool(const std::string& section, const std::string& entry) const
{
	return Utility::stringToBool(get(section, entry));
}



// LIST ACCESS:

std::vector<std::string> ParameterHandler::get_list(const std::string& entry) const
{
	return Utility::split_string_list(get(entry), ",");
}

std::vector<Point<2> > ParameterHandler::get_point_list(const std::string& entry) const
{
	std::vector<std::string> locations = Utility::split_string_list(get(entry), "");
	std::vector<Point<2> > points;

	for(unsigned int i = 0; i < locations.size(); ++i)
	{
		if( !boost::iequals(locations[i],",") )
		{
			std::vector<std::string> numbers = Utility::split_string_list(locations[i],",");
			if(numbers.size() == 2)
			{
				Point<2> p;
				for(unsigned int i = 0; i < numbers.size(); ++i)
					p[i] = std::stod(numbers[i]);
				points.emplace_back(p);
			}
		} // if not empty token
	} // for given points

	return points;
}

std::vector<std::string> ParameterHandler::get_string_vector(const std::string& entry) const
{
	return get_list(entry);
}

std::vector<double> ParameterHandler::get_double_vector(const std::string& entry) const
{
	std::vector<std::string> value_list = get_list(entry);

	std::vector<double> values;
	values.reserve(value_list.size());

	for(unsigned int i = 0 ; i < value_list.size(); ++i)
		values.emplace_back(std::stod(value_list[i]));

	return values;
}

// list access by section:
std::vector<std::string> ParameterHandler::get_list(const std::string& section, 
		const std::string& entry) const
{
	return Utility::split_string_list(get(section, entry), ",");
}

std::vector<Point<2> > ParameterHandler::get_point_list(const std::string& section, 
		const std::string& entry) const
{
	std::vector<std::string> locations = Utility::split_string_list(get(section, entry), "");
	std::vector<Point<2> > points;

	for(unsigned int i = 0; i < locations.size(); ++i)
	{
		if( !boost::iequals(locations[i],",") )
		{
			std::vector<std::string> numbers = Utility::split_string_list(locations[i],",");
			if(numbers.size() == 2)
			{
				Point<2> p;
				for(unsigned int i = 0; i < numbers.size(); ++i)
					p[i] = std::stod(numbers[i]);
				points.emplace_back(p);
			}
		} // if not empty token
	} // for given points

	return points;
}

std::vector<std::string> ParameterHandler::get_string_vector(const std::string& section, 
		const std::string& entry) const
{
	return get_list(section, entry);
}

std::vector<double> ParameterHandler::get_double_vector(const std::string& section, 
	const std::string& entry) const
{
	std::vector<std::string> value_list = get_list(section, entry);

	std::vector<double> values;
	values.reserve(value_list.size());

	for(unsigned int i = 0 ; i < value_list.size(); ++i)
		values.emplace_back(std::stod(value_list[i]));

	return values;
}

// HELPER METHODS:
// -------------------------------------------------------------------------------------------


std::string ParameterHandler::get_current_path() const
{
	std::string path;

	if(current_path.size() > 0)
	{
		for(unsigned int i = 0; i < current_path.size()-1; ++i)
		{
			path += current_path[i] + path_separator;
		}
		path += current_path[current_path.size()-1];
	}
	else
		return "";

	return path;
}

std::string ParameterHandler::get_current_full_path(const std::string &name) const
{
	return get_current_path() + path_separator + name;
}

void ParameterHandler::parse_parameter_file()
{
	std::ifstream infile(parameter_file);
	parse_file(infile);

	assign_looped_parameters();
}


void ParameterHandler::parse_file(std::ifstream& infile)
{
	std::cout << "parsing file..." << std::endl;

	std::string line;

	unsigned int count = 0;
	// read in each line:
	while(std::getline(infile,line)){
		++count;
		// std::cout << "processing line: " << count << std::endl;
		// tokenize and process each line:
		process_line(line);
	}
}

void ParameterHandler::process_line(std::string line)
{
	std::vector<std::string> tokens;

	// remove leading whitespace:
	std::string trimmed_line = line; 
	boost::trim(trimmed_line);

	// ignore comments:
	if(trimmed_line.size() == 0)
		return;
	else if(trimmed_line.size() >= 2)
		if(trimmed_line[0]=='/' && trimmed_line[1]=='/')
			return;

	// tokenize line:
	boost::split(tokens, trimmed_line, boost::is_any_of(" "), boost::token_compress_on); 

	if(!tokens.empty())
	{
		if(tokens[0].compare("subsection") == 0)
		{
			// get subsection string by combining remaining tokens:
			std::string subsection = 
				Utility::combine_vector_string(tokens,"subsection","//"); 
			enter_subsection(subsection);

			// std::cout << "current_path is " << get_current_path() << std::endl;
		}
		else if(tokens[0].compare("set") == 0)
		{
			std::string variable_name = 
				Utility::combine_vector_string(tokens,"set","=");
			std::string value = 
				Utility::combine_vector_string(tokens,"=","//");

			set(variable_name,value);
		}
		else if(tokens[0].compare("end") == 0)
		{
			leave_subsection();
		}
		else
		{
			std::cout << "invalid keyword <" << tokens[0] << "> " << std::endl;
			assert(false);
 		}
	}
} // process_line()


std::string indent(int level) {
  std::string s; 
  for (int i=0; i<level; i++) s += "  ";
  return s; 
} 

void ParameterHandler::print(std::ostream& out) const
{
	out << "\n\n" << Utility::medium_line << std::endl 
		<< "\t\t PARAMETERS:" << std::endl
		<< Utility::medium_line << std::endl;

	int level = 0;

	print_recursive(out,*entries,level);
}

void ParameterHandler::print_simple(std::ostream& out) const
{
	out << "\n\n" << Utility::medium_line << std::endl 
		<< "\t\t PARAMETERS:" << std::endl
		<< Utility::medium_line << std::endl;

	int level = 0;

	print_recursive_simple(out,*entries,level);
	out << Utility::medium_line
		<< std::endl << std::endl << std::endl;
}

void ParameterHandler::print_recursive(std::ostream& out, 
	 const pt::ptree &tree, int level) const
{

	if (tree.empty()) {
		out << "\""<< tree.data()<< "\"";
	}
	else {
		if(level)
			 out << std::endl; 

		out << indent(level) << "{" << std::endl;     

		for (pt::ptree::const_iterator pos = tree.begin(); 
				pos != tree.end(); ) {
			out << indent(level+1) << "\"" << pos->first << "\": "; 

			print_recursive(out, pos->second, level + 1); 
			++pos; 
		
			if (pos != tree.end()) {
				out << ","; 
			}
			
			out << std::endl;
		} 

		out << indent(level) << " }";     
	}

	return; 
}

void ParameterHandler::print_recursive_simple(std::ostream& out, 
	 const pt::ptree &tree, int level) const
{
	if (tree.empty()) 
	{
		/** @todo adjust precision appropriately... //std::setprecission(10) << */
		out << "\""<<  tree.data()<< "\"" << std::endl;
	}
	else 
	{
		if(!(tree.begin()->second).empty())
			out << std::endl;
		else
			out << "\t\t\t";

		for (pt::ptree::const_iterator pos = tree.begin(); 
				pos != tree.end(); ) 
		{
			out << indent(2*level+2) << pos->first << ": "; 
			print_recursive_simple(out, pos->second, level + 1); 
		
			// skip over non-value print-outs:
			pt::ptree::const_iterator pos_value = pos;
			++(++(++(++pos_value)));
			if((pos->second).empty())
				if (pos_value == tree.end()) 
					break;

			++pos; 
		} 
	}

	return; 
} // print_recursive_simple()

} // CLOSE NAMESPACE
#endif