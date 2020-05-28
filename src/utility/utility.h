#ifndef MICROBE_SIMULATOR_UTILITY_H
#define MICROBE_SIMULATOR_UTILITY_H

// #define NDEBUG
#include <assert.h>

#include <chrono>
#include <string>

// boost library:
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

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

namespace Utility{
	static constexpr double PI =  3.14159265358979323846; 
	static constexpr double PI_2 = 1.57079632679489661923;
	static const std::string short_line = 
		"-----------------------------------------------------";
	static const std::string medium_line = 
		"-----------------------------------------------------------------------------";
	static const std::string long_line = 
		"-----------------------------------------------------"
				"-----------------------------------------------------";


	struct LoopedParameter{
		double start_value;
		double step_size;
		double end_value;
		std::string name;
	};

	// struct Timer{
	// 	std::chrono::time_point<std::chrono::steady_clock> start, end;
	// 	std::chrono::duration<float> duration;

	// 	Timer()
	// 	{
	// 		start = std::chrono::high_resolution_clock::now();
	// 	}

	// 	~Timer()
	// 	{
	// 		end = std::chrono::high_resolution_clock::now();
	// 		duration = end-start;

	// 		std::cout << "Timer took " << duration << "s" << std::endl;
	// 	}
	// };

	double getRand()
	{
		return ((double) rand() / (RAND_MAX));
	}

	// DATE-TIME FUNCTION:
	std::string currentDateTime() 
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
		if( (boost::iequals(value, "TRUE")) || (value.compare("1") == 0) )
			return true;
		else if( (boost::iequals(value, "FALSE")) || (value.compare("0") == 0) )
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

	std::string combine_vector_string(const std::vector<std::string>& raw_tokens,
									const std::string& begin_token ) //,
									// const std::string& end_token = '\n')
	{
		std::string combined_string;

		unsigned int i = 0; 

		// trim tokens:
		std::vector<std::string> tokens = raw_tokens;
		for(unsigned int i = 0; i < tokens.size(); ++i)
			boost::trim(tokens[i]);

		// move up until begin token found:
		while( tokens[i].compare(begin_token) != 0  && i < tokens.size())
			++i;
		++i; // skip the begin token itself

		for(; i < tokens.size()-1; ++i)
			combined_string += tokens[i] + " ";
		combined_string += tokens[tokens.size()-1];

		return combined_string;
	}

	std::string combine_vector_string(const std::vector<std::string>& raw_tokens,
									const std::string& begin_token,
									const std::string& end_token)
	{
		std::string combined_string;

		unsigned int i = 0; 

		// trim tokens:
		std::vector<std::string> tokens = raw_tokens;
		for(unsigned int i = 0; i < tokens.size(); ++i)
			boost::trim(tokens[i]);

		// move up until begin token found:
		while( tokens[i].compare(begin_token) != 0  && i < tokens.size())
			++i;
		++i; // skip the begin token itself

		for(; i < tokens.size(); ++i)
		{
			if(tokens[i].compare(end_token) == 0)
				break;

			combined_string += tokens[i] + " ";
		}

		boost::trim(combined_string);
		return combined_string;
	}

	bool has_only_whitespace(std::istream &in)
	{
		while (in)
		{
			char c;

			// skip if we've reached the end of
			// the line
			if (!(in >> c))
				break;

			if ((c != ' ') && (c != '\t'))
				return false;
		}
		return true;
	}

	std::vector<std::string> split_string_list(
		const std::string& test_string_list, const std::string& separator)
	{
		std::string trimmed_line(test_string_list);
		boost::trim(trimmed_line);

		std::vector<std::string> split_list;

		std::string splitters = separator + "{}";
		boost::split(split_list, 
			trimmed_line, 
			boost::is_any_of(splitters), 
			boost::token_compress_on); 

		// remove padded whitespace and empty tokens:
		for(auto it = split_list.begin(); it != split_list.end(); )
		{
			boost::trim(*it);

			if(it->empty())
				it = split_list.erase(it);
			else
				++it;
		}

		return split_list;	
	}

} // CLOSE NAMESPACE Utility

} // CLOSE NAMESPACE
#endif