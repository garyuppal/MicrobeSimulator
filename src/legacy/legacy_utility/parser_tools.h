#ifndef MICROBESIMULATOR_PARSER_TOOLS_H
#define MICROBESIMULATOR_PARSER_TOOLS_H

#include <ctime>
#include <stdlib.h>
#include <string>
#include <iostream>

namespace MicrobeSimulator{ namespace ParserTools{

enum class ParameterIntialization : int
{
	NONE = 0, DIR, GEO_DIR, FILE
}

// for parameters we're looping over:
struct LoopedParameter{
	double start_value;
	double step_size;
	double end_value;
	std::string name;
};

// Common Parser Methods:
std::string currentDateTime();
bool stringToBool(std::string value);
std::string boolToString(bool value);

// want mutliple tokenizer methods depending on line...
// write tokenizeLine to take a function as an argument?

// can either: assign locations, assign single parameter, or 
// create looped parameter ... but does this method need to be local to class?
	// might be the easiest to just make it local...
// maybe one class is best, and perhaps organize parts with structs....

// IMPLEMENTATION:
// ---------------------------------------------------------------------------
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

}} // close namespaces
#endif