#ifndef MICROBESIMULATOR_PARAMETER_PATTERNS_H
#define MICROBESIMULATOR_PARAMETER_PATTERNS_H

#include <memory>

#include "./utility.h"

namespace MicrobeSimulator{ namespace Patterns{

// namespace
// {

// } // namespace

// INTERFACE CLASS:
class PatternBase{
public:
	PatternBase() {}
	virtual ~PatternBase() {}

	virtual std::unique_ptr<PatternBase> clone() const=0;
	virtual bool match(const std::string& parameter) const=0;
protected:
};

// --- ANYTHING ---
// ----------------------------------------------------
class Anything : public PatternBase{
public:
	Anything() {}
	std::unique_ptr<PatternBase> clone() const override;
	bool match(const std::string& parameter) const override;
};

// IMPL
// -------------------------------------------------
std::unique_ptr<PatternBase> Anything::clone() const 
{
	return std::unique_ptr<PatternBase>(new Anything());
}

bool Anything::match(const std::string& /* parameter */) const 
{
	return true;
}

// --- SELECTION TYPE ---
// ----------------------------------------------------
class Selection : public PatternBase{
public:
	Selection(const std::string& c)
	: choices(c) {}

	std::unique_ptr<PatternBase> clone() const override;
	bool match(const std::string& parameter) const override;
private:
	static const char separator = '|';
	std::string choices;
};

// IMPL
// -------------------------------------------------
std::unique_ptr<PatternBase> Selection::clone() const 
{
	return std::unique_ptr<PatternBase>(new Selection(choices));
}

bool Selection::match(const std::string &test_string) const
{
	std::string tmp(choices);

	// remove whitespace at beginning
	while ((tmp.length() != 0) && (std::isspace(tmp[0])))
		tmp.erase(0, 1);

	// check the different possibilities
	while (tmp.find(separator) != std::string::npos)
	{
		if (test_string == std::string(tmp, 0, tmp.find(separator)))
			return true;

		tmp.erase(0, tmp.find(separator) + 1);
	}

	// remove whitespace at the end
	while ((tmp.length() != 0) && (std::isspace(*(tmp.end() - 1))))
		tmp.erase(tmp.end() - 1);

	// check last choice, not finished by |
	if (test_string == tmp)
		return true;

	// not found
	return false;
}

// --- DOUBLE ---
// ----------------------------------------------------------------
class Double : public PatternBase{
public:
	Double() {}

	std::unique_ptr<PatternBase> clone() const override;
	bool match(const std::string& parameter) const override;
};

// IMPL
// --------------------------------------------------

std::unique_ptr<PatternBase> Double::clone() const 
{
	return std::unique_ptr<PatternBase>(new Double());
}

bool Double::match(const std::string &test_string) const
{
	std::istringstream str(test_string);

	double d;
	str >> d;
	if (str.fail())
		return false;

	if (!Utility::has_only_whitespace(str))
		return false;
	// // check whether valid bounds
	// // were specified, and if so
	// // enforce their values
	// if (lower_bound <= upper_bound)
	// 	return ((lower_bound <= d) && (upper_bound >= d));
	// else 
	/// not enforcing bounds
	return true;
}

// --- INTEGER ---
// ----------------------------------------------------------------
class Integer : public PatternBase{
public:
	Integer() {}

	std::unique_ptr<PatternBase> clone() const override;
	bool match(const std::string& parameter) const override;
};

// IMPL
// --------------------------------------------
std::unique_ptr<PatternBase> Integer::clone() const 
{
	return std::unique_ptr<PatternBase>(new Integer());
}

bool Integer::match(const std::string &test_string) const
{
	std::istringstream str(test_string);

	int i;
	if (!(str >> i))
		return false;

	if (!Utility::has_only_whitespace(str))
		return false;

	return true;
}

// --- UNSIGNED ---
// ----------------------------------------------------------------
class Unsigned : public PatternBase{
public:
	Unsigned() {}

	std::unique_ptr<PatternBase> clone() const override;
	bool match(const std::string& parameter) const override;
};

// IMPL
// --------------------------------------------
std::unique_ptr<PatternBase> Unsigned::clone() const 
{
	return std::unique_ptr<PatternBase>(new Unsigned());
}

bool Unsigned::match(const std::string &test_string) const
{
	std::istringstream str(test_string);

	unsigned int i;
	if (!(str >> i))
		return false;

	if (!Utility::has_only_whitespace(str))
		return false;

	return true;
}

// --- BOOL ---
// ----------------------------------------------------------------
class Bool : public PatternBase{
public:
	Bool() {}

	std::unique_ptr<PatternBase> clone() const override;
	bool match(const std::string& parameter) const override;
};

// IMPL
// --------------------------------------------
std::unique_ptr<PatternBase> Bool::clone() const 
{
	return std::unique_ptr<PatternBase>(new Bool());
}

bool Bool::match(const std::string &test_string) const
{
	std::string tmp(test_string);
	boost::trim(tmp);

	if( boost::iequals(tmp, "TRUE") || (tmp.compare("1") == 0) ||
			 boost::iequals(tmp, "FALSE") || (tmp.compare("0") == 0) )
		return true;

	// not found
	return false;
}


// --- LIST ---
// ----------------------------------------------------------------
class List : public PatternBase{
public:
	List(const PatternBase& p, 
		const std::string& sep = ",") 
	: base_pattern(p.clone()), separator(sep) {}

	std::unique_ptr<PatternBase> clone() const override;
	bool match(const std::string& parameter) const override;

	const std::string& get_separator() const;  // can use to return vector of strings
private:
	std::unique_ptr<PatternBase> base_pattern;
	std::string separator;
};

// IMPL
// --------------------------------------------
std::unique_ptr<PatternBase> List::clone() const 
{
	return std::unique_ptr<PatternBase>(new List(*base_pattern, separator));
}

bool List::match(const std::string &test_string_list) const
{
	const std::vector<std::string> split_list =
		Utility::split_string_list(test_string_list, separator);

	// if ((split_list.size() < min_elements) ||
	// 		(split_list.size() > max_elements))
	// 	return false;

	// check the different possibilities
	for (const std::string &string : split_list)
		if (base_pattern->match(string) == false)
			return false;

	return true;
}
 

}} // CLOSE NAMESPACE
#endif