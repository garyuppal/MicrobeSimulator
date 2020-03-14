#ifndef MICROBESIMULATOR_BOOL_FUNCTIONS_H
#define MICROBESIMULATOR_BOOL_FUNCTIONS_H

/**
* RETURN BOOLEAN INSTEAD OF DOUBLE !!!
* NEED to also know number of switches ....
*/

#include <iostream> 
#include <cmath>

namespace MicrobeSimulator{ namespace BoolFunctions{

// VIRTUAL BASE FOR TIME FUNCTIONS:
class BoolFunction{
public:
	BoolFunction() {}
	virtual ~BoolFunction() {}

	virtual bool value(double time) const = 0;
};


// --------------------------------------------------------------------
// SQUARE PULSES:
class SquareWave : public BoolFunction{
public:
	SquareWave();
	bool value(double time) const override;

	double getOnPeriod() const;
	double getOffPeriod() const;

	void setOnPeriod(double on);
	void setOffPeriod(double off);
private:
	double on_period;
	double off_period;
};

SquareWave::SquareWave()
	:
	on_period(0),
	off_period(0)
{}


bool SquareWave::value(double time) const
{
	// start at time 0

	// ``pull back'' time to 1st total period:
	double t = std::fmod(time,  on_period+off_period );

	if ( t > on_period)
		return false;	// false if inactive
	else
		return true;	// true if secreting
}


double SquareWave::getOnPeriod() const
{
	return on_period;
}


double SquareWave::getOffPeriod() const
{
	return off_period;
}


void SquareWave::setOnPeriod(double on)
{
	on_period = on;
}


void SquareWave::setOffPeriod(double off)
{
	off_period = off;
}

}} // close namespace
#endif