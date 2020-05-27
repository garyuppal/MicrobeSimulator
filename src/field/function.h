#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

#include "../vector/vect.h"

namespace Math{

template<int dim, typename Number>
class Function{
public:
	Function() {}
	virtual ~Function() {}

	virtual double value(const Vect<dim,Number>& ) const=0;
};