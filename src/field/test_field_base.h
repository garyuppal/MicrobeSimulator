#ifndef MATH_TEST_FIELD_BASE_H
#define MATH_TEST_FIELD_BASE_H

#include "./field_base.h"

namespace Chemotaxis{ namespace Math{

class Test_FieldBase{
public:
	Test_FieldBase();
	void run();

private:

};

// IMPL

Test_FieldBase::Test_FieldBase()
{}

void
Test_FieldBase::run()
{
	FieldBase<2,int> field;
	field.printTest();

	int x = 5;
	int y = 4;

	assert((x+y)==9);

	std::cout << "got passed assertion" << std::endl;

}

}} // CLOSE NAMSPACE
#endif