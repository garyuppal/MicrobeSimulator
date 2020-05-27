#ifndef MATH_TEST_FDM_FIELD_H
#define MATH_TEST_FDM_FIELD_H

#include "./fdm_field.h"

namespace Chemotaxis{ namespace Math{

class Test_FDMField{
public:
	Test_FDMField();
	void run();

private:
	void test_intialization_2d();
	void test_intialization_3d();
};

// IMPL
Test_FDMField::Test_FDMField()
{}

void
Test_FDMField::run()
{
	test_intialization_2d();
	test_intialization_3d();
}


// TESTS:
void 
Test_FDMField::test_intialization_2d()
{
	Vect<2> lower(0,0);
	Vect<2> upper(10,10);

	std::array<unsigned int, 2> disc = {10,10};

	FDM_Field<2> field(lower, upper, 
		disc, 3.4);

	field.print(std::cout);

	std::cout << "--------------------------------------" << std::endl;
}

void 
Test_FDMField::test_intialization_3d()
{
	Vect<3> lower(0,0,0);
	Vect<3> upper(10,10,10);

	std::array<unsigned int, 3> disc = {5,5,3};

	FDM_Field<3> field(lower, upper, 
		disc, -1.2);

	field.print(std::cout);

	std::cout << "--------------------------------------" << std::endl;
}

}} // close namespace
#endif