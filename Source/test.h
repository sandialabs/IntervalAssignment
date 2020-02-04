// test.h
#include "IA.h"
#include <iostream>

void test_result(IA &ia)
{
	if (ia.get_result().error)
	{
		std::cout << "error, ";
	}
	if (!ia.get_result().solved)
	{
		std::cout << "failed to solve, ";
	}
 	if (!ia.get_result().constraints_satisfied())
	{
		std::cout << "failed to satisfy constraints, ";
	}
 	if (!ia.get_result().bounds_satisfied())
	{
		std::cout << "failed to satisfy bounds, ";
	}
	std::cout << "\nlog\n";
	for (auto &message : log)
	{
		std::cout << message << "\n";
	}
	std::cout << std::endl;
}

void test_solution(IA &ia, std::vector<int> expected_solution)
{
	for (int i=0; i < expected_solution.size(); ++i)
	{
		auto actual = ia.get_solution(i);
		auto expected = expected_solution[i];
		if (actual != expected_solution)
		{
			std::cout << "solution x[" << i << "] = " << actual << ", NOT " << expected << " as expected." << std::endl;
		}
	}
}

void test_problem_1()
{
	// x0 - x1 = 0
	// g[0]=1
	// g[1]=4
	// expect solution [2,2]
	std::cout << "test_problem_1 started: ";

	IA ia;
	ia.resize(1,2);
	std::vector<int> cols = {0, 1};
	std::vector<int> vals = {1, -1};
	ia.set_row(0,cols,vals);
	ia.set_goal(0,1);
	ia.set_goal(1,4);
	ia.solve();

	test_result(ia);
	std::vector<int> expected_solution = {2,2};
	test_solution(ia,expected_solution);

}

main
{

	test_problem_1();
	
}