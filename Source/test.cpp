// test.cpp
#include "IA.h"
#include <iostream>

void test_result(IIA::IA &ia)
{
  if (ia.get_result()->error)
	{
		std::cout << "error, ";
	}
  if (!ia.get_result()->solved)
	{
		std::cout << "failed to solve, ";
	}
  if (!ia.get_result()->constraints_satisfied)
	{
		std::cout << "failed to satisfy constraints, ";
	}
  if (!ia.get_result()->bounds_satisfied)
	{
		std::cout << "failed to satisfy bounds, ";
	}
	std::cout << std::endl;
}

void test_solution(IIA::IA &ia, std::vector<int> expected_solution)
{
	for (int i=0; i < expected_solution.size(); ++i)
	{
		const auto actual = ia.get_solution(i);
		const auto expected = expected_solution[i];
		if (actual != expected)
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

	IIA::IA ia;
  ia.get_result()->message_log = &std::cout;
  
  ia.get_result()->log_info = true;
  ia.get_result()->log_warning = true;
  ia.get_result()->log_error = true;
  ia.get_result()->log_debug = true;
  
	ia.resize(1,2);
	std::vector<int> cols = {0, 1};
	std::vector<int> vals = {1, -1};
	ia.set_row(0,cols,vals);
  ia.set_goal(0,1);
  ia.set_goal(1,4);
  ia.set_goal(1,4);
	ia.solve();

	test_result(ia);
	std::vector<int> expected_solution = {2,2};
	test_solution(ia,expected_solution);

  std::cout << "test_problem_1 finished" << std::endl;
}

int main(int argc, const char * argv[])
{

	test_problem_1();
	
}
