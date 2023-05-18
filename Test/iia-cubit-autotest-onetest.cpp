#include "IA.h"
#include "IncrementalIntervalAssignment.h"
#include "IAResultImplementation.h"
#include <iostream>

void setup_io(const IIA::IAResult *result);
void test_result(IIA::IA &ia);
void test_solution_autogen(IIA::IA &ia, std::vector<int> expected_solution, int dx=2)
{
  bool discrepancy = false;
  bool first_problem = true;
  int max_col = ia.used_cols();
  if (max_col > (int) expected_solution.size())
    max_col = (int) expected_solution.size();
  if (max_col > (int) ia.get_solution().size())
    max_col = (int) ia.get_solution().size();
  for (int c=0; c < max_col; ++c)
  {
    const auto actual = ia.get_solution(c);
    const auto expected = expected_solution[c];
    if (actual == expected)
      continue;
    
    if (actual+1 < ia.get_bound_hi(c) && actual+1 == expected)
      continue;
    if (actual-1 > ia.get_bound_lo(c) && actual-1 == expected)
      continue;
    
    // is it OK to be off by 2? yes if the variable is connected to a variable with a "2" coefficient somehow    
    /*
    int dx=2;
    const std::vector<int> *rows;
    ia.get_col(i, rows);
    for (auto r : *rows)
    {
      auto m = ia.get_row_col(r,i);
      if (abs(m)>1)
      {
        dx = std::max(dx,abs(m));
      }
    }
    */
    if (dx>1)
    {
      auto lo = std::max(ia.get_bound_lo(c), actual-dx );
      auto hi = std::min(ia.get_bound_hi(c), actual+dx );
      if (lo <= expected && expected <= hi)
        continue;
    }
    
    // off by more than allowed
    if (first_problem)
    {
      std::cout << "\n";
      first_problem = false;
    }
    std::cout << "solution x[" << c << "] = " << actual << ", NOT " << expected << " as expected." << std::endl;
    discrepancy=true;
  }
  if (discrepancy)
  {
    std::cout << "solution we found =";
    for (auto s : ia.get_solution())
      std::cout << " " << s;
    std::cout << std::endl;
    std::cout << "solution expected =";
    for (auto s : expected_solution)
      std::cout << " " << s;
    std::cout << std::endl;

  }
}


void iia_cubit_test_problem_1583198028742() 
{
  std::cout << "iia_cubit_test_problem_1583198028742   start: ";

  IIA::IA ia;
  setup_io(ia.get_result());
  ia.get_result()->log_debug=true;
  ia.get_result()->log_info=true;

  ia.resize(6, 30);
  // constraints
  {
    std::vector<int> cols = {0, 1, 2, 3, 4, 5, 19, 20, 21, 22, 23, 24};
    std::vector<int> vals = {1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1};
    ia.set_row(0,cols,vals);

    ia.set_rhs(0, 2);
  }
  {
    std::vector<int> cols = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 25};
    std::vector<int> vals = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1};
    ia.set_row(1,cols,vals);
  }
  {
    std::vector<int> cols = {1, 2, 3, 20, 21, 22, 26};
    std::vector<int> vals = {-1, -1, -1, 1, 1, 1, -1};
    ia.set_row(2,cols,vals);

    ia.set_rhs(2, 1);
  }
  {
    std::vector<int> cols = {2, 3, 21, 22, 27};
    std::vector<int> vals = {-1, -1, 1, 1, -1};
    ia.set_row(3,cols,vals);

    ia.set_rhs(3, 1);
  }
  {
    std::vector<int> cols = {3, 22, 28};
    std::vector<int> vals = {-1, 1, -1};
    ia.set_row(4,cols,vals);

    ia.set_rhs(4, 1);
  }
  {
    std::vector<int> cols = {0, 1, 2, 3, 4, 5, 19, 20, 21, 22, 23, 24, 29};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1};
    ia.set_row(5,cols,vals);

    ia.set_rhs(5, 1);
  }
  // bounds
  ia.set_bounds(26, 0, 0);
  ia.set_bounds(27, 0, 0);
  ia.set_bounds(28, 0, 0);
  ia.set_bounds(29, 0, 0);
  // goals
  ia.set_goal(0, 22.6103);
  ia.set_goal(1, 22.6103);
  ia.set_goal(2, 22.6103);
  ia.set_goal(3, 22.6103);
  ia.set_goal(4, 22.6103);
  ia.set_goal(5, 22.6103);
  ia.set_goal(6, 7.82423);
  ia.set_goal(7, 2.97172);
  ia.set_goal(8, 7.82423);
  ia.set_goal(9, 2.97172);
  ia.set_goal(10, 7.82423);
  ia.set_goal(11, 2.97172);
  ia.set_goal(12, 3.91212);
  ia.set_goal(13, 3.91212);
  ia.set_goal(14, 2.97172);
  ia.set_goal(15, 7.82423);
  ia.set_goal(16, 2.97172);
  ia.set_goal(17, 7.82423);
  ia.set_goal(18, 2.97172);
  ia.set_goal(19, 22.6103);
  ia.set_goal(20, 22.6103);
  ia.set_goal(21, 22.6103);
  ia.set_goal(22, 23.1103);
  ia.set_goal(23, 22.6103);
  ia.set_goal(24, 22.6103);
  ia.set_goal(25, 64.7757);
  ia.set_no_goal(26);
  ia.set_no_goal(27);
  ia.set_no_goal(28);
  ia.set_no_goal(29);

  ia.solve();

  // cubit version solved, did we?
  test_result(ia);
  // hand edited. for some reason cubit assigned the last dummy var to -3 instead of the correct value of 0
  std::vector<int> expected_solution = {23, 23, 23, 22, 23, 23, 8, 3, 8, 3, 8, 3, 4, 4, 3, 8, 3, 8, 3, 22, 23, 23, 23, 22, 22, 66, 0, 0, 0, 0}; 
  test_solution_autogen(ia,expected_solution);
  std::cout << "iia_cubit_test_problem_1583198028742   end" << std::endl;
}


void iia_cubit_test_problem_1583197819390()
{
  std::cout << "iia_cubit_test_problem_1583197819390   start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  ia.get_result()->log_debug=true;
  ia.get_result()->log_info=true;

  ia.resize(4, 4);
  // constraints
  {
    std::vector<int> cols = {1};
    std::vector<int> vals = {1};
    ia.set_row(0,cols,vals);
    
    ia.set_rhs(0, 14);
  }
  {
    std::vector<int> cols = {2};
    std::vector<int> vals = {1};
    ia.set_row(1,cols,vals);
    
    ia.set_rhs(1, -6);
  }
  {
    std::vector<int> cols = {3};
    std::vector<int> vals = {1};
    ia.set_row(2,cols,vals);
    
    ia.set_rhs(2, 22);
  }
  {
    std::vector<int> cols = {0};
    std::vector<int> vals = {2};
    ia.set_row(3,cols,vals);
    
    ia.set_rhs(3, 36);
  }
  // bounds
  ia.set_bound_lo(0, 2);
  ia.set_bound_lo(1, 0);
  ia.set_bound_lo(2, 0);
  ia.set_bound_lo(3, 0);
  // goals
  ia.set_no_goal(0);
  ia.set_no_goal(1);
  ia.set_no_goal(2);
  ia.set_no_goal(3);
  
  ia.solve();
  
  // cubit version solved, did we?
  test_result(ia);
  std::vector<int> expected_solution = {18, 0, 0, 0};
  test_solution_autogen(ia,expected_solution);
  std::cout << "iia_cubit_test_problem_1583197819390   end" << std::endl;
}

void iia_cubit_test_problem_1583197986982() 
{
  std::cout << "iia_cubit_test_problem_1583197986982   start: ";

  IIA::IA ia;
  setup_io(ia.get_result());
  ia.get_result()->log_debug=true;
  ia.get_result()->log_info=true;

  ia.resize(59, 183);
  // constraints
  {
    std::vector<int> cols = {0, 1, 2, 3, 4, 124};
    std::vector<int> vals = {-1, -1, -1, -1, -1, 2};
    ia.set_row(0,cols,vals);
  }
  {
    std::vector<int> cols = {1, 5, 6, 7, 125};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(1,cols,vals);
  }
  {
    std::vector<int> cols = {6, 8, 9, 10, 126};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(2,cols,vals);
  }
  {
    std::vector<int> cols = {9, 11, 12, 13, 127};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(3,cols,vals);
  }
  {
    std::vector<int> cols = {12, 14, 15, 16, 128};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(4,cols,vals);
  }
  {
    std::vector<int> cols = {15, 17, 18, 19, 129};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(5,cols,vals);
  }
  {
    std::vector<int> cols = {18, 20, 21, 22, 130};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(6,cols,vals);
  }
  {
    std::vector<int> cols = {21, 23, 24, 25, 131};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(7,cols,vals);
  }
  {
    std::vector<int> cols = {24, 26, 27, 28, 29, 132};
    std::vector<int> vals = {-1, -1, -1, -1, -1, 2};
    ia.set_row(8,cols,vals);
  }
  {
    std::vector<int> cols = {8, 11, 14, 17, 20, 30, 133};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, 2};
    ia.set_row(9,cols,vals);
  }
  {
    std::vector<int> cols = {3, 29, 31, 32, 33, 34, 134};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, 2};
    ia.set_row(10,cols,vals);
  }
  {
    std::vector<int> cols = {5, 23, 30, 35, 135};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(11,cols,vals);
  }
  {
    std::vector<int> cols = {4, 28, 31, 36, 136};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(12,cols,vals);
  }
  {
    std::vector<int> cols = {0, 27, 35, 36, 137};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(13,cols,vals);
  }
  {
    std::vector<int> cols = {37, 38, 39, 40, 138};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(14,cols,vals);
  }
  {
    std::vector<int> cols = {38, 41, 42, 43, 139};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(15,cols,vals);
  }
  {
    std::vector<int> cols = {42, 44, 45, 46, 140};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(16,cols,vals);
  }
  {
    std::vector<int> cols = {46, 47, 48, 49, 141};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(17,cols,vals);
  }
  {
    std::vector<int> cols = {37, 41, 45, 47, 50, 51, 52, 53, 54, 55, 56, 57, 142};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2};
    ia.set_row(18,cols,vals);
  }
  {
    std::vector<int> cols = {53, 58, 59, 60, 143};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(19,cols,vals);
  }
  {
    std::vector<int> cols = {52, 58, 61, 62, 144};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(20,cols,vals);
  }
  {
    std::vector<int> cols = {51, 61, 63, 64, 145};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(21,cols,vals);
  }
  {
    std::vector<int> cols = {16, 44, 49, 59, 65, 66, 67, 68, 69, 70, 146};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2};
    ia.set_row(22,cols,vals);
  }
  {
    std::vector<int> cols = {13, 62, 70, 71, 147};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(23,cols,vals);
  }
  {
    std::vector<int> cols = {72, 148};
    std::vector<int> vals = {-1, 2};
    ia.set_row(24,cols,vals);
  }
  {
    std::vector<int> cols = {2, 7, 10, 32, 63, 71, 73, 74, 75, 76, 77, 78, 79, 80, 149};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2};
    ia.set_row(25,cols,vals);
  }
  {
    std::vector<int> cols = {81, 82, 83, 84, 150};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(26,cols,vals);
  }
  {
    std::vector<int> cols = {84, 85, 86, 87, 88, 151};
    std::vector<int> vals = {-1, -1, -1, -1, -1, 2};
    ia.set_row(27,cols,vals);
  }
  {
    std::vector<int> cols = {82, 89, 90, 91, 92, 152};
    std::vector<int> vals = {-1, -1, -1, -1, -1, 2};
    ia.set_row(28,cols,vals);
  }
  {
    std::vector<int> cols = {83, 86, 92, 93, 153};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(29,cols,vals);
  }
  {
    std::vector<int> cols = {81, 85, 89, 94, 154};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(30,cols,vals);
  }
  {
    std::vector<int> cols = {95, 155};
    std::vector<int> vals = {-1, 2};
    ia.set_row(31,cols,vals);
  }
  {
    std::vector<int> cols = {19, 43, 65, 96, 156};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(32,cols,vals);
  }
  {
    std::vector<int> cols = {48, 57, 66, 97, 157};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(33,cols,vals);
  }
  {
    std::vector<int> cols = {55, 56, 67, 68, 97, 98, 158};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, 2};
    ia.set_row(34,cols,vals);
  }
  {
    std::vector<int> cols = {54, 60, 69, 98, 159};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(35,cols,vals);
  }
  {
    std::vector<int> cols = {99, 100, 101, 160};
    std::vector<int> vals = {-1, -2, -1, 2};
    ia.set_row(36,cols,vals);
  }
  {
    std::vector<int> cols = {101, 161};
    std::vector<int> vals = {-1, 2};
    ia.set_row(37,cols,vals);
  }
  {
    std::vector<int> cols = {72, 102, 103, 162};
    std::vector<int> vals = {-1, -2, -1, 2};
    ia.set_row(38,cols,vals);
  }
  {
    std::vector<int> cols = {103, 163};
    std::vector<int> vals = {-1, 2};
    ia.set_row(39,cols,vals);
  }
  {
    std::vector<int> cols = {95, 164};
    std::vector<int> vals = {-1, 2};
    ia.set_row(40,cols,vals);
  }
  {
    std::vector<int> cols = {104, 165};
    std::vector<int> vals = {-1, 2};
    ia.set_row(41,cols,vals);
  }
  {
    std::vector<int> cols = {79, 105, 106, 107, 166};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(42,cols,vals);
  }
  {
    std::vector<int> cols = {78, 106, 108, 109, 167};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(43,cols,vals);
  }
  {
    std::vector<int> cols = {77, 108, 110, 111, 168};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(44,cols,vals);
  }
  {
    std::vector<int> cols = {76, 110, 112, 113, 169};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(45,cols,vals);
  }
  {
    std::vector<int> cols = {114, 170};
    std::vector<int> vals = {-1, 2};
    ia.set_row(46,cols,vals);
  }
  {
    std::vector<int> cols = {40, 50, 64, 73, 115, 116, 171};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, 2};
    ia.set_row(47,cols,vals);
  }
  {
    std::vector<int> cols = {74, 87, 91, 93, 115, 117, 118, 119, 172};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, -1, -1, 2};
    ia.set_row(48,cols,vals);
  }
  {
    std::vector<int> cols = {114, 173};
    std::vector<int> vals = {-1, 2};
    ia.set_row(49,cols,vals);
  }
  {
    std::vector<int> cols = {120, 174};
    std::vector<int> vals = {-1, 2};
    ia.set_row(50,cols,vals);
  }
  {
    std::vector<int> cols = {120, 175};
    std::vector<int> vals = {-1, 2};
    ia.set_row(51,cols,vals);
  }
  {
    std::vector<int> cols = {121, 176};
    std::vector<int> vals = {-1, 2};
    ia.set_row(52,cols,vals);
  }
  {
    std::vector<int> cols = {121, 177};
    std::vector<int> vals = {-1, 2};
    ia.set_row(53,cols,vals);
  }
  {
    std::vector<int> cols = {104, 178};
    std::vector<int> vals = {-1, 2};
    ia.set_row(54,cols,vals);
  }
  {
    std::vector<int> cols = {99, 179};
    std::vector<int> vals = {-1, 2};
    ia.set_row(55,cols,vals);
  }
  {
    std::vector<int> cols = {22, 25, 26, 34, 39, 96, 107, 109, 111, 113, 116, 118, 122, 123, 180};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2};
    ia.set_row(56,cols,vals);
  }
  {
    std::vector<int> cols = {33, 80, 105, 122, 181};
    std::vector<int> vals = {-1, -1, -1, -1, 2};
    ia.set_row(57,cols,vals);
  }
  {
    std::vector<int> cols = {75, 88, 90, 94, 112, 117, 119, 123, 182};
    std::vector<int> vals = {-1, -1, -1, -1, -1, -1, -1, -1, 2};
    ia.set_row(58,cols,vals);
  }
  // bounds
  ia.set_bound_lo(72, 3);
  ia.set_bound_lo(95, 3);
  ia.set_bound_lo(99, 3);
  ia.set_bound_lo(101, 3);
  ia.set_bound_lo(103, 3);
  ia.set_bound_lo(104, 3);
  ia.set_bound_lo(114, 3);
  ia.set_bound_lo(120, 3);
  ia.set_bound_lo(121, 3);
  ia.set_bound_lo(124, 2);
  ia.set_bound_lo(125, 2);
  ia.set_bound_lo(126, 2);
  ia.set_bound_lo(127, 2);
  ia.set_bound_lo(128, 2);
  ia.set_bound_lo(129, 2);
  ia.set_bound_lo(130, 2);
  ia.set_bound_lo(131, 2);
  ia.set_bound_lo(132, 2);
  ia.set_bound_lo(133, 2);
  ia.set_bound_lo(134, 2);
  ia.set_bound_lo(135, 2);
  ia.set_bound_lo(136, 2);
  ia.set_bound_lo(137, 2);
  ia.set_bound_lo(138, 3);
  ia.set_bound_lo(139, 2);
  ia.set_bound_lo(140, 2);
  ia.set_bound_lo(141, 2);
  ia.set_bound_lo(142, 2);
  ia.set_bound_lo(143, 2);
  ia.set_bound_lo(144, 2);
  ia.set_bound_lo(145, 3);
  ia.set_bound_lo(146, 2);
  ia.set_bound_lo(147, 2);
  ia.set_bound_lo(148, 2);
  ia.set_bound_lo(149, 2);
  ia.set_bound_lo(150, 2);
  ia.set_bound_lo(151, 2);
  ia.set_bound_lo(152, 2);
  ia.set_bound_lo(153, 2);
  ia.set_bound_lo(154, 3);
  ia.set_bound_lo(155, 2);
  ia.set_bound_lo(156, 2);
  ia.set_bound_lo(157, 2);
  ia.set_bound_lo(158, 2);
  ia.set_bound_lo(159, 2);
  ia.set_bound_lo(160, 2);
  ia.set_bound_lo(161, 2);
  ia.set_bound_lo(162, 2);
  ia.set_bound_lo(163, 2);
  ia.set_bound_lo(164, 2);
  ia.set_bound_lo(165, 2);
  ia.set_bound_lo(166, 2);
  ia.set_bound_lo(167, 2);
  ia.set_bound_lo(168, 2);
  ia.set_bound_lo(169, 2);
  ia.set_bound_lo(170, 2);
  ia.set_bound_lo(171, 2);
  ia.set_bound_lo(172, 2);
  ia.set_bound_lo(173, 2);
  ia.set_bound_lo(174, 2);
  ia.set_bound_lo(175, 2);
  ia.set_bound_lo(176, 2);
  ia.set_bound_lo(177, 2);
  ia.set_bound_lo(178, 2);
  ia.set_bound_lo(179, 2);
  ia.set_bound_lo(180, 2);
  ia.set_bound_lo(181, 2);
  ia.set_bound_lo(182, 2);
  // goals
  ia.set_goal(0, 0.621693);
  ia.set_goal(1, 0.101369);
  ia.set_goal(2, 0.661972);
  ia.set_goal(3, 0.0167756);
  ia.set_goal(4, 0.093824);
  ia.set_goal(5, 0.195669);
  ia.set_goal(6, 0.101369);
  ia.set_goal(7, 0.260891);
  ia.set_goal(8, 0.151023);
  ia.set_goal(9, 0.101369);
  ia.set_goal(10, 0.151023);
  ia.set_goal(11, 0.202738);
  ia.set_goal(12, 0.101369);
  ia.set_goal(13, 0.304106);
  ia.set_goal(14, 0.903468);
  ia.set_goal(15, 0.101369);
  ia.set_goal(16, 0.903468);
  ia.set_goal(17, 0.202738);
  ia.set_goal(18, 0.101369);
  ia.set_goal(19, 0.304106);
  ia.set_goal(20, 0.151023);
  ia.set_goal(21, 0.101369);
  ia.set_goal(22, 0.151023);
  ia.set_goal(23, 0.195669);
  ia.set_goal(24, 0.101369);
  ia.set_goal(25, 0.260891);
  ia.set_goal(26, 0.661972);
  ia.set_goal(27, 0.621693);
  ia.set_goal(28, 0.093824);
  ia.set_goal(29, 0.0167756);
  ia.set_goal(30, 1.1616);
  ia.set_goal(31, 1.2875);
  ia.set_goal(32, 1.37446);
  ia.set_goal(33, 1.29067);
  ia.set_goal(34, 1.37446);
  ia.set_goal(35, 1.1616);
  ia.set_goal(36, 1.1616);
  ia.set_goal(37, 3.34502);
  ia.set_goal(38, 0.101369);
  ia.set_goal(39, 3.15718);
  ia.set_goal(40, 0.207529);
  ia.set_goal(41, 0.208293);
  ia.set_goal(42, 0.103018);
  ia.set_goal(43, 0.304348);
  ia.set_goal(44, 0.129067);
  ia.set_goal(45, 0.129067);
  ia.set_goal(46, 0.103018);
  ia.set_goal(47, 1.8317);
  ia.set_goal(48, 1.92773);
  ia.set_goal(49, 0.13754);
  ia.set_goal(50, 1.1616);
  ia.set_goal(51, 3.34503);
  ia.set_goal(52, 0.208293);
  ia.set_goal(53, 0.129067);
  ia.set_goal(54, 1.83171);
  ia.set_goal(55, 0.322666);
  ia.set_goal(56, 2.12649e-05);
  ia.set_goal(57, 0.322666);
  ia.set_goal(58, 0.103018);
  ia.set_goal(59, 0.129067);
  ia.set_goal(60, 0.103018);
  ia.set_goal(61, 0.101369);
  ia.set_goal(62, 0.304348);
  ia.set_goal(63, 3.15718);
  ia.set_goal(64, 0.207529);
  ia.set_goal(65, 0.329971);
  ia.set_goal(66, 0.322666);
  ia.set_goal(67, 2.12649e-05);
  ia.set_goal(68, 0.322666);
  ia.set_goal(69, 0.13754);
  ia.set_goal(70, 0.329971);
  ia.set_goal(71, 0.320515);
  ia.set_goal(72, 5.37245);
  ia.set_goal(73, 0.902447);
  ia.set_goal(74, 0.78396);
  ia.set_goal(75, 0.774401);
  ia.set_goal(76, 0.0420516);
  ia.set_goal(77, 0.17761);
  ia.set_goal(78, 0.960115);
  ia.set_goal(79, 0.170899);
  ia.set_goal(80, 0.00532931);
  ia.set_goal(81, 0.391433);
  ia.set_goal(82, 0.724657);
  ia.set_goal(83, 1.87054);
  ia.set_goal(84, 0.724655);
  ia.set_goal(85, 0.000842076);
  ia.set_goal(86, 0.0137935);
  ia.set_goal(87, 0.4936);
  ia.set_goal(88, 0.491451);
  ia.set_goal(89, 0.000841721);
  ia.set_goal(90, 0.491451);
  ia.set_goal(91, 0.4936);
  ia.set_goal(92, 0.0137922);
  ia.set_goal(93, 1.91628);
  ia.set_goal(94, 0.620118);
  ia.set_goal(95, 0.608213);
  ia.set_goal(96, 0.320515);
  ia.set_goal(97, 1.92773);
  ia.set_goal(98, 1.92773);
  ia.set_goal(99, 5.37243);
  ia.set_goal(100, 0.516976);
  ia.set_goal(101, 5.20245);
  ia.set_goal(102, 0.516976);
  ia.set_goal(103, 5.20245);
  ia.set_goal(104, 0.608213);
  ia.set_goal(105, 1.29067);
  ia.set_goal(106, 1.29067);
  ia.set_goal(107, 0.170899);
  ia.set_goal(108, 1.29067);
  ia.set_goal(109, 0.960115);
  ia.set_goal(110, 1.29067);
  ia.set_goal(111, 0.17761);
  ia.set_goal(112, 1.29067);
  ia.set_goal(113, 0.0420516);
  ia.set_goal(114, 3.29299);
  ia.set_goal(115, 2.02738);
  ia.set_goal(116, 0.902447);
  ia.set_goal(117, 0.0555492);
  ia.set_goal(118, 0.78396);
  ia.set_goal(119, 0.0555492);
  ia.set_goal(120, 2.43285);
  ia.set_goal(121, 2.43285);
  ia.set_goal(122, 0.00532931);
  ia.set_goal(123, 0.774401);
  ia.set_no_goal(124);
  ia.set_no_goal(125);
  ia.set_no_goal(126);
  ia.set_no_goal(127);
  ia.set_no_goal(128);
  ia.set_no_goal(129);
  ia.set_no_goal(130);
  ia.set_no_goal(131);
  ia.set_no_goal(132);
  ia.set_no_goal(133);
  ia.set_no_goal(134);
  ia.set_no_goal(135);
  ia.set_no_goal(136);
  ia.set_no_goal(137);
  ia.set_no_goal(138);
  ia.set_no_goal(139);
  ia.set_no_goal(140);
  ia.set_no_goal(141);
  ia.set_no_goal(142);
  ia.set_no_goal(143);
  ia.set_no_goal(144);
  ia.set_no_goal(145);
  ia.set_no_goal(146);
  ia.set_no_goal(147);
  ia.set_no_goal(148);
  ia.set_no_goal(149);
  ia.set_no_goal(150);
  ia.set_no_goal(151);
  ia.set_no_goal(152);
  ia.set_no_goal(153);
  ia.set_no_goal(154);
  ia.set_no_goal(155);
  ia.set_no_goal(156);
  ia.set_no_goal(157);
  ia.set_no_goal(158);
  ia.set_no_goal(159);
  ia.set_no_goal(160);
  ia.set_no_goal(161);
  ia.set_no_goal(162);
  ia.set_no_goal(163);
  ia.set_no_goal(164);
  ia.set_no_goal(165);
  ia.set_no_goal(166);
  ia.set_no_goal(167);
  ia.set_no_goal(168);
  ia.set_no_goal(169);
  ia.set_no_goal(170);
  ia.set_no_goal(171);
  ia.set_no_goal(172);
  ia.set_no_goal(173);
  ia.set_no_goal(174);
  ia.set_no_goal(175);
  ia.set_no_goal(176);
  ia.set_no_goal(177);
  ia.set_no_goal(178);
  ia.set_no_goal(179);
  ia.set_no_goal(180);
  ia.set_no_goal(181);
  ia.set_no_goal(182);

  ia.solve();

  // cubit version solved, did we?
  test_result(ia);
  std::vector<int> expected_solution = {1, 1, 2, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 4, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 3, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 2, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 4, 1, 2, 2, 6, 1, 6, 1, 6, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 4, 4, 1, 1, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 2, 2, 2, 5, 2, 2, 3, 10, 3, 2, 4, 6, 3, 3, 9, 4, 3, 3, 4, 3, 2, 2, 3, 4, 3, 7, 3, 7, 3, 2, 2, 2, 2, 2, 2, 2, 4, 6, 2, 2, 2, 2, 2, 2, 3, 9, 2, 5};
  test_solution_autogen(ia,expected_solution);
  std::cout << "iia_cubit_test_problem_1583197986982   end" << std::endl;
}

void iia_cubit_autotest()
{
  iia_cubit_test_problem_1583198028742(); // the one with a warning
  iia_cubit_test_problem_1583197819390(); // the one a bad solution 2, instead of 18
  // iia_cubit_test_problem_1583197986982();
}
