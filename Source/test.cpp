// test.cpp
#include "IA.h"
#include "IncrementalIntervalAssignment.h"
#include "IAResultImplementation.h"
#include <iostream>

// tests
// HNF, RREF on a few small matrices
// opposite sides equal
// sum-even
// swept volume: sum-even top and bottom, mapped sides
// submap

// tests todo:

// radish test
// triangle prism where lower bounds on sum-even variables matter
// triangle primitive inequality constraints
// triangle primitive mixed with other constraints
// some infeasible prolems, like the slanted L
// u sub-map, dynamic constraints

void setup_io(const IIA::IAResult *result)
{
  result->message_log = &std::cout;
  result->log_info = false;
  result->log_warning = true;
  result->log_error = true;
  result->log_debug = false;
}

void test_result(IIA::IA &ia)
{
  bool something_printed = false;
  if (ia.get_result()->error)
  {
    std::cout << "error, ";
    something_printed = true;
  }
  if (!ia.get_result()->solved)
  {
    std::cout << "failed to solve, ";
    something_printed = true;
  }
  if (!ia.get_result()->constraints_satisfied)
  {
    std::cout << "failed to satisfy constraints, ";
    something_printed = true;
  }
  if (!ia.get_result()->bounds_satisfied)
  {
    std::cout << "failed to satisfy bounds, ";
    something_printed = true;
  }
  if (something_printed)
  {
    std::cout << std::endl;
  }
}

void test_solution(IIA::IA &ia, std::vector<int> expected_solution)
{
  bool discrepancy = false;
  for (int i=0; i < expected_solution.size(); ++i)
  {
    const auto actual = ia.get_solution(i);
    const auto expected = expected_solution[i];
    if (actual != expected)
    {
      std::cout << "solution x[" << i << "] = " << actual << ", NOT " << expected << " as expected." << std::endl;
      discrepancy=true;
    }
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

// solution is symmetric, x 2,4,6,8,10, are equivalent, 3 should be 1, 3 should be 2
void shuffle_solution(IIA::IA &ia, std::vector<int> &expected_solution, const std::vector<int> &equivalent_indices)
{
  std::vector<int> new_solution = expected_solution;
  std::multiset<int> sol_set;
  for (auto c : equivalent_indices)
  {
    sol_set.insert(expected_solution[c]);
  }
  for (auto c : equivalent_indices)
  {
    auto s = ia.get_solution(c);
    auto i = sol_set.find(s);
    if (i==sol_set.end())
    {
      // error, nobody left is supposed to have that value
      return;
    }
    sol_set.erase(i);
    new_solution[c]=s;
  }
  // if got to here, all is well
  expected_solution.swap(new_solution);
}


void test_problem_A()
{
  // test solving a single equality constraint with goals
  // x0 - x1 = 0
  // g[0]=1
  // g[1]=4
  // expect solution [2,2]
  std::cout << "test_problem_A start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
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
  
  std::cout << "test_problem_A end." << std::endl;
}


void test_problem_AA()
{
  // dynamic sizing
  // test solving a single equality constraint with goals
  // x0 - x1 = 0
  // g[0]=1
  // g[1]=4
  // expect solution [2,2]
  std::cout << "test_problem_AA start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
  // ia.resize(1,2);
  std::vector<int> vals = {1, -1};
  std::vector<double> goals = {1, 4};
  int r = ia.next_row();
  for (int i = 0; i < vals.size(); ++i)
  {
    int c = ia.next_col();
    ia.set_row_col(r,c,vals[i]);
    ia.set_goal(c,goals[i]);
  }

  ia.solve();
  
  test_result(ia);
  std::vector<int> expected_solution = {2,2};
  test_solution(ia,expected_solution);
  
  std::cout << "test_problem_AA end." << std::endl;
}

void test_problem_B()
{
  // test solving a single sum-even with goals
  // x0 + x1 + x2 - 2x3= 0
  // g[0]=1
  // g[1]=3
  // g[2]=5
  // g[3]=0, a dummy sum-even variable
  // expect solution [1,3,6,5]
  std::cout << "test_problem_B start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
  ia.resize(1,4); // 1,2 causes error
  std::vector<int> cols = {0, 1, 2, 3};
  std::vector<int> vals = {1, 1, 1, -2};
  ia.set_row(0,cols,vals);
  ia.set_goal(0,1);
  ia.set_goal(1,3);
  ia.set_goal(2,5);
  ia.set_no_goal(3);
  ia.set_bound_lo( 3, 2 );
  
  ia.solve();
  
  test_result(ia);
  std::vector<int> expected_solution = {1,3,6,5};
  test_solution(ia,expected_solution);
  
  std::cout << "test_problem_B end." << std::endl;
}

void test_problem_BB()
{
  // test solving a single sum-even with goals, but using IIA::EVEN
  // x0 + x1 + x2 = Even
  // g[0]=1
  // g[1]=3
  // g[2]=5
  // expect solution [1,3,6]
  std::cout << "test_problem_BB start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
  ia.resize(1,3);
  std::vector<int> cols = {0, 1, 2};
  std::vector<int> vals = {1, 1, 1};
  ia.set_row(0,cols,vals);
  ia.set_goal(0,1);
  ia.set_goal(1,3);
  ia.set_goal(2,5);
  ia.set_constraint(0, IIA::EVEN);
  
  ia.solve();
  
  test_result(ia);
  std::vector<int> expected_solution = {1,3,6,5};
  test_solution(ia,expected_solution);
  
  std::cout << "test_problem_BB end." << std::endl;
}


void test_problem_BC()
{
  // dynamic sizing
  // test solving a single sum-even with goals, but using IIA::EVEN
  // x0 + x1 + x2 = Even
  // g[0]=1
  // g[1]=3
  // g[2]=5
  // expect solution [1,3,6]
  std::cout << "test_problem_BC start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
  std::vector<int> cols = {0, 1, 2};
  std::vector<int> vals = {1, 1, 1};
  std::vector<double> goals = {1, 3, 5};
  int r = ia.next_row();
  ia.set_constraint(r, IIA::EVEN);
  for (int i = 0; i < vals.size(); ++i)
  {
    int c = ia.next_col();
    ia.set_row_col(r,c,vals[i]);
    ia.set_goal(c,goals[i]);
  }

  ia.solve();
  
  test_result(ia);
  std::vector<int> expected_solution = {1,3,6,5};
  test_solution(ia,expected_solution);
  
  std::cout << "test_problem_BC end." << std::endl;
}



void test_problem_C()
{
  // test a cube with top and bottom paved (sum-even) and sides mapped (equal)
  //  1 -1
  //     1 -1
  //        1 -1
  // -1        1
  //             1      -1
  //               1      -1
  //                 1      -1
  //                   1      -1
  //             1 1 1 1        -2
  //                     1 1 1 1  -2
  // just set goals to some arbitrary floating point numbers
  // g[0]=1
  // g[1]=3
  // g[2]=5
  // g[3]=0, a dummy sum-even variable
  // expect solution [1,3,6,5]
  std::cout << "test_problem_C start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
  ia.resize(10,14);
  // linked
  {
    std::vector<int> cols = {0,  1};
    std::vector<int> vals = {1, -1};
    ia.set_row(0,cols,vals);
  }
  {
    std::vector<int> cols = {1,  2};
    std::vector<int> vals = {1, -1};
    ia.set_row(1,cols,vals);
  }
  {
    std::vector<int> cols = {2,  3};
    std::vector<int> vals = {1, -1};
    ia.set_row(2,cols,vals);
  }
  {
    std::vector<int> cols = {3,  0};
    std::vector<int> vals = {1, -1};
    ia.set_row(3,cols,vals);
  }
  
  // up/down
  {
    std::vector<int> cols = {4,  8};
    std::vector<int> vals = {1, -1};
    ia.set_row(4,cols,vals);
  }
  {
    std::vector<int> cols = {5,  9};
    std::vector<int> vals = {1, -1};
    ia.set_row(5,cols,vals);
  }
  {
    std::vector<int> cols = {6, 10};
    std::vector<int> vals = {1, -1};
    ia.set_row(6,cols,vals);
  }
  {
    std::vector<int> cols = {7, 11};
    std::vector<int> vals = {1, -1};
    ia.set_row(7,cols,vals);
  }
  
  // sum-even tops
  {
    std::vector<int> cols = {4, 5, 6, 7, 12};
    std::vector<int> vals = {1, 1, 1, 1, -2};
    ia.set_row(8,cols,vals);
    ia.set_bound_lo( 12, 2 );
  }
  {
    std::vector<int> cols = {8, 9,10,11, 13};
    std::vector<int> vals = {1, 1, 1, 1, -2};
    ia.set_row(9,cols,vals);
    ia.set_bound_lo( 13, 2 );
  }
  
  ia.set_goal (0,  0.1);
  ia.set_goal (1,  1.2);
  ia.set_goal (2,  2.1);
  ia.set_goal (3,  3.2);
  ia.set_goal (4,  3.2);
  ia.set_goal (5,  5.2);
  ia.set_goal (6,  6.2);
  ia.set_goal (7,  7.2);
  ia.set_goal (8,  8.2);
  ia.set_goal (9,  9.4);
  ia.set_goal(10, 10.1);
  ia.set_goal(11, 11.3);
  ia.set_no_goal(12);
  ia.set_no_goal(13);
  
  ia.solve();
  
  test_result(ia);
  std::vector<int> expected_solution = {1,1,1,1,5,7,8,10,5,7,8,10,15,15};
  test_solution(ia,expected_solution);
  
  std::cout << "test_problem_C end." << std::endl;
}


void test_problem_D()
{
  // submap test
  // 0 opposite 2,4,6,8,10
  // alternating 1,3,5,7,9,11
  std::cout << "test_problem_D start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
  ia.resize(2,12);
  {
    std::vector<int> cols = { 0, 2, 4, 6, 8, 10};
    std::vector<int> vals = {-1, 1, 1, 1, 1,  1};
    ia.set_row(0, cols, vals);
  }
  {
    std::vector<int> cols = { 1, 3,  5, 7,  9, 11};
    std::vector<int> vals = {-1, 1, -1, 1, -1,  1};
    ia.set_row(1, cols, vals);
  }
  ia.set_goal(0,8);
  ia.set_goal(2,1);
  ia.set_goal(4,1);
  ia.set_goal(6,1);
  ia.set_goal(8,1);
  ia.set_goal(10,1);
  
  ia.set_goal(1,4);
  ia.set_goal(3,2);
  ia.set_goal(5,3);
  ia.set_goal(7,4);
  ia.set_goal(9,5);
  ia.set_goal(11,2);
  
  ia.solve();
  
  test_result(ia);
  
  // solution is symmetric, x 2,4,6,8,10, are equivalent, but they are all expected to be 1
  std::vector<int> expected_solution = {5,3,1,2,1,3,1,6,1,4,1,2};
  test_solution(ia,expected_solution);
  
  std::cout << "test_problem_D end." << std::endl;
}

void test_problem_E()
{
  // submap test
  // 0 opposite 2,4,6,8,10
  // alternating 1,3,5,7,9,11
  std::cout << "test_problem_E start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
  ia.resize(2,12);
  {
    std::vector<int> cols = { 0, 2, 4, 6, 8, 10};
    std::vector<int> vals = {-1, 1, 1, 1, 1,  1};
    ia.set_row(0, cols, vals);
  }
  {
    std::vector<int> cols = { 1, 3,  5, 7,  9, 11};
    std::vector<int> vals = {-1, 1, -1, 1, -1,  1};
    ia.set_row(1, cols, vals);
  }
  ia.set_goal(0,16);
  ia.set_goal(2,1);
  ia.set_goal(4,1);
  ia.set_goal(6,1);
  ia.set_goal(8,1);
  ia.set_goal(10,1);
  
  ia.set_goal(1,4);
  ia.set_goal(3,2);
  ia.set_goal(5,3);
  ia.set_goal(7,4);
  ia.set_goal(9,5);
  ia.set_goal(11,2);
  
  ia.solve();
  
  test_result(ia);
  
  std::vector<int> expected_solution = {8,3,1,2,2,3,1,6,2,4,2,2};
  // solution is symmetric, x 2,4,6,8,10, are equivalent, 3 should be 1, 3 should be 2
  shuffle_solution(ia,expected_solution, {2,4,6,8,10} );
  test_solution(ia,expected_solution);
  
  std::cout << "test_problem_E end." << std::endl;
}


void test_problem_F()
{
  // U submap test with resolving.
  // 0 opposite 2,4,6
  // alternating 1,3,5,7
  std::cout << "test_problem_F start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
  ia.resize(2,8);
  {
    std::vector<int> cols = { 0, 2, 4, 6};
    std::vector<int> vals = {-1, 1, 1, 1};
    ia.set_row(0, cols, vals);
  }
  {
    std::vector<int> cols = { 1, 3,  5, 7};
    std::vector<int> vals = {-1, 1, -1, 1};
    ia.set_row(1, cols, vals);
  }
  ia.set_goal(0,6);
  ia.set_goal(2,2);
  ia.set_goal(4,2);
  ia.set_goal(6,2);
  
  ia.set_goal(1,4);
  ia.set_goal(3,7);
  ia.set_goal(5,8);
  ia.set_goal(7,4);
  
  
  ia.solve();

  test_result(ia);
  
  if (ia.is_solved() && !ia.get_result()->error)
  {
    {
      // std::vector<int> expected_solution = {6,4,2,8,2,8,2,4};
      std::vector<int> expected_solution = {6,4,2,7,2,7,2,4};
      test_solution(ia,expected_solution);
    }
    
    // expect this to occur once
    // if x3 >= x1
    if (ia.get_solution(3) >= ia.get_solution(1))
    {
      // add constraint x3 < x1  <==>  -x1 + x3 < 0   <==>   -x1 + x3 <= -1
      auto r = ia.new_row();
      std::vector<int> cols = { 1, 3};
      std::vector<int> vals = {-1, 1};
      ia.set_row(r, cols, vals);
      ia.set_rhs(r,-1);
      ia.set_constraint(r, IIA::LE);
      
      ia.solve(); // re-solve
      
      test_result(ia);

      std::vector<int> expected_solution = {6,6,2,5,2,5,2,6};
      test_solution(ia,expected_solution);
    }
  }
  
  std::cout << "test_problem_F end." << std::endl;
}

void test_problem_G()
{
  // U submap test with resolving. multiple legs
  // 0 opposite 2,4,6,8,10,12,14
  // alternating 1,3,5,7,9,11,13,15
  std::cout << "test_problem_G start: ";
  
  IIA::IA ia;
  setup_io(ia.get_result());
  
  // ia.get_result()->log_info=true;
  // ia.get_result()->log_debug=true;
  
  ia.resize(2,16);
  {
    std::vector<int> cols = { 0, 2, 4, 6, 8, 10, 12, 14};
    std::vector<int> vals = {-1, 1, 1, 1, 1,  1,  1,  1};
    ia.set_row(0, cols, vals);
  }
  std::vector<int> cren_cols = { 1, 3,  5, 7,  9, 11, 13, 15};
  std::vector<int> cren_vals = {-1, 1, -1, 1, -1,  1, -1,  1};
  {
    std::vector <int> cols = cren_cols;
    std::vector <int> vals = cren_vals;
    ia.set_row(1, cols, vals);
  }

  ia.set_goal(0,6);
  ia.set_goal(2,2);
  ia.set_goal(4,2);
  ia.set_goal(6,2);
  ia.set_goal(8,2);
  ia.set_goal(10,2);
  ia.set_goal(12,2);
  ia.set_goal(14,2);

  // increasing
  ia.set_goal(1,  2);
  ia.set_goal(3,  4);
  ia.set_goal(5,  6);
  ia.set_goal(7,  8);
  ia.set_goal(9, 10);
  ia.set_goal(11,16);
  ia.set_goal(13,22);
  ia.set_goal(15, 1);
  
  int iter=0;
  bool resolve=true;
  while (resolve)
  {
    ia.solve();
    resolve=false;
    
    test_result(ia);
    
    if (iter++ < 5 && ia.is_solved() && !ia.get_result()->error)
    {
      int sum = 0;
      for (size_t i=0; i<cren_cols.size()-1; ++i)
      {
        sum += cren_vals[i]*ia.get_solution(cren_cols[i]);
        if (sum >= 0)
        {
          auto r = ia.new_row();
          int v = 1;
          for (size_t j = 0; j <= i; ++j)
          {
            ia.set_row_col(r, cren_cols[j], v);
            v = -v;
          }
          ia.set_rhs(r,1);
          ia.set_constraint(r, IIA::GE);
          resolve=true;
          break;
        }
      }
    }
  }

  std::vector<int> expected_solution = {12, 3, 1, 2, 1, 7, 2, 7, 2, 13, 2, 13, 2, 4, 2, 5};
  shuffle_solution(ia,expected_solution, {2,4,6,8,10,12,14} );
  test_solution(ia,expected_solution);
  
  std::cout << "test_problem_G end." << std::endl;
}

namespace IIA_Internal
{
  // tester class is derived so we can test protected members
  class IIATester : public IncrementalIntervalAssignment
  {
  public:
    IIATester(IAResultImplementation *result_ptr) : IncrementalIntervalAssignment(result_ptr) {setup_io(result_ptr);}
    
    bool test_hnf_solver(std::vector<int> &expected_solution, std::vector<int> *nullspace_vector=nullptr, bool expect_fail=false);
    bool test_rref_constraints();
    bool test_rref_improve();
    
    // regression tests
    static bool test_problem0();
    static bool test_problem1();
    static bool test_problem2();
    static bool test_problem3();
  };
  
  bool IIATester::test_rref_constraints()
  {
    // test putting the matrix into rref form, choosing pivots in the context of setting up to solve constraints.
    std::vector<int> rref_col_order;
    bool rref_OK = rref_constraints(rref_col_order);
    if (!rref_OK)
    {
      std::cout << "ERROR: Failed to put matrix into RREF form for constraints!\n";
    }
    
    std::vector<int> cols_dep, cols_ind, rows_dep;
    categorize_vars(rref_col_order, cols_dep, cols_ind, rows_dep);
    MatrixSparseInt M(result);
    int MaxMrow(0);
    bool nullspace_OK = create_nullspace(cols_dep, cols_ind, rows_dep, M, MaxMrow);
    // M.print_matrix("nullspace");
    if (!nullspace_OK)
    {
      std::cout << "ERROR: Failed to find nullspace of RREF for constraints!\n";
    }
    
    return (rref_OK && nullspace_OK);
  }
  
  bool IIATester::test_rref_improve()
  {
    // test putting the matrix into rref form, choosing pivots in the context of setting up to solve constraints.
    std::vector<int> rref_col_order;
    bool rref_OK = rref_improve(rref_col_order);
    if (!rref_OK)
    {
      std::cout << "ERROR: Failed to put matrix into RREF form for improve!\n";
    }
    
    std::vector<int> cols_dep, cols_ind, rows_dep;
    categorize_vars(rref_col_order, cols_dep, cols_ind, rows_dep);
    MatrixSparseInt M(result);
    int MaxMrow(0);
    bool nullspace_OK = create_nullspace(cols_dep, cols_ind, rows_dep, M, MaxMrow);
    // M.print_matrix("nullspace");
    if (!nullspace_OK)
    {
      std::cout << "ERROR: Failed to find nullspace of RREF for improve!\n";
    }
    
    return (rref_OK && nullspace_OK);
  }
  
  bool IIATester::test_hnf_solver(std::vector<int> &expected_solution, std::vector<int> *nullspace_vector, bool expect_fail)
  {
    IIA_Internal::MatrixSparseInt B(result), U(result);
    std::vector<int> hnf_col_order;
    auto OK = HNF(B,U,hnf_col_order);
    if (!OK)
    {
      std::cout << "ERROR: HNF failed\n";
      return false;
    }
    OK = HNF_satisfy_constraints(B,U,hnf_col_order);
    if (!OK)
    {
      if (!expect_fail)
      {
        std::cout << "ERROR: HNF satisfy_constraints failed\n";
      }
      return false;
    }
    
    // vectors might be different length, this isn't the sort of expected failure that's allowed
    // if ( col_solution.size() == expected_solution.size() &&  col_solution == expected_solution )
    if (used_col >= (int) expected_solution.size())
    {
      result->error_message("HNF solution is length %d > expected %d\n", used_col, (int) expected_solution.size());
      return false;
    }
    bool expected=true;
    if (used_col>=0)
    {
      // if there is a nullspace vector, add it so the first-index values match
      if (nullspace_vector && col_solution[0]!=expected_solution[0])
      {
        const int k = (expected_solution[0]-col_solution[0])/((*nullspace_vector)[0]);
        for (int i = 0; i <= used_col; ++i)
        {
          col_solution[i] += k * (*nullspace_vector)[i];
        }
      }
      // test solution == expected
      for (int i = 0; i <= used_col; ++i)
      {
        if (col_solution[i]!=expected_solution[i])
          expected = false;
      }
    }
    if (expected)
      return true;
    
    int degrees_of_freedom = (int)B.col_rows.size() - (int)B.rows.size();
    result->warning_message("HNF solution satisfies the constraints, but is a different solution than expected. Number of rows and cols suggest %d degrees of freedom, but it may be rank or column deficient.\n",degrees_of_freedom);
    result->info_message("%s.\n",
                         (nullspace_vector!=nullptr ? "One nullspace vector was given and already taken into account" : "No nullspace vector was specified"));
    auto old_size = col_solution.size();
    col_solution.resize(expected_solution.size());
    result->info_message("Found solution");  print_vec(result,col_solution);
    result->info_message("Expected");        print_vec(result, expected_solution);
    col_solution.resize(old_size);
    return false;
  }
  
  bool IIATester::test_problem0()
  {
    IAResultImplementation result;
    IIATester IIA0(&result);
    std::cout << "IIATester:test_problem0 start: ";
    result.info_message("one degree of freedom, multiple solutions, 4x3.\n");
    //
    //  [ -1 -1  0  2 ]          -1
    //  |  3  1  1  0 | = A,  b=  6  from x=[1 2 1 1]
    //  [  4  2 -1 -1 ]           6
    //
    const int nr=3;
    const int nc=4;
    IIA0.number_of_rows=nr;
    IIA0.number_of_cols=nc;
    IIA0.freeze_problem_size();
    IIA0.used_row=nr-1;
    IIA0.used_col=nc-1;
    
    IIA0.rows[0].cols = {0, 1, 3};    IIA0.rows[0].vals = {-1,-1, 2};
    IIA0.rows[1].cols = {0, 1, 2};    IIA0.rows[1].vals = {3,1,1,};
    IIA0.rows[2].cols = {0, 1, 2, 3}; IIA0.rows[2].vals = {4,2,-1,-1};
    IIA0.MatrixSparseInt::fill_in_cols_from_rows();
    IIA0.rhs = {-1, 6, 6};
    
    std::vector<int> x0 = {1, 2, 1, 1};
    
    // skip test_rref_constraints for now
    
    bool OK0 = IIA0.test_hnf_solver(x0);
    if (!OK0)
    {
      result.error_message("This problem has many solutions but the solver didn't find any or didn't find the expected one!\n");
    }
    bool rref_OK1 = IIA0.test_rref_constraints();
    
    bool rref_OK2 = IIA0.test_rref_improve();
    
    std::cout << "IIATester:test_problem0 end." << std::endl;
    return rref_OK1 && rref_OK2 && OK0;
  }
  
  
  bool IIATester::test_problem1()
  {
    IAResultImplementation result;
    IIATester IIA0(&result);
    std::cout << "IIATester:test_problem1 start: ";
    result.info_message("one degree of freedom, multiple solutions, 4x4 but one redundant row.\n");
    //
    // redundant constraint, as with 4 sides of a brick being mapped, but with a different b
    //  [  1 -1       ]          1
    //  |  1    -1    | = A,  b= 2  from x=[4 3 2 1], not a mapped solution
    //  |     1    -1 |          2
    //  [        1 -1 ]          1
    //
    const int nr=4;
    const int nc=4;
    IIA0.number_of_rows=nr;
    IIA0.number_of_cols=nc;
    IIA0.freeze_problem_size();
    IIA0.used_row=nr-1;
    IIA0.used_col=nc-1;
    
    IIA0.rows[0].cols = {0, 1};  IIA0.rows[0].vals = {1,-1};
    IIA0.rows[1].cols = {0, 2};  IIA0.rows[1].vals = {1,-1};
    IIA0.rows[2].cols = {1, 3};  IIA0.rows[2].vals = {1,-1};
    IIA0.rows[3].cols = {2, 3};  IIA0.rows[3].vals = {1,-1};
    IIA0.MatrixSparseInt::fill_in_cols_from_rows();
    IIA0.rhs = {1, 2, 2, 1};
    
    // expected solution
    std::vector<int> x0 = {4, 3, 2, 1};
    // nullspace row
    std::vector<int> m0 = {1, 1, 1, 1};

    bool rref_OK1 = IIA0.test_rref_constraints();
    
    bool HNF_OK = IIA0.test_hnf_solver(x0,&m0);
    if (!HNF_OK)
    {
      std::cout << "ERROR: This problem has many solutions but the solver didn't find any or didn't find the expected one!\n";
    }
    
    bool rref_OK2 = IIA0.test_rref_improve();
    
    std::cout << "IIATester:test_problem1 end." << std::endl;
    return rref_OK1 && rref_OK2 && HNF_OK;
  }
  
  
  bool IIATester::test_problem2()
  {
    IAResultImplementation result;
    IIATester IIA0(&result);
    std::cout << "IIATester:test_problem2 start: ";
    result.info_message("fully constrained, one solution, 4x4\n");
    //
    // modified last row to make it non-redundant
    //  [  1 -1       ]          1
    //  |  1    -1    | = A,  b= 2  from x=[4 3 2 1], not a mapped solution
    //  |     1    -1 |          2
    //  [       1  -2 ]          0
    //
    const int nr=4;
    const int nc=4;
    IIA0.number_of_rows=nr;
    IIA0.number_of_cols=nc;
    IIA0.freeze_problem_size();
    IIA0.used_row=nr-1;
    IIA0.used_col=nc-1;
    
    IIA0.rows[0].cols = {0, 1};  IIA0.rows[0].vals = {1,-1};
    IIA0.rows[1].cols = {0, 2};  IIA0.rows[1].vals = {1,-1};
    IIA0.rows[2].cols = {1, 3};  IIA0.rows[2].vals = {1,-1};
    IIA0.rows[3].cols = {2, 3};  IIA0.rows[3].vals = {1,-2};
    IIA0.MatrixSparseInt::fill_in_cols_from_rows();
    IIA0.rhs = {1, 2, 2, 0};
    
    std::vector<int> x0 = {4, 3, 2, 1};
    
    bool rref_OK1 = IIA0.test_rref_constraints();
    
    bool HNF_OK = IIA0.test_hnf_solver(x0);
    if (!HNF_OK)
    {
      std::cout << "ERROR: This problem has a solution but the solver didn't find it!\n";
    }
    
    bool rref_OK2 = IIA0.test_rref_improve();
    
    std::cout << "IIATester:test_problem2 end." << std::endl;
    return rref_OK1 && rref_OK2 && HNF_OK;
  }
  
  bool IIATester::test_problem3()
  {
    IAResultImplementation result;
    IIATester IIA0(&result);
    std::cout << "IIATester:test_problem3 start: ";
    result.info_message("over-constrained, no solution, 5x4\n");
    //
    // redundant constraint, as with 4 sides of a brick being mapped, but with a different b
    //  [  1 -1       ]          1
    //  |  1    -1    | = A,  b= 2  from x=[4 3 2 1], not a mapped solution
    //  |     1    -1 |          2
    //  [        1 -2 ]          0
    //  [  1  1  1  1 ]          3
    //
    const int nr=5;
    const int nc=4;
    IIA0.number_of_rows=nr;
    IIA0.number_of_cols=nc;
    IIA0.freeze_problem_size();
    IIA0.used_row=nr-1;
    IIA0.used_col=nc-1;
    
    IIA0.rows[0].cols = {0, 1};  IIA0.rows[0].vals = {1,-1};
    IIA0.rows[1].cols = {0, 2};  IIA0.rows[1].vals = {1,-1};
    IIA0.rows[2].cols = {1, 3};  IIA0.rows[2].vals = {1,-1};
    IIA0.rows[3].cols = {2, 3};  IIA0.rows[3].vals = {1,-1};
    IIA0.rows[4].cols = {0, 1, 2, 3};  IIA0.rows[4].vals = {1,1,1,1};
    IIA0.MatrixSparseInt::fill_in_cols_from_rows();
    IIA0.rhs = {1, 2, 2, 0, 3};
    
    std::vector<int> x0 = {4, 3, 2, 1};
    
    bool rref_OK1 = IIA0.test_rref_constraints();
    
    result.log_warning = false;
    if (result.log_warning)
    {
      std::cout << "\nThe problem coefficients are designed to have no integer solution."
      " A warning like\n'WARNING: HNF no integer solution because row_sum[0] %% coeff[0] != 0:  -9 %% 4 = -1' is expected.\n";
    }
    bool HNF_OK = IIA0.test_hnf_solver(x0,nullptr,true);
    if (HNF_OK)
    {
      std::cout << "ERROR: This problem has no solution but the solver thought it found one!\n";
    }
    
    result.log_warning = true;
    bool rref_OK2 = IIA0.test_rref_improve();
    
    std::cout << "IIATester:test_problem3 end." << std::endl;
    return rref_OK1 && rref_OK2 && !HNF_OK;
  }
} // namespace IIA_Internal

int main(int argc, const char * argv[])
{
  
  // test Reduce Row Echelon RREF and Hermite Normal Form HNF matrix routines
  using IIA_Internal::IIATester;
  IIATester::test_problem0();
  IIATester::test_problem1();
  IIATester::test_problem2();
  IIATester::test_problem3();
  
  // test solving some interval assignment problem
  test_problem_A();
  test_problem_AA();
  test_problem_B();
  test_problem_BB();
  test_problem_BC();
  test_problem_C();
  test_problem_D();
  test_problem_E();
  
  test_problem_F();
  test_problem_G();

}
