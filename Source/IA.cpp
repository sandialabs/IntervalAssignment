// IA.cpp

#include "IA.h"
#include "IncrementalIntervalAssignment.h"
#include "IAResultImplementation.h"
#include <limits>
#include <iostream>

namespace IIA
{
  using IIA_Internal::IncrementalIntervalAssignment;
  using IIA_Internal::IAResultImplementation;

  // can't be inline because we need the non-stubbed version of IncrementalIntervalAssignment
  IA::IA() : free_row(0)
  {
    auto r = new IAResultImplementation();
    result = r;
    ia = new IncrementalIntervalAssignment(r);
  }
  
  IA::~IA()
  {
    delete result;
    result = nullptr;
    delete ia;
    ia = nullptr;
  }
  
  
  //== space for problem
  void IA::clear()
  {
    // could do this more efficiently...
    delete ia;
    ia = new IncrementalIntervalAssignment(&result);
    free_row = 0;
  }
  
  void IA::reserve(size_t nrows, size_t ncols)
  {
    ia->reserve_rows(nrows);
    ia->reserve_cols(ncols);
  }
  void IA::reserve_rows(size_t nrows)
  {
    ia->reserve_rows(nrows);
  }
  void IA::reserve_cols(size_t ncols)
  {
    ia->reserve_cols(ncols);
  }
  
  void IA::resize(size_t nrows, size_t ncols)
  {
    ia->reserve_rows(nrows);
    ia->reserve_cols(ncols);
    ia->freeze_problem_size();
  }
  void IA::resize_rows(size_t nrows)
  {
    ia->reserve_rows(nrows);
    ia->freeze_problem_size();
  }
  void IA::resize_cols(size_t ncols)
  {
    ia->reserve_cols(ncols);
    ia->freeze_problem_size();
  }
  
  std::pair<size_t, size_t> IA::size()
  {
    return std::make_pair(size_rows(), size_cols());
  }
  size_t IA::size_rows()
  {
    return ia->num_rows();
  }
  size_t IA::size_cols()
  {
    return ia->num_cols();
  }
  int IA::next_row()
  {
    return ia->next_available_row(EQ);
  }
  void IA::set_row(int row, std::vector<int> &cols, std::vector<int> &vals)
  {
    ia->set_M(row, cols, vals);
  }
  void IA::set_row_col(int row, int col, int val)
  {
    ia->set_M(row,col,val);
  }
  void IA::set_constraint(int row, ConstraintType constraint_type)
  {
    ia->set_constraint(row, constraint_type);
  }
  
  void IA::set_rhs(int row, int val)
  {
    ia->set_B(row,val);
  }
  
  // x
  void IA::set_bounds(int col, int lo, int hi)
  {
    ia->set_lower(col,lo);
    ia->set_upper(col,hi);
  }
  void IA::set_bound_lo(int col, int lo)
  {
    ia->set_lower(col,lo);
  }
  void IA::set_bound_hi(int col, int hi)
  {
    ia->set_upper(col,hi);
  }
  void IA::set_no_bound_lo(int col)
  {
    ia->set_lower(col,std::numeric_limits<int>::lowest());
  }
  void IA::set_no_bound_hi(int col)
  {
    ia->set_upper( col, std::numeric_limits<int>::max() );
  }
  void IA::set_no_bounds(int col)
  {
    ia->set_lower( col, std::numeric_limits<int>::lowest() );
    ia->set_upper( col, std::numeric_limits<int>::max() );
  }
  // goals must be in range (0,infinity) (int max in practice)
  // goals are scaled to be in this range and the modified value is returned
  //   if goals outside this range are required, you must scale the col variable instead
  double IA::set_goal(int col, double goal)
  {
    return ia->set_goal(col,goal);
  }
  void IA::set_no_goal(int col)
  {
    return ia->set_no_goal(col);
  }
  
  // get versions of the set methods
  void IA::get_row(int row, const std::vector<int> *&cols, const std::vector<int> *&vals) const
  {
    return ia->get_M(row,cols,vals);
  }
  int IA::get_row_col(int row, int col) const
  {
    return ia->get_M(row,col);
  }
  ConstraintType IA::get_constraint(int row) const
  {
    return ia->get_constraint(row);
  }
  int IA::get_rhs(int row) const
  {
    return ia->get_B(row);
  }
  
  void IA::get_bounds(int col, int &lo, int &hi) const
  {
    lo = ia->get_lower(col);
    hi = ia->get_upper(col);
  }
  int IA::get_bound_lo(int col) const
  {
    return ia->get_lower(col);
  }
  int IA::get_bound_hi(int col) const
  {
    return ia->get_upper(col);
  }
  bool IA::has_goal(int col) const
  {
    return ia->get_goal(col) != 0.;
  }
  double IA::get_goal(int col) const
  {
    return ia->get_goal(col);
  }
  
  // solve
  bool IA::solve()
  {
    return ia->solve(false, true, 3, true); // to do, 4th argument is not always true...
  }
  bool IA::solve_feasible()
  {
    assert(0); // not implemented yet, we don't have the flags for skipping improvement
    return ia->solve(false, true, 3, true); // to do, 4th argument is not always true...
  }
  bool IA::is_solved()
  {
    return ia->get_is_solved();
  }

  int IA::get_solution(int col) const
  {
    return ia->get_solution(col);
  }
  
  const std::vector<int> &IA::get_solution() const
  {
    return ia->get_solution();
  }
  
} // namespace

