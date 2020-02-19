// IA.h
// Public Interface to Interval Assignment (IA)
// Interval assignment is deciding the number of mesh edges (intervals) on model curves to be 
//   * compatible with the constraints imposed by quad and hex meshes, and the meshing algorithm 
//   * close to the user-desired sizes (goals)

#ifndef INTERVAL_ASSIGNMENT_H
#define INTERVAL_ASSIGNMENT_H

#include <vector>
#include <string>

#include "IAEnums.h"
#include "IAResult.h"

namespace IIA_Internal
{
  class IncrementalIntervalAssignment;
}

namespace IIA
{
  
  class IA
  {
  public:
    //== Problem
    // lexicographic min max f(x) :
    //    Ax=b
    //    x \in [lo,hi]
    // and f(x) is approximately x/g if x>g, else g/x
    // A is a matrix, x,b,lo,hi,g are vectors.
    
    IA();
    IA(IA &copy_me);
    void copy_me(IA &target);
    virtual ~IA();    
    
    //== space for problem
    void clear();
    
    void reserve(size_t nrows, size_t ncols);
    void reserve_rows(size_t nrows);
    void reserve_cols(size_t ncols);
    
    void resize(size_t nrows, size_t ncols);
    void resize_rows(size_t nrows);
    void resize_cols(size_t ncols);
    
    std::pair<size_t, size_t> size() const;
    size_t size_rows() const;
    size_t size_cols() const;
    
    //== Define Problem
    
    // A
    int next_row(); // optional, returns index of next unused row 0..
    int next_col(); // optional, returns index of next unused col 0..
    
    // swaps vectors with current contents, doesn't make a copy
    void set_row(int row, std::vector<int> &cols, std::vector<int> &vals);
    void set_row_col(int row, int col, int val); // A(row,col) = val
                                                 // =
    void set_constraint(int row, ConstraintType constraint_type);
    // b
    void set_rhs(int row, int val); // b(row) = val
    
    // x
    void set_bounds(int col, int lo, int hi);
    void set_bound_lo(int col, int lo);
    void set_bound_hi(int col, int hi);
    void set_no_bound_lo(int col);
    void set_no_bound_hi(int col);
    void set_no_bounds(int col);
    
    // goals must be in range (0,infinity) (int max in practice)
    // goals are scaled to be in this range and the modified value is returned
    //   if goals outside this range are required, you must scale the col variable instead
    double set_goal(int col, double goal);
    void set_no_goal(int col);
    
    //== Defaults
    //   rows are sparse: vectors of column indices and coefficient values
    //   constraint default to equality
    //   rhs b defaults to 0
    //   bounds default to [1,std::numeric_limits<int>::max()]
    //   goals default to 1.
    //   EVEN constraint rows are constrained to be an even number >= 4.
    //     To change the lower bound to a bigger/smaller value, set the rhs. Or
    //     Set up the sum-even variable xE yourself and add an EQ constraint x1 + x2 + ... - 2xE = b
    
    // get versions of the set methods
    void get_row(int row, const std::vector<int> *&cols, const std::vector<int> *&vals) const;
    int get_row_col(int row, int col) const;
    ConstraintType get_constraint(int row) const;
    int get_rhs(int row) const;
    void get_bounds(int col, int &lo, int &hi) const;
    int get_bound_lo(int col) const;
    int get_bound_hi(int col) const;
    bool has_goal(int col) const;
    double get_goal(int col) const;
    
    //== Solve
    
    // sets x
    bool solve_from_scratch(); // solves, overwriting any prior solution
    bool solve();              // solves, if there is a prior solution, we continue from there
    bool solve_feasible();     // solves constraints and bounds only, ignores goals.
    bool is_solved();
    
    // If you solve the problem, then add rows, the current solution is not thrown out.
    // Instead, (equality) rows are added with slack variables whose current value is out of bounds.
    // To get everything in bounds, call solve.
    
    // after solving, add some more rows (or cols). Returns the index of the first new one.
    int new_row(int num_rows=1);
    int new_col(int num_cols=1);
    // Remove a row by clearing it
    //   otherwise, it is *not* OK to modify a row after solving and then resolve; create a new row instead.
    void clear_row(int row);

    int get_solution(int col) const;
    const std::vector<int> &get_solution() const;
    
    // can set which types of messages are logged
    // can check this for "error" after each operation
    // see IAResult.h for details
    const IAResult *get_result() const;
    
  private:
    IIA_Internal::IncrementalIntervalAssignment *ia;
    int free_row;
    IAResult *result;
  };
  
#include "IAInline.h"
  
} // namespace
#endif
