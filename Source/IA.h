// IA.h
// Public Interface to Interval Assignment (IA)
// Interval assignment is deciding the number of mesh edges (intervals) on model curves to be 
//   * compatible with the constraints imposed by quad and hex meshes, and the meshing algorithm 
//   * close to the user-desired sizes (goals)
//
//  min lex f(x), where f(x) = (x>g ? x/g : g/x)
//  s.t. Ax = b
//       x in [lo,hi]
//       A, x, b integer
//
// Inputs are matrix A, and vectors b, g, lo, hi, and for each row whether "=" in Ax=b is equality or inequality or "sum-is-even".
// Outputs are vector x.
//
// Intended usage
//   Define the problem using set_row, set_bounds, set_goals, set_constraint
//   Solve
//   If needed (e.g. check whether the solution is good for your intended application)
//      add new constraints
//      modify old constraints, goals
//      solve again (uses prior solution as a warm-start if things were only added, not modified)
//
//   Space management
//     Caller can "resize" the problem to be big enough then select variable and row indices directly, and/or
//     call "next_row" and "next_col" to dynamically increase the problem size as needed (reserve may improve performance).
//
// Defaults
//   A starts as all zeros
//   constraints default to equality
//     EVEN constraint rows are constrained to be an even number >= 4.
//       To change the lower bound to a bigger/smaller value, set the rhs. Or
//       Set up the sum-even variable xE yourself and add an EQ constraint x1 + x2 + ... - 2xE = b
//   rhs b defaults to 0
//   [lo,hi] bounds default to [1,std::numeric_limits<int>::max()]
//   goals default to 1.

#ifndef INTERVAL_ASSIGNMENT_H
#define INTERVAL_ASSIGNMENT_H

#include <vector>

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
    
    // constructor
    IA();
    IA(IA &copy_me);
    void copy_me(IA &target);
    virtual ~IA();    
    
    //== space management
    //   analogous to std::vector methods
    void clear();
    
    void reserve     (size_t nrows, size_t ncols);
    void reserve_rows(size_t nrows);
    void reserve_cols(size_t ncols);
    
    void resize     (size_t nrows, size_t ncols);
    void resize_rows(size_t nrows);
    void resize_cols(size_t ncols);
        
    // A indices
    int next_row(); // index of next unused row 0..
    int next_col(); // index of next unused col 0..

    //== define problem
    
    // A
    //   rows are sparse: vectors of column indices and coefficient values
    void set_row    (int row, std::vector<int> &cols, std::vector<int> &vals); // swaps vectors, doesn't make a copy
    void set_row_col(int row, int col, int val); // A(row,col) = val
          
    // Ax "=" b.  What do you mean by "="?
    //            enum ConstraintType {EQ, LE, GE, EVEN};
    //                                  =  <=  >=   =2k 
    void set_constraint(int row, ConstraintType constraint_type);

    // b
    void set_rhs(int row, int val); // b(row) = val

    // x in [lo, hi]
    void set_bounds     (int col, int lo, int hi);
    void set_bound_lo   (int col, int lo);
    void set_bound_hi   (int col, int hi);
    void set_no_bounds  (int col);
    void set_no_bound_lo(int col);
    void set_no_bound_hi(int col);

    // g objective for x
    //   goals must be in range (0,infinity), where double infinity is limited to INT_MAX in practice.
    //   goals are scaled to be in this range and the modified value is returned.
    //     if goals outside this range are required, you must scale the variable instead.
    double set_goal   (int col, double goal);
    void   set_no_goal(int col); // "don't care" variable
        
    //== get methods
    
    //    size
    std::pair<size_t, size_t> size()      const;
    size_t                    size_rows() const;
    size_t                    size_cols() const;
    std::pair<int, int> used()      const;
    int                 used_rows() const;
    int                 used_cols() const;
    
    //    A
    void get_row    (int row, const std::vector<int> *&cols, const std::vector<int> *&vals) const;
    int  get_row_col(int row, int col) const;

    // get_col only works after "solve" or "fill_in_cols_from_rows" has been called. Use get_row_col instead.
    void get_col    (int col, const std::vector<int> *&rows) const;
    void fill_in_cols_from_rows();

    //   "="
    ConstraintType get_constraint(int row) const;
    
    //   b
    int get_rhs(int row) const;
    
    //   [lo,hi]
    void get_bounds(int col, int &lo, int &hi) const;
    int  get_bound_lo(int col) const;
    int  get_bound_hi(int col) const;
    
    //    g
    bool   has_goal(int col) const; // false for a "don't care" variable
    double get_goal(int col) const;
    
    //== Solve
    
    // solve for x
    bool solve();              // Recommended. If there is a prior solution, we continue from there.
    bool solve_from_scratch(); // Ignores any prior solution and starts from the beginning.
    void solve_from_scratch_next_time(); // next call to "solve" will not continue from a prior solution.
    bool solve_feasible();     // solves constraints and bounds only, ignores goals.

    // x
    bool is_solved() const;
    int                     get_solution(int col) const;
    const std::vector<int> &get_solution()        const;

    // Add more constraints
    //   after solving, add some more rows (or cols). Returns the index of the first new one.
    int new_row(int num_rows=1);
    int new_col(int num_cols=1);
    //   remove a row by clearing it
    void clear_row(int row);
    
    // i/o
    //   see IAResult.h
    //   set which types of messages are logged
    //   check for "error" after each operation
    const IAResult *get_result() const;
    
  private:
    // data
    IIA_Internal::IncrementalIntervalAssignment *ia;
    IAResult *result;
  };
  
#include "IAInline.h"
  
} // namespace
#endif
