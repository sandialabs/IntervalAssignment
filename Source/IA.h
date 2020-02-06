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

namespace IIA_Internal
{
  class IncrementalIntervalAssignment;
}

namespace IIA
{

class IAResult
{
public:
	IAResult();
	virtual ~IAResult();

	bool solved;
	bool constraints_satisfied;
	bool bounds_satisfied;
	bool optimized;
	bool error;
	int log_start; // first index of log entries for latest method call
	std::vector<std::string> log;
};

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
	virtual ~IA();

	//== space for problem
	void clear();

	void reserve(size_t nrows, size_t ncols);
	void reserve_rows(size_t nrows);
	void reserve_cols(size_t ncols);

	void resize(size_t nrows, size_t ncols);
	void resize_rows(size_t nrows);
	void resize_cols(size_t ncols);

	size_t size();
	size_t size_rows();
	size_t size_cols();

	//== Define Problem

	// =
	enum ConstraintType {EQ, LE, GE, EVEN, BAD};
	//                   ==, <=, >=, =2k,

	// A
	int next_row(); // optional, returns index of next unused row 0..	
	// swaps vectors with current contents, doesn't make a copy
	void set_row(int row, std::vector<int> &cols, std::vector<int> &vals);
	void set_row_col(int row, int col, int val); // A(row,col) = val
	// =
	void set_constraint(int row, ConstraintType constraint_type);
	// b
	void set_rhs(int row, int val); // b(row) = val

	// x
	// to set no bounds
	//   set_bound_hi( c, std::numeric_limits<int>::max() );
    //   set_bound_lo( c, std::numeric_limits<int>::lowest() );
	void set_bounds(int col, int lo, int hi); 
	void set_bound_lo(int col, int lo);
	void set_bound_hi(int col, int hi);
	// goals must be in range (0,infinity) (int max in practice)
	// goals are scaled to be in this range and the modified value is returned
	//   if goals outside this range are required, you must scale the col variable instead
	double set_goal(int col, double goal);
	double set_no_goal(int col);

	//== Defaults
	//   rows are sparse: vectors of column indices and coefficient values
    //   constraint default to equality 
    //   rhs b defaults to 0
    //   bounds default to [1,std::numeric_limits<int>::max()]
    //   goals default to 1.

	// get versions of the set methods
	void get_row(int row, const std::vector<int> &cols, const std::vector<int> &vals);
	int get_row_col(int row, int col);
	ConstraintType get_constraint(int row);
	int get_rhs(int row);
	void get_bounds(int col, int &lo, int &hi); 
	int get_bound_lo(int col);
	int get_bound_hi(int col);
	bool has_goal(int col);
	double get_goal(int col);

	//== Solve

	// sets x
	void solve();
	void solve_feasible(); // solves constraints and bounds only, ignores goals.

	// If you solve the problem, then add rows, the current solution is not thrown out.
	// Instead, (equality) rows are added with slack variables whose current value is out of bounds.
	// To get everything in bounds, call solve.

	int get_solution(int col);
	const std::vector<int> &get_solution();

	// can check this for "error" after each operation
	const IAResult &get_result();

private:
  IIA_Internal::IncrementalIntervalAssignment *ia;
	int free_row;

};

} // namespace
#endif
