//- Class: IncrementalIntervalAssignment
//- Owner: Scott Mitchell
//- Description: Implementation of setting up and solving interval assigment.
//- Checked By: 
//- Version: $Id:

#ifndef INTERVAL_ASSIGNMENT_IMPLEMENTATION_H
#define INTERVAL_ASSIGNMENT_IMPLEMENTATION_H

#include <cassert>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <limits>
#include <string>

using std::vector;
using std::set;
using std::map;
using std::priority_queue;
using std::numeric_limits;
using std::string;

#include "IAEnums.h"

namespace IIA
{
  class IAResult;
}

namespace IIA_Internal
{
  using IIA::IAResult;
  using IIA::ConstraintType;
  
// column, <negative blocking coeff, positive blocking coeff>
// many map entries are undefined
// if negative is -1 and positive is +1, then we can't change the variable at all
typedef std::map<int, std::pair<int,int> > BlockingCols;
const auto block_min = std::numeric_limits<int>::min();
const auto block_max = std::numeric_limits<int>::max();

struct RowSparseInt
{
  std::vector<int> cols, vals;
  
  RowSparseInt() {}
  RowSparseInt(std::vector<int> &c, std::vector<int> &v)   : cols(c), vals(v) {}
  RowSparseInt(std::vector<int> &&c, std::vector<int> &&v) : cols(c), vals(v) {}
  RowSparseInt(const RowSparseInt &copy) : cols(copy.cols), vals(copy.vals) {}

  void print_row() const;
  int get_val(int c) const; // return the value of the row in column c
  int get_col_index(int c) const; // return the index of c in cols
  void sort(); // sort by increasing column index. Also updates the vals.
};
struct MatrixSparseInt
{
  std::vector<RowSparseInt> rows;
  std::vector< std::vector<int> > col_rows;
  
  // nontrivial, adds up the number of non-zeros in all rows
  size_t num_nonzeros() const;

  // build the columns from the rows
  // Call once after rows are assigned and sorted
  void fill_in_cols_from_rows();

  void print_matrix_summary(std::string prefix=std::string());
  void print_matrix(std::string prefix=std::string());
  // s = u*s - v*r, with u and v chosen so that s[c]==0
  // standard call
  void matrix_row_us_minus_vr(int r, int s, int u, int v) { matrix_row_us_minus_vr(rows[r], rows[s], s, &col_rows, u, v); }
  // s = u*s - v*r, with u and v chosen so that s[c]==0
  void matrix_row_eliminate_column(int r, RowSparseInt &s_row, int c, int &u, int &v); // s_row not in matrix
  void matrix_row_eliminate_column(int r, RowSparseInt &s_row, int c) {int u,v; matrix_row_eliminate_column(r, s_row, c, u, v);} // s_row not in matrix
  void matrix_row_eliminate_column(int r, int s, int c, int &u, int &v);
  void matrix_row_eliminate_column(int r, int s, int c) {int u,v; matrix_row_eliminate_column(r,s,c,u,v);}
  void matrix_allrows_eliminate_column(int r, int c); // eliminate c from all rows except r
  
  // Gaussian elimination, specialized for creating a new_row that isn't blocked (that we know of).
  // starting from row 0, as if blocking_cols[0] was the first column etc.
  // return true if we are able to have row[unblocked_row] with a non-zero coefficient for col, and zeros for next_cols
  // and a col coefficient < max_coeff
  bool gaussian_elimination(int col, const BlockingCols &blocking_cols, int max_coeff, std::map<int,int> &improved_cols,
                            BlockingCols &all_blocking_cols, int &r, int &unblocked_row);
  // modify row (by row operations) so the coefficient of col < max_coeff
  //   here r is the last row that gaussian elimination diagonalized
  bool reduce_coefficient(int row, int col, int r, int max_coeff);
  // return true if row has a blocking coefficient when trying to increment col,
  // could modify to return the blocking cols
  bool is_row_blocked(int row, int col, const BlockingCols &blocking_cols );
  static
  bool is_blocked(int v, std::pair<int,int> block); // true if v is outside (int,int) as an integer interval
  // true if the coefficients for column col are such that their gcd >= max_coeff
  bool coefficient_irreducible(int col, int max_coeff);
protected:
  // unreduce c by row swaps and updating r; and remove c from all_blocking_cols
  void unreduce(int &r, int c, BlockingCols &all_blocking_cols);
  
  // counts the number of columns blocking bounds, i.e. a positive and negative block adds 2 to the count
  static
  size_t blocking_size(BlockingCols &blocking_cols);

  // workhorse implementation
  void matrix_row_us_minus_vr(const RowSparseInt &r_row, RowSparseInt &s_row, int s, std::vector< std::vector<int> > *col_rows_loc, int u, int v);

private:
  // re-used workspace for matrix_row_us_minus_vr
  std::vector<int> new_cols, new_vals;
protected:
  std::vector<int> kill_rows; 

public:

  // swap rows, update columns
  void matrix_row_swap(int r, int s);
  // swap cols, update the rows
  // void matrix_col_swap(int c, int d); // don't use, instead we keep a column order vector
  
  // T = transpose of this matrix if num_rows and num_cols are not given;
  //     else transpose of this's submatrix rows [0,num_rows), [0,num_cols)
  void transpose(MatrixSparseInt &T, size_t num_rows=0, size_t num_cols=0);
  void identity(MatrixSparseInt &I, size_t sz);
  
  // y = Ax, treating x and y as column vectors
  // return false if the problem is illformed
  bool multiply(const std::vector<int> &x, std::vector<int> &y);
  
  // set rsc to be the columns in row r or row s, in sorted order.
  void row_union(int r, int s, std::vector<int> &rsc);
  // set cdr to be the rows in column c or column d, in sorted order.
  void col_union(int c, int d, std::vector<int> &cdr);

  // create a copy of row as the new last row, return the row number
  int push_row(const RowSparseInt &row);
  int push_row(int r);
  // throw away the last row
  void pop_row();
  
  // push column col of this to matrix A
  // void push_column(MatrixSparseInt &A, int col);

  void copy(MatrixSparseInt&M,int MaxMRow);
  // copy constructor
  MatrixSparseInt(MatrixSparseInt&M,int MaxMRow);
  MatrixSparseInt(){}
};
struct EqualsB
{
  std::vector<int>rhs;
  std::vector<ConstraintType>constraint;
};

class IncrementalIntervalAssignment : 
protected MatrixSparseInt, protected EqualsB
{
protected: // data
  
  bool hasBeenSolved=false;
  int number_of_rows=0, used_row=-1, solved_used_row=-1;
  int number_of_cols=0, used_col=-1;
  
  std::vector<int> col_lower_bounds, col_upper_bounds, col_solution;
  std::vector<double> goals;
  // bounds are in [ std::numeric_limits<int>::lowest(), std::numeric_limits<int>::max() ]

  // names for debugging only. Usually these aren't assigned so incur minimal overhead
  bool names_exist();
  std::vector<std::string> row_names, col_names;
  std::string problem_name;

// from IntervalProblem
protected:

  IncrementalIntervalAssignment *parentProblem=nullptr;
  // the problem I'm a subproblem of

  std::vector<int> parentRows;
  // my row x is parentProblem's row parentRows[x]

  std::vector<int> parentCols;
  // std::vector<TDILPIntegerVariable*> intVars;
    int numSoftEdges=0;
  //- number of edges figuring into objective function
  
  int lastCopiedCol=0;
  //- If this is a subproblem, then lastCopiedCol = the number of
  //- copied columns from big ilp. Otherwise 0.
  
  int numDummies=0;
  //- number of dummy variables associated with this LP
  
  int sumEvenDummies=0;
  //- number of dummies that should be set to integer, usually for
  //- sum(I(e)) = even constraints.
  
  int usedDummies=0;
  
  int num_dummies();
  //- number of dummy variables in the IncrementalIntervalAssignment, for e.g. evenality constraints
  
  // the design is all interval matching implementations will have two passes
  //   map means solve everything except the sum-even constraints (possibly return a non-integer solution)
  //   even means solve all constraints, return an integer solution
  //   if row_min and row_max are passed in, then we only try to solve the subproblems involving those rows and ignore ones not containing those rows.
  int solve_map (int create_groups, int report_only, int create_infeasible_groups, int do_map, int do_pave, int row_min = -1, int row_max = -1 );
  int solve_even(int create_groups, int report_only, int create_infeasible_groups, int do_map, int do_pave, int row_min = -1, int row_max = -1 );

  // protected virtual methods
protected:
  // Add a new row in the middle of a solve, after the size has been frozen.
  //   Typical usage is for enforcing submap non-overlap "U" constraints, because there are two constraints we could add,
  //   and which one is best is not determined until we have a rough solution.
  // Return the index of the row
  // int new_row(MRow &Mrow);
  
  void print_associated_int_vars();

  int num_int_vars() { return intVars.size(); }
  //- number of TDILPIntegerVariables associated with this LP
  
  int add_dummy_variable(int num_variables, bool is_sum_even );
  //- add one or more dummy variables to the LP. Before freeze_size.
  
  int next_dummy();
  //- return next unused dummy variable. Use after freeze_size.
  //- Make sure all of the sum-even dummies are obtained first.



public:
  int parent_col(int c) {return parentCols[c];}


private:
    void recursively_add_edge( int int_var_column,
                             int do_sum_even,
                             std::vector <int> &sub_rows,
                             std::vector <int> &sub_cols,
                             std::vector<int> &sub_row_array,
                             std::vector<int> &sub_col_array );
  //- Add a column, and all non-sum-even rows for that column, and
  //- recursively add the other columns in those rows.
  //- Don't add a row or recurse if a row is a sum-even row and
  //- include_sum_even is false.
  //- sub_row_array[i] and sub_col_array[i] are (set to) TRUE if
  //- that row or column has been recusively added to the problem.
  //- sub_rows and sub_cols hold the index of the rows and columns,
  //- and are not sorted. Both structures are around to keep the
  //- running time proportional to subproblem size.

  //   if row_min and row_max are passed in, then we generate subproblems involving those rows and ignore ones not containing those rows.
  void subdivide_problem( std::vector<IncrementalIntervalAssignment*> &sub_problems,
                          bool do_sum_even = false, int row_min = -1, int row_max = -1 );
  //- Subdivide the equality (sum-even) constraints of this IncrementalIntervalAssignment into
  //- independent sub-problems, if do_sum_even is false (true).

  IncrementalIntervalAssignment* new_sub_problem();
  void delete_subproblems( std::vector<IncrementalIntervalAssignment*> &sub_problems );

  //- gather the solutions of the sub-problems into this problem
  void gather_solutions( std::vector<IncrementalIntervalAssignment*> &sub_problems );
  void gather_solution( IncrementalIntervalAssignment* sub_problem );

  int solve_sub(int create_infeasible_groups);
  int solve_sub_even();
  void copy_submatrix(std::vector <int> *rows, std::vector <int> *columns,
                              int *row_map, int *column_map, IncrementalIntervalAssignment *target );
  void copy_solution_to_sub  ( int i, IncrementalIntervalAssignment *sub_problem, int sub_i );
  void copy_solution_from_sub( int i, IncrementalIntervalAssignment *sub_problem, int sub_i );
  void copy_bounds_to_sub( IncrementalIntervalAssignment *sub_problem );


// original IIA
protected: // methods
  
  // layout of variable blocks:
  //   intVars, sumEvenDummies, convert inequality to equality dummies
  size_t intVar_start()         {return 0;}
  size_t intVar_end()           {return (size_t) intVars.size();}
  size_t allDummies_start()     {return (size_t) intVars.size();}
  size_t allDummies_end()       {return (size_t) intVars.size()+usedDummies;} // or numDummies
  // this
  // size_t nonevenDummies_start() {return (size_t) intVars.size()+sumEvenDummies;}
  // size_t nonevenDummies_end()   {return (size_t) intVars.size()+sumEvenDummies+usedDummies;} // or numDummies
  size_t sumEvenDummies_start() {return (size_t) intVars.size();}
  size_t sumEvenDummies_end()   {return (size_t) intVars.size()+sumEvenDummies;}
  //
  // not this
  // size_t nonevenDummies_start() {return (size_t) intVars.size();}
  // size_t nonevenDummies_end()   {return (size_t) intVars.size()+numDummies-sumEvenDummies;}
  // size_t sumEvenDummies_start() {return (size_t) intVars.size()+numDummies-sumEvenDummies;}
  // size_t sumEvenDummies_end()   {return (size_t) intVars.size()+numDummies;}
  
  // return the values of the non-zeros in column c
  //   we don't store these, so lookup takes a little bit of time
  std::vector<int> column_values(int c);
  
  // build the columns from the row data.
  // Call once after the problem is fixed and rows are sorted
  void fill_in_cols_from_rows();

  // sort each row's vector of cols and vals, by increasing cols (column index)
  void sort_rows();

  void print_rowIIA(size_t r);
  void print_colIIA(size_t c);

  
  // HNF Hermite Normal Form. Like RREF but using column operations. For solving initial constraints
  // This is A.  AU = (B 0) for column operations U.
  bool HNF(MatrixSparseInt &B, MatrixSparseInt &U, std::vector<int> &hnf_col_order);
  // use the HNF to solve Ax=b for x, i.e. satisfy the constraints but not the bounds
  bool HNF_satisfy_constraints(MatrixSparseInt &B, MatrixSparseInt &U, std::vector<int> &hnf_col_order);
  
  // transform data to reduced row echelon form
  //   return true if sucessful
  //   two ways of selecting pivots: pick coeff 2 sum-even dummies for bounds, but not constraints
  bool rref_constraints(std::vector<int> &rref_col_order); // for assigning dependent variables to satisfy constraints
  bool rref_improve    (std::vector<int> &rref_col_order); // for setting up nullspace to both satisfy bounds and improving quality
private:
  
  // utility for rref
  //   swap row rr into r and reduce column c
  //   i.e. pivot (rr,c) to (r) and add c to the col_order. remove c from all other rows. increment r
  bool rref_elim(int &r, int rr, int c, std::vector<int> &rref_col_order, std::vector<int> &rref_col_map);
  // various priorities for choosing next row to call rref_elim
  bool rref_step0(int &rref_r, std::vector<int> &rref_col_order, std::vector<int> &rref_col_map);
  bool rref_step1(int &rref_r, std::vector<int> &rref_col_order, std::vector<int> &rref_col_map);
  bool rref_step2(int &rref_r, std::vector<int> &rref_col_order, std::vector<int> &rref_col_map);
  bool rref_stepZ(int &rref_r, std::vector<int> &rref_col_order, std::vector<int> &rref_col_map);
  bool rref_step_numrows(int &rref_r, std::vector<int> &rref_col_order, std::vector<int> &rref_col_map);

protected:
  
  void row_swap(int r, int s);
  // void col_swap(int c, int d); // this messes up the layout of IntVars and dummies, instead keep a new order in a vector

  // multiply the row by some integer so that
  //   leading entry is positive and as small as possible
  // return the multiplier
  void row_simplify(int r);
  // get the greatest common divisor of matrix entries in this column, using only rows r and following
  int gcd_col(int r, int c);

  // replace row s with some multiple of row s and row r,
  //   so that in row s the column c entry is zero.
  void row_eliminate_column(int r, int s, int c);
  // row[s] = u*s - v*r
  void row_us_minus_vr(int r, int s, int u, int v);
private:
  // do the remaining stuff that the matrix versions don't
  void multiply_names_rhs(int r, int s, int u, int v);
protected:

  // The three main functions are
  //   satisfy_constraints:  initial feasible solution satisfying the constraints, Ax=b
  //   satisfy_bounds:       improve initial feasible solution so that all the variables are within their allowed ranges x \in [lo,hi]
  //   improve_solution:     improve the bounded-feasible solution so intervals are closer to the goal
  // satisfy_bounds and improve_solution work by finding a set of spanning vectors {x} for the nullspace of A, i.e. Ax=0, and using
  // those to move in a "downhill" direction
  
  //*** start: the following only work correctly after we've called rref

  // return true if calculated solution actually satisfies the constraints
  bool verify_constraints_satisfied(bool print=false);
  
  // find initial feasible solution, satisfying constraints
  // return true if successful
  bool satisfy_constaints(const std::vector<int>&cols_dep,
                          const std::vector<int>&cols_ind,
                          const std::vector<int>&rows_dep);
  
  // find initial feasible solution, satisfying constraints + variable bounds
  //   on input, solution already satisfies constraints, here we use nullspace M to increment intervals to satisfy bounds
  bool satisfy_bounds(MatrixSparseInt&M, int MaxMrow);

  // partition columns into dependent(leading) and independent(trailing) variables
  // returns true if some of the columns were swapped
  void categorize_vars(std::vector<int>&col_order,
                       std::vector<int>&cols_dep,
                       std::vector<int>&cols_ind,
                       std::vector<int>&rows_dep);
  
  // find all the rows containing the seed cols, then all the columns containing those rows
  void relevant_rows_cols(const std::vector<int>&seed_cols,
                          std::vector<int>&relevant_cols,
                          std::vector<int>&relevant_rows);
  
  // assign vars to their closest integer goal
  // only assigns intVars columns
  // if overwrite is false, then only variables who's value is 0 (uninitialized) are assigned.
  void assign_vars_goals(bool overwrite);
  // assign the dummy variables the closest value that would make the constraints feasible
  // returns true if the assignment is feasible
  bool assign_dummies_feasible();
  // adjust independent vars to try to make the dependent ones feasible
  bool adjust_independent_vars(std::vector<int>&cols_fail);
  // if the variables weren't assignable as integers, then return false and save the unsatisfied rows in rows_fail
  bool assign_dependent_vars(const std::vector<int>&cols_dep,
                             const std::vector<int>&rows_dep,
                             std::vector<int>&cols_fail);
  // return sum of coeff*solution - rhs
  static
  int row_sum(const std::vector<int>&cols, const std::vector<int>&coeff, const std::vector<int> &sol, int rhs);
  int row_sum(const std::vector<int>&cols, const std::vector<int>&coeff, int rhs);
  int row_sum(int r);
  bool bounds_satisfied(std::vector<int>&cols_fail);
  
  // set M's rows to the nullspace basis vectors of this
  bool create_nullspace(const std::vector<int>&cols_dep,
                        const std::vector<int>&cols_ind,
                        const std::vector<int>&rows_dep,
                        MatrixSparseInt&M,
                        int &MaxMrow);

  // make a feasible solution better, using rows of nullspace matrix M.
  // Searches for a local min of the max deviation from goals.
  void improve_solution(MatrixSparseInt&M, int MaxMrow, MatrixSparseInt&N);
  
  // make a feasible solution better, using rows of nullspace matrix M.
  // call after improve_solution is done,
  // see if there is anything else we can do to improve the lexicographic solution, even if the max cannot be improved
  void fine_tune(MatrixSparseInt&M, MatrixSparseInt&N);

  // create a copy of this problem into M
  // bool copy_matrix(MatrixSparseInt&M, int MaxMrow);

  //*** end: only works after we've called rref


protected: 

  void solution_init();
    
public:

  IncrementalIntervalAssignment( IAResult *result_ptr = nullptr );
  virtual ~IncrementalIntervalAssignment();  
  
  IAResult *result=nullptr;
  
  static void initialize_settings();
  // Tie global parameters to the UI

  void print_problem_summary(std::string prefix);
  void print_problem(std::string prefix);
  void print_problem(std::string prefix, const std::vector<int> &col_order);
  void print_solution(std::string prefix);
  void print_solution_row(int r, std::string prefix="");
  void print_solution_col(int c);
  void print_solution_summary(std::string prefix);
  void print();
  void dump_var(int c); // print debug info about the column and rows containing that column
  //- Debug. Print out the lp.

  void freeze_problem_size();
  void add_more_columns();
  void add_more_rows();


  int next_available_row( ConstraintType constraint_type);

  // int set_row_is_sum_even(int row, MRefEntity *ref_entity, double dummy_coeff, double dummy_bound);

  int solve( int create_groups, int report_flag, int scheme_flag, bool first_time );

  // return true if a row is something like "0 = 3" or some other obviously bad constraint
  bool infeasible_constraints();

  int get_is_solved() {return hasBeenSolved;}
  void set_is_unsolved() {hasBeenSolved=false;}


  //- returns the index = column of the first variable (if num_variables>1, rest are consecutive.
  int add_variable( int num_variables ); 
  //- returns the row of the first constraint. (if num_constraints>1, rest are consecutive.)
  int add_constraint( int num_constraints );

  int   get_M( int row, int col );
  void  set_M( int row, int col, int val ); // used?
    
  double set_goal(int c, double goal);
  void set_no_goal(int c);
  
protected:
  
  // for the tied-to columns we keep this
  struct tied_datum
  {
    int orig_bound_lo=0, orig_bound_hi=0;
    double goal_lo=0., goal_hi=0.;
  };
  std::map<int,tied_datum> tied_data;
  
  // variables that are constrained to be identical by the matrix
  //   index of columns that are constrained to be identical
  //   one variable keeps a vector of all other variables that are identical to it; its column is retained in the matrix
  //   the other variables get a vector of that master variable; their columns entries are zeroed out
  bool has_tied_variables = false;
  std::vector< std::vector<int> > tied_variables;
  bool find_tied_variables();
  void remove_tied_variables();
  // for tied variables, set the new col_solution to be the geometric average of the min and max tied variable
  bool should_adjust_solution_tied_variables = false;
  void adjust_solution_tied_variables();
  // if tied variables, then the nullspace will have spurious elementary rows for each removed variable
  void cull_nullspace_tied_variables(MatrixSparseInt&M, int &MaxMrow);
  // generate datum, returns false if the bounds are non-overlapping so infeasible
  bool generate_tied_data();
  // assign tied-variable to the same value as what they are tied-to.
  void assign_tied_variables();

  // amount we need to change the intervals to make it satisfy its variable bounds.
  // + = need to increase, - = need to decrease. Uses passed in value v rather than IIA's solution value
  int compute_out_of_bounds(int c, int v);
  // compute the number of goal intervals
  // int compute_goal(int c); // rounded to nearest int
  double compute_goal_float(int c);

  bool verify_full_solution(bool print_unsatisfied_constraints);

public:
  // checks if the variable is a tied-to variable, and uses the datum if so.
  // returns true if goal_lowest != goal_highest
  bool compute_tied_goals(int c, double &goal_lowest, double &goal_highest);

public: // public so MatrixSparseInt can use it
  // an element that we put into a priority queue, inherently tied to a variable column
  // to do, make a variant for constraint rows...
  class QElement
  {
  public:
    int c = -1; // column index
    
    bool solution_set = false;
    int solution = -1;
    
    // selection criteria
    double valueA = 0; // primary criteria
    double valueB = 0; // secondary
    double valueC = 0; // tertiary
    // int valueD = 0.; // 4th
    // int valueE = 0.; // 5th
    // c is last-resort arbitrary tiebreaker
                        // biggest elements are the ones at the top of the queue, end of the set
    
    int dx = 0; // desire increment +1, or decrement -1

    bool operator<(QElement const& rhs) const
    {
      if (valueA < rhs.valueA)
        return true;
      if (valueA > rhs.valueA)
        return false;
      if (valueB < rhs.valueB)
        return true;
      if (valueB > rhs.valueB)
        return false;
      if (valueC < rhs.valueC)
        return true;
      // we do not want a tiebreaker when assessing quality, but we do need one for sets
      // if (valueC > rhs.valueC)
      //  return false;
      // if (c < rhs.c)
      //   return true;
      // if (c > rhs.c), or if equal
      return false;
//      return
//      ( valueA <  rhs.valueA ) ||
//      ( valueA == rhs.valueA   && valueB <  rhs.valueB ) ||
//      ( valueA == rhs.valueA   && valueB == rhs.valueB   && valueC <  rhs.valueC ); // ||
      // ( value == rhs.value   && valueB == rhs.valueB   && valueC == rhs.valueC   && valueD <  rhs.valueD) ||
      // ( value == rhs.value   && valueB == rhs.valueB   && valueC == rhs.valueC   && valueD == rhs.valueD && valueE < rhs.valueE);
    }

    void print();
  };
  
  // ensure elements are unique if they refer to a different column
  struct QElement_less_for_sets
  {
    bool operator() (QElement const& lhs, QElement const& rhs) const
    {
      if (lhs.valueA < rhs.valueA)
        return true;
      if (lhs.valueA > rhs.valueA)
        return false;
      if (lhs.valueB < rhs.valueB)
        return true;
      if (lhs.valueB > rhs.valueB)
        return false;
      if (lhs.valueC < rhs.valueC)
        return true;
      if (lhs.valueC > rhs.valueC)
        return false;
      // we do not want a tiebreaker when assessing quality, but we do need one for sets
      if (lhs.c < rhs.c)
        return true;
      // if (c > rhs.c), or if equal
      return false;
    }
  };
  
protected:
  
  // set_values function. Updates value, valueB, valueC if intervals have changed.
  // return true if the new values are less than (lower priority than) the old values
  class SetValuesFn
  {
  public:
    // sets values based on current IIA.col_solution
    //   updates solution and returns true if solution has changed or was unset
    bool update_values(IncrementalIntervalAssignment &IIA, QElement &qe);
    // sets values based on current IIA.col_solution
    void set_values(IncrementalIntervalAssignment &IIA, QElement &qe);
    // uses the passed in solution instead of IIA's solution
    void set_values(IncrementalIntervalAssignment &IIA, QElement &qe, int solution);
  protected:
    // assumes qe.solution is already set, uses that instead of col_solution
    virtual void set_values_implementation(IncrementalIntervalAssignment &IIA, QElement &qe ) = 0;
    // virtual void print( const QElement &qe );
  };
  // specializations set the values based on differing criteria
  
protected:
  
  // ratio  (current-increment)/goal                  for current>goal
  //                       goal/(current+increment)   for goal>current
  class SetValuesRatioR: public SetValuesFn
  {
  public:
    SetValuesRatioR() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &IIA, QElement &qe ) override;
  private:
    void set_values_goal(IncrementalIntervalAssignment &IIA, double g, QElement &qe );
  };

  // ratio  current/goal     for current>goal
  //           goal/current  for goal>current
  class SetValuesRatio: public SetValuesFn
  {
  public:
    SetValuesRatio() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &IIA, QElement &qe ) override;
    void set_values_goal(IncrementalIntervalAssignment &IIA, double g, QElement &qe );
  };
  
  class SetValuesBounds: public SetValuesFn
  {
  public:
    SetValuesBounds() {}
  protected:
    void set_values_implementation(IncrementalIntervalAssignment &IIA, QElement &qe ) override;
  };
  // also want a queue element for a constraint?
  
  // for picking a column to eliminate in rref that has a single 1
  class SetValuesOneCoeff: public IncrementalIntervalAssignment::SetValuesFn
  {
  public:
    SetValuesOneCoeff() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &IIA, QElement &qe ) override;
  };

  // in order of increasing number of coefficients = number of rows the variable appears in
  class SetValuesNumCoeff: public IncrementalIntervalAssignment::SetValuesFn
  {
  public:
    SetValuesNumCoeff() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &IIA, QElement &qe ) override;
  };

  
  // pick column to eliminate in rref, that has a small coeff, is in few rows, and is long
  class SetValuesCoeffRowsGoal: public IncrementalIntervalAssignment::SetValuesFn
  {
  public:
    SetValuesCoeffRowsGoal() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &IIA, QElement &qe ) override;
  };



  // get the top of the Q, after ensuring it is up to date
  // if true, Q is empty and the returned QElement is just the default, and shouldn't be processed
   bool tip_top( std::priority_queue<QElement> &Q, SetValuesFn &val_fn, double threshold, QElement &t );

  // build a priority_queue, containing qcol and/or qrow, with priority set by set_val_fn
  // queue elements with a value not less than threshold are not put on the queue
  void build_Q(std::priority_queue< QElement > &Q,
               SetValuesFn &set_val_fn,
               double threshold,
               const std::vector<int> &qcol);
  
  // accept_strict_improvement
  //   if the solution strictly improves
  //      solution+=dx*mrow
  //      update Q
  //      return true
  //   else
  //      blocking_cols = vars that would get worse (accept self), return false
  
  class QWithReplacement
  {
  public:
    
    QWithReplacement(IncrementalIntervalAssignment *IIA_in, SetValuesFn &f, double f_threshold)
    : val_fn(&f), threshold(f_threshold), IIA(IIA_in) {}

    void build(const std::vector<int> &cols);

    // returns true if done
    bool tip_top( QElement &t );

    // update the Q
    // cols are the indices of the columns that have changed
    // col_solution is the new solution of the IncrementalIntervalAssignment
    void update(const std::vector<int> &cols);
    
    bool empty() {return Q.empty();}
    
    // debug
    void print();
    
    size_t size() {return Q.size();}
    
  private:
    std::map<int,QElement> elements;
    std::set<QElement,QElement_less_for_sets> Q;
    SetValuesFn *val_fn = nullptr;
    double threshold = 0.0;
    IncrementalIntervalAssignment *IIA = nullptr;
    
    // conditionally add a q element for column c
    void add(int c);
  };

  // used by both satisy_bounds and improve_solution
  //   bounds: define "improvement" by the priority decreasing
  //   improve: define "improvement" by the quality function, and satisfying the constraint function
  bool accept_strict_improvement(QWithReplacement &Q,
                                 SetValuesFn &quality_fn,
                                 double quality_threshold,
                                 SetValuesFn *constraint_fn, // satisfying_bounds if this is null, "improving" otherwise
                                 double constraint_threshold,
                                 int tc,
                                 RowSparseInt &R,
                                 int dx,
                                 bool &self_block,
                                 int &self_block_coeff,
                                 BlockingCols &blocking_cols,
                                 int &num_block,
                                 std::map<int,int> &improved_cols);
  
  bool is_improvement(QElement &s,
                      QElement &t,
                      SetValuesFn *constraint_fn,
                      double constraint_threshold );
  
  void compute_quality_ratio(std::vector<QElement> &q, std::vector<int> &cols);
  // increase the solution by dx times R. Returns true if the solution is (still) within bounds
  bool increment_solution(RowSparseInt &R, int dx);
  // true if A<B by lexicographic min max
  bool is_better( std::vector<QElement> &qA, std::vector<QElement> &qB);

  
  // io
public:
  enum MessageType {DEBUG_MSG, INFO_MSG, WARNING_MSG, ERROR_MSG };
  void log_message(MessageType message_type, const char* format, ...);
protected:
  void print_vec(const std::vector<int> &vec, bool lf = true);
  void print_vec(const std::vector< std::pair<int,int> > &vec);
  void print_set(const std::set< std::pair<int,int> > &aset);
  void print_map(const BlockingCols &amap);
  void print_map(const std::map<int,int> &amap);
};

inline
IncrementalIntervalAssignment::IncrementalIntervalAssignment( IAResult *result_ptr ) : result(result_ptr)
{
}

inline
IncrementalIntervalAssignment* IncrementalIntervalAssignment::new_sub_problem()
{
  auto np = new IncrementalIntervalAssignment(result);
  np->should_adjust_solution_tied_variables = should_adjust_solution_tied_variables;
  return np;
}

}

#endif  // INCREMENTALINTERVALASSIGNMENT_HPP
