//- Class: IncrementalIntervalAssignment
//- Description: Implementation of setting up and solving interval assigment.

#ifndef INTERVAL_ASSIGNMENT_IMPLEMENTATION_H
#define INTERVAL_ASSIGNMENT_IMPLEMENTATION_H

#include <vector>
#include <set>
#include <map>
#include <queue>
#include <limits>
#include <string>

#include "IAEnums.h"

namespace IIA_Internal
{
  using std::pair;
  using std::vector;
  using std::set;
  using std::map;
  using std::priority_queue;
  using std::numeric_limits;
  using std::string;
  
  using IIA::ConstraintType;

  // forward
  class IAResultImplementation;
  class QWithReplacement;
  class SetValuesFn;
  class QElement;
  
  // column, <negative blocking coeff, positive blocking coeff>
  // many map entries are undefined
  // if negative is -1 and positive is +1, then we can't change the variable at all
  typedef map<int, pair<int,int> > BlockingCols;
  const auto block_min = numeric_limits<int>::min();
  const auto block_max = numeric_limits<int>::max();
  
  struct RowSparseInt
  {
    vector<int> cols, vals;
    
    RowSparseInt() {}
    RowSparseInt(vector<int> &c, vector<int> &v)   : cols(c), vals(v) {}
    RowSparseInt(vector<int> &&c, vector<int> &&v) : cols(c), vals(v) {}
    RowSparseInt(const RowSparseInt &copy) : cols(copy.cols), vals(copy.vals) {}
    
    void multiply(int k);     // row = k*row
    int dot(const RowSparseInt &b); // return this * b as a dot product

    void print_row(IAResultImplementation *result) const;
    int get_val(int c) const; // return the value of the row in column c
    int get_col_index(int c) const; // return the index of c in cols
    void sort(); // sort by increasing column index. Also updates the vals.
    
    void sparse_to_full(vector<int> &full) const;
  };
  struct MatrixSparseInt
  {
  public: // methods

    // constructor
    MatrixSparseInt(IAResultImplementation *result_ptr) : result(result_ptr){}
    // copy constructor, copies up to row MaxMrow
    MatrixSparseInt(const MatrixSparseInt&M, int MaxMrow);
    // copy constructor, concatenate M[1:MaxMrow,:] and N
    MatrixSparseInt(const MatrixSparseInt&M, int MaxMrow, const MatrixSparseInt&N, bool no_duplicate_rows);

    const IAResultImplementation *get_result() const {return result;}

    // nontrivial, adds up the number of non-zeros in all rows
    size_t num_nonzeros() const;
    
    // build the columns from the rows
    // Call once after rows are assigned and sorted
    void fill_in_cols_from_rows();
    
    void print_matrix_summary(string prefix=string()) const;
    void print_matrix(string prefix=string()) const;
    
    // r += m*s
    // conceptually static: "this" unchanged, only used for "workspace"
    void row_r_add_ms(RowSparseInt &r, int m, const RowSparseInt &s )
    {
      matrix_row_us_minus_vr( s, r, -1, nullptr, 1, -m);
    }

    // s = u*s - v*r
    void matrix_row_us_minus_vr(int r, int s, int u, int v) { matrix_row_us_minus_vr(rows[r], rows[s], s, &col_rows, u, v); }
    // s = u*s - v*r, with u and v chosen so that s[c]==0, for all rows s containing c except r
    void matrix_allrows_eliminate_column(int r, int c, vector<int> *rhs = nullptr); // eliminate c from all rows except r
    virtual void allrows_eliminate_column(int r, int c); // with rhs=nullptr

    
    // Gaussian elimination, specialized for creating a new_row that isn't blocked (that we know of).
    // starting from row 0, as if blocking_cols[0] was the first column etc.
    // return true if we are able to have row[unblocked_row] with a non-zero coefficient for col, and zeros for next_cols
    // and a col coefficient < max_coeff
    bool gaussian_elimination(int col, const BlockingCols &blocking_cols, int max_coeff, map<int,int> &improved_cols,
                              BlockingCols &all_blocking_cols, int &r, int &unblocked_row);

    
    // swap rows, update columns
    virtual void row_swap(int r, int s);
    virtual void row_negate(int r); // row = -row
    // swap cols, update the rows
    // void matrix_col_swap(int c, int d); // don't use, instead we keep a column order vector
    
    // T = transpose of this matrix if num_rows and num_cols are not given;
    //     else transpose of this's submatrix rows [0,num_rows), [0,num_cols)
    void transpose(MatrixSparseInt &T, size_t num_rows=0, size_t num_cols=0) const;
    static void identity(MatrixSparseInt &I, size_t sz);
    
    // y = Ax, treating x and y as column vectors
    // return false if the problem is illformed
    bool multiply(const vector<int> &x, vector<int> &y) const;

    // set rsc to be the columns in row r or row s, in sorted order.
    void row_union(int r, int s, vector<int> &rsc);
    
    // create a copy of row as the new last row, return the row number
    int push_row(const RowSparseInt &row);
    
    // construct HermiteNormalForm for this
    // specifically, for the submatrix this[0..this_rows-1][0..this_cols-1]
    // assumes there are no zero rows
    bool HNF(int this_rows, int this_cols, MatrixSparseInt &B, MatrixSparseInt &U, vector<int> &hnf_col_order, vector<int> &g, bool use_HNF_goals) const;
    // find the submatrix this[0..r, 0..c] that prunes off non-zero rows and columns from this
    void nonzero(int &r, int &c) const;
    // true if Ax=b, where this==A
    bool verify_Ax_b(const vector<int> &x, const vector<int> &b) const;
    // least common multiple of non-zero
    int lcm_entries() const;
    
    // multiply the row by some integer so that
    //   leading entry is positive and as small as possible
    virtual void row_simplify(int r);
    
  public: // data
    IAResultImplementation *result;
    
    // public so IncrementalIntervalAssignment methods can access and manipulate them directly for nullspace matrices, etc.
    vector<RowSparseInt> rows;
    vector< vector<int> > col_rows;

  protected:
    
    void print_map(const BlockingCols &amap) const;

    // modify row (by row operations) so the coefficient of col < max_coeff
    //   here r is the last row that gaussian elimination diagonalized
    bool reduce_coefficient(int row, int col, int r, int max_coeff);

    // return true if row has a blocking coefficient when trying to increment col,
    // could modify to return the blocking cols
    bool is_row_blocked(int row, int col, const BlockingCols &blocking_cols );
    static
    bool is_blocked(int v, pair<int,int> block); // true if v is outside (int,int) as an integer interval
    
    // true if the coefficients for column col are such that their gcd >= max_coeff
    bool coefficient_irreducible(int col, int max_coeff);
    
    // unreduce c by row swaps and updating r; and remove c from all_blocking_cols
    void unreduce(int &r, int c, BlockingCols &all_blocking_cols);
    
    // counts the number of columns blocking bounds, i.e. a positive and negative block adds 2 to the count
    static
    size_t blocking_size(BlockingCols &blocking_cols);
    
    // workhorse implementation
    void matrix_row_us_minus_vr(const RowSparseInt &r_row, RowSparseInt &s_row, int s, vector< vector<int> > *col_rows_loc, int u, int v);
    
  private:
    // re-used workspace for matrix_row_us_minus_vr
    vector<int> new_cols, new_vals;
    
  };
  struct EqualsB
  {
    vector<int>rhs;
    vector<ConstraintType>constraint;
  };
  
  class IncrementalIntervalAssignment : 
  protected MatrixSparseInt, protected EqualsB
  {
    
  public: // interface
    
    IncrementalIntervalAssignment( IAResultImplementation *result_ptr );
    virtual ~IncrementalIntervalAssignment();
    
    const IAResultImplementation *get_result() const {return result;}

    void copy_me( IncrementalIntervalAssignment *target );

    // solve interval assignment Mx = B for x
    //   also works to re-solve a problem after adding new rows, continuing from the prior solution
    //   if first_time is set to true, then we solve from scratch and not re-solve no matter what.
    bool solve( bool do_improve = true, bool first_time = false );

    // Compare quality ratio of solution X to solution Y, lex min max, using is_better
    //   1 if qX<qY
    //  -1 if qY<qX
    //   0 if qX==qY
    // X and Y are conceptually const
    int solution_X_is_better_than_Y( vector<int> &X, vector<int> &Y, bool print_summary, bool print_detail );
    // X == col_solution
    int solution_is_better_than_Y( /*col_solution,*/ vector<int> &Y, bool print_summary, bool print_detail );

    // has the problem been solved. Changing the problem will set this flag to false
    int get_is_solved() const;
    // set_is_unsolved will also ensure than next call to solve will start from scratch
    void set_is_unsolved();

    // problem size
    int num_rows() const {return number_of_rows;}
    int num_cols() const {return number_of_cols;}
    int num_used_rows() const {return used_row+1;}
    int num_used_cols() const {return used_col+1;}
    // keeps track of how many rows and columns we will allocate when we freeze size
    void reserve_cols( int num_cols );
    void reserve_rows( int num_rows );
    // actually allocate storage for the size we reserved
    void freeze_problem_size();
    void resize_cols( int num_cols );
    void resize_rows( int num_rows );
    
    // return next unused row, so we can fill it in
    int next_available_row( ConstraintType constraint_type = IIA::EQ);
    
    // return next unused variable, as a dummy variable or an interval variable
    int next_dummy();
    int next_intvar();

    // constraint matrix values Mx = B
    int    get_M_unsorted( int row, int col ) const;
    int    get_M( int row, int col ) const;  // caution: relies on entries being sorted, be careful when we allow queries
    void   get_M( int row, const vector<int> *&cols, const vector<int> *&vals ) const;
    void   set_M( int row, int col, int val ); // doesn't rely on entries being sorted, may be slow
    void   set_M( int row,  vector<int> &cols, vector<int> &vals ); // modifes rows and cols
    void clear_M(int row); // clear the entries of this row. col_rows is invalid

    // only works after solve or fill_in_cols_from_rows has been called
    void get_rows(int col, const vector<int> *&rows) const { rows = &col_rows[col]; };

    // build the columns from the row data.
    //   It is called once automatically near the start of solve, or can be called explicitly if the caller needs get_rows often
    void fill_in_cols_from_rows();

    // type of constraint for each row, =, <=, >=, EVEN,...
    void set_constraint(int row, ConstraintType constraint_type);
    ConstraintType get_constraint(int row) const {return constraint[row];}

    // right hand side
    int  get_B( int row ) const;
    void set_B( int row, int val );
    
    // ideal value for each variable
    double  get_goal(int c) const; // no goal is internally represented by 0.
    double  set_goal(int c, double goal);
    void set_no_goal(int c); // dummy variables have no goal

    // variable bounds
    int  get_lower ( int col ) const {return col_lower_bounds[col];}
    int  get_upper ( int col ) const {return col_upper_bounds[col];}
    void set_lower ( int col, int bound );
    void set_upper ( int col, int bound );

    // to set no bounds
    //   set_lower( c, numeric_limits<int>::lowest() );
    //   set_upper( c, numeric_limits<int>::max() );

    // x
    int get_solution( int col )       const { return col_solution[col]; }
    const vector<int> &get_solution() const { return col_solution; }
    
    // Add a new row after the problem has already been solved.
    //   Typical usage is for enforcing submap non-overlap "U" constraints, because there are two constraints we could add,
    //   and which one is best is not determined until we have a rough solution.
    // Return the index of the row/col
    int new_row(int num_rows=1);
    int new_col(int num_cols=1);

    // print and debug
    void print_problem_summary(string prefix) const;
    void print_problem(string prefix) const;
    void print_problem(string prefix, const vector<int> &col_order) const;
    void print_solution(string prefix) const;
    void print_solution_row(int r, string prefix="") const;
    void print_solution_col(int c) const;
    void print_solution_summary(string prefix) const;
    void print() const;
    
    // research options
    bool satisfy_constraints_only = false; // Default false. If true, don't attempt to satisfy_bounds or improve_solution. For debugging.
    const bool turn_on_pave_research_code = true; // Default true. If true, then copy mapping nullspace to create local/intuitive paving nullspace vectors.

    // turn on or off some algorithm options to test output quality
    bool use_HNF_always = false; // false then try rref solution first
    bool use_HNF_goals = true;  // false c=1
    bool use_map_nullspace = true; // false use M only during sum-even phase
    bool use_best_improvement = true; // false use first improvement vector

    bool use_better_rref_step1_pivots_in_sumeven_phase = true; // Default true. If false, use the same pivot rules as in the mapping phase
    bool use_fixed_RatioR_straddle = true; // Default true. If false, reintroduce a bug for calculating RatioR for solutions less than 1 away from goals.
    bool use_smart_adjust_for_fractions = true; // Default true. If false, then we always round upwards
    bool use_Ratio_for_improvement_queue = true; // Default true. If false, then use RatioR for selecting which variable
    bool use_Nullspace_sumeven_to_flip_sign_of_blocker = true; // Default true. If true, then when trying to increase x0, if "x0 - x1" is blocked by x1, then add a "2x1 + y" Nullspace row to flip the sign of x1
    bool debug_Nullspace_sumeven_to_flip_sign_of_blocker = false; // output debugging print statements

  protected: // data
    
    bool hasBeenSolved=false;
    int number_of_rows=0, used_row=-1, solved_used_row=-1;
    int number_of_cols=0, used_col=-1, solved_used_col=-1;
    
    vector<int> col_lower_bounds, col_upper_bounds, col_solution;
    vector<double> goals;
    // bounds are in [ numeric_limits<int>::lowest(), numeric_limits<int>::max() ]
    // goals are in (0,numeric_limits<int>::max())
    
    // INT_VAR = interval variables, representing some number of intervals
    //   other vars are sum-even or slack dummy variables
    enum VarType {INT_VAR, EVEN_VAR, DUMMY_VAR, UNKNOWN_VAR};
    vector<VarType> col_type;
    
    // == used by sub_problems only
    
    IncrementalIntervalAssignment *parentProblem=nullptr;
    // the problem I'm a subproblem of
    
    vector<int> parentRows;
    // my row x is parentProblem's row parentRows[x]
    
    vector<int> parentCols;
    // my col y is parentProblem's column parentCols[y]
    
    int lastCopiedCol=0;
    //- If this is a subproblem, then lastCopiedCol = the number of
    //- copied columns from big ilp. Otherwise 0.
    
    MatrixSparseInt Nullspace;
    
  protected: // methods
    
    // the design is all interval matching implementations will have two passes
    //   map means solve everything except the sum-even constraints (possibly return a non-integer solution)
    //   even means solve all constraints, return an integer solution
    //   if row_min and row_max are passed in, then we only try to solve the subproblems involving those rows and ignore ones not containing those rows.
    bool solve_phase( bool map_only_phase, bool do_improve, int row_min, int row_max);

    // convert parent's mapping-phase nullspace vectors to sub_problem's sumeven-phase nullspace vectors, possibly containing dummy variables
    void transfer_parent_nullspace_to_subs(vector<IncrementalIntervalAssignment*> &sub_problems);

    bool tiny_subspace(int c, MatrixSparseInt &M, MatrixSparseInt &MN);
    // find a tiny subset of the matrix that contains the column c, and some some columns of the rows with c,
    //   and enough rows that the nullspace has dimension >= 1
    // add those new nullspace rows to M and MN
    // return false if we couldn't find a reasonable subspace

    void recursively_add_edge( int int_var_column,
                              int do_sum_even,
                              vector <int> &sub_rows,
                              vector <int> &sub_cols,
                              vector<int> &sub_row_array,
                              vector<int> &sub_col_array );
    //- Add a column, and all non-sum-even rows for that column, and
    //- recursively add the other columns in those rows.
    //- Don't add a row or recurse if a row is a sum-even row and
    //- include_sum_even is false.
    //- sub_row_array[i] and sub_col_array[i] are (set to) TRUE if
    //- that row or column has been recusively added to the problem.
    //- sub_rows and sub_cols hold the index of the rows and columns,
    //- and are not sorted. Both structures are around to keep the
    //- running time proportional to subproblem size.
    
    //   if row_min and row_max are >-1, then we generate subproblems involving those rows and ignore ones not containing those rows.
    void subdivide_problem(vector<IncrementalIntervalAssignment*> &sub_problems,
                           vector<int> &orphan_cols,
                           bool do_sum_even, int row_min, int row_max );
    //- Subdivide the equality (sum-even) constraints of this IncrementalIntervalAssignment into
    //- independent sub-problems, if do_sum_even is false (true).
    
    IncrementalIntervalAssignment* new_sub_problem(const vector <int> &sub_rows, const vector <int> &sub_cols);
    void delete_subproblems( vector<IncrementalIntervalAssignment*> &sub_problems );
    
    //- gather the solutions of the sub-problems into this problem
    void gather_solutions( vector<IncrementalIntervalAssignment*> &sub_problems, bool want_sub_nullspaces, vector<int> &orphan_cols );
    
    int solve_sub(bool do_improve, bool map_only_phase);
    void copy_submatrix(const vector <int> &rows, const vector <int> &cols,
                        IncrementalIntervalAssignment *target ) const;
    void copy_bounds_to_sub( IncrementalIntervalAssignment *sub_problem ) const;
    
    // force equality row r to be satisfied by introducing a slack variable r
    void force_satisfied_with_slack(int r);
    
    // return the non-zeros values of the matrix in column c
    //   we don't store these, so lookup takes a little bit of time
    vector<int> column_coeffs(int c) const;
        
    // sort each row's vector of cols and vals, by increasing cols (column index)
    //   skips rows <row_start
    void sort_rows(int row_start);
    
    void print_row_iia(size_t r) const;
    void print_col_iia(size_t c) const;
    
    // transform data to reduced row echelon form
    //   return true if sucessful
    //   two ways of selecting pivots: pick coeff 2 sum-even dummies for bounds, but not constraints
    bool rref_constraints(vector<int> &rref_col_order, bool map_only_phase); // for assigning dependent variables to satisfy constraints
    bool rref_improve    (vector<int> &rref_col_order); // for setting up nullspace to both satisfy bounds and improving quality
    
    // utility for rref
    //   swap row rr into r and reduce column c
    //   i.e. pivot (rr,c) to (r) and add c to the col_order. remove c from all other rows. increment r
    bool rref_elim(int &r, int rr, int c, vector<int> &rref_col_order, vector<int> &rref_col_map);
    // various priorities for choosing next row to call rref_elim
    bool rref_step0(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map);
    bool rref_step1(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map, bool map_only_phase);
    bool rref_step2(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map);
    bool rref_stepZ(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map);
    bool rref_step_numrows(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map);
    bool rref_step_numrowsV2(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map);

    // Put matrix into HNF Hermite Normal Form.
    //   Like RREF but using column operations. For solving initial constraints.
    //   "this" is A.  AU = (B 0) for column operations U.
    bool HNF(MatrixSparseInt &B, MatrixSparseInt &U, vector<int> &hnf_col_order, vector<int> &g);
    // use the HNF to solve Ax=b for x, i.e. satisfy the constraints but not the bounds
    bool HNF_satisfy_constraints(MatrixSparseInt &B, MatrixSparseInt &U, vector<int> &hnf_col_order, const vector<int> &g);
    // implementation, useful for other contexts
    static bool HNF_satisfy_rhs(MatrixSparseInt &B, MatrixSparseInt &U, vector<int> &hnf_col_order,
                                const vector<int> &rhs, vector<int> &col_solution, int col_size, const vector<int> &g);

    virtual void row_swap(int r, int s) override;
    virtual void row_negate(int r) override;
    // void col_swap(int c, int d); // this messes up the layout of IntVars and dummies, instead keep a new order in a vector

    virtual void allrows_eliminate_column(int r, int c) override; // with rhs

    virtual void row_simplify(int r) override;
    
    // get the greatest common divisor of matrix entries in this column, using only rows r and following
    int gcd_col(int r, int c);
    
    // replace row s with some multiple of row s and row r,
    //   so that in row s the column c entry is zero.
    // void row_eliminate_column(int r, int s, int c);
    // row[s] = u*s - v*r
    void row_us_minus_vr(int r, int s, int u, int v);
  
    // do the remaining stuff that the matrix versions don't
    void multiply_rhs(int r, int s, int u, int v);
    
    // The three main methods of the overall solve are
    //   satisfy_constraints:  initial feasible solution satisfying the constraints, Ax=b
    //   satisfy_bounds:       improve initial feasible solution so that all the variables are within their allowed ranges x \in [lo,hi]
    //   improve_solution:     improve the bounded-feasible solution so intervals are closer to the goal
    // satisfy_bounds and improve_solution work by finding a set of spanning vectors {x} for the nullspace of A, i.e. Ax=0, and using
    // those to move in a "downhill" direction
    
    //*** start: the following only work correctly after we've called rref
    
    // return true if a row is something like "0 = 3" or some other obviously bad constraint
    bool infeasible_constraints() const;

    // return true if calculated solution actually satisfies the constraints
    bool verify_constraints_satisfied(bool print=false);

    // find initial feasible solution, satisfying constraints
    // return true if successful
    bool satisfy_constraints(const vector<int>&cols_dep,
                            const vector<int>&cols_ind,
                            const vector<int>&rows_dep);
    
    // find initial feasible solution, satisfying constraints + variable bounds
    //   On input, solution already satisfies constraints,
    //   Use nullspace MN to increment intervals to satisfy bounds
    //     do Gaussian Elimination on M if blocked
    bool satisfy_bounds(const MatrixSparseInt&M, MatrixSparseInt&MN);

    // used by satisfy_bounds
    bool unbounded_bounds_improvement(int tc, int dx, map<int,int> &improved_cols,
                                      QWithReplacement &Q, SetValuesFn &val_bounds_Q, const double threshold,
                                      const MatrixSparseInt &MN,
                                      BlockingCols &blocking_cols, vector<int> &deferred_rows);
    bool any_bounds_improvement(int tc, int dx, map<int,int> &improved_cols,
                                QWithReplacement &Q, SetValuesFn &val_bounds_Q, const double threshold,
                                const MatrixSparseInt &MN, BlockingCols &blocking_cols, vector<int> &deferred_rows,
                                int &smallest_self_block_coeff, bool &give_up);
    

    // partition columns into dependent(leading) and independent(trailing) variables
    // returns true if some of the columns were swapped
    void categorize_vars(vector<int>&col_order,
                         vector<int>&cols_dep,
                         vector<int>&cols_ind,
                         vector<int>&rows_dep);
    
    // decide whether to round the solution of variable c up or down,
    //   for a variable with a non-unity coefficient in its row such that any integer value will not make the row feasible.
    //   based on bounds and goals
    bool decide_roundup(int r, int c) const;
    // assign vars to their closest integer goal
    //   only assigns INT_VARs
    //   if overwrite is false, then only variables who's value is 0 (uninitialized) are assigned.
    void assign_vars_goals(bool overwrite);
    // assign the dummy variables the closest value that would make the constraints feasible
    // returns true if the assignment is feasible
    void assign_dummies_feasible();
    // adjust independent vars to try to make the dependent ones feasible
    bool adjust_independent_vars(vector<int>&cols_fail);
    // if the variables weren't assignable as integers, then return false and save the unsatisfied rows in rows_fail
    bool assign_dependent_vars(const vector<int>&cols_dep,
                               const vector<int>&rows_dep,
                               vector<int>&cols_fail);
    // return sum of coeff*solution - rhs
    int row_sum(int r);
    int row_sum(const vector<int>&cols, const vector<int>&coeff, int rhs);
    static
    int row_sum(const vector<int>&cols, const vector<int>&coeff, const vector<int> &sol, int rhs); // implementation

    bool bounds_unsatisfied(int c);
    bool bounds_satisfied(vector<int>&cols_fail);
    
    // set M's rows to the nullspace basis vectors of this
    bool create_nullspace(const vector<int>&cols_dep,
                          const vector<int>&cols_ind,
                          const vector<int>&rows_dep,
                          MatrixSparseInt&M,
                          int &MaxMrow);
    
    // make a feasible solution better, using rows of nullspace matrix M.
    // Searches for a local min of the max deviation from goals.
    void improve_solution(const MatrixSparseInt&M, MatrixSparseInt&MN);
    
    // make a feasible solution better, using rows of nullspace matrix M.
    // call after improve_solution is done,
    // see if there is anything else we can do to improve the lexicographic solution, even if the max cannot be improved
    void fine_tune(const MatrixSparseInt&MN);
    
    //*** end: only works after we've called rref
    
    void identify_col_type();
    
    // convert inequalities to equalities
    //   if re-solving, add slack variables so we start with an initial feasible soution
    bool convert_inequalities();

    void count_used_rows_cols();

    void add_more_columns();
    void add_more_rows();
    
    // actually resize the vectors to number_of_rows and number_of_cols
    void resize_rows();
    void resize_cols();

    //== tied variables start
    // for the tied-to columns we keep this
    struct TiedDatum
    {
      int orig_bound_lo=0, orig_bound_hi=0;
      double goal_lo=0., goal_hi=0.;
    };
    map<int,TiedDatum> tied_data;
    
    // variables that are constrained to be identical by the matrix
    //   index of columns that are constrained to be identical
    //   one variable keeps a vector of all other variables that are identical to it; its column is retained in the matrix
    //   the other variables get a vector of that master variable; their columns entries are zeroed out
    bool has_tied_variables = false;
    vector< vector<int> > tied_variables;
    bool find_tied_variables();
    // replace variables with the ones they are tied to (i.e. constrained equal to)
    //  for each row, add up any identical-variable coefficients, fix rows and col_rows
    void replace_tied_variables();
    // replace tied variables with zero, e.g. for the nullspace
    void cull_tied_variables(MatrixSparseInt &M);
    // for tied variables, set the new col_solution to be the geometric average of the min and max tied variable
    bool should_adjust_solution_tied_variables = false;
    void adjust_solution_tied_variables();
    // remove all spurious elementary rows consisting of a single tied variable
    //   if tied variables, then the nullspace will have spurious elementary rows, one for each tied variable
    void cull_nullspace_tied_variable_rows(MatrixSparseInt&M, int &MaxMrow);
    // generate datum, returns false if the bounds are non-overlapping so infeasible
    bool generate_tied_data();
    // assign tied-variable to the same value as what they are tied-to.
    void assign_tied_variables();
    // checks if the variable is a tied-to variable, and uses the datum if so.
    // returns true if goal_lowest != goal_highest
    bool get_tied_goals(int c, double &goal_lowest, double &goal_highest) const; 
    //== tied variables end

    // amount we need to change the intervals to make it satisfy its variable bounds.
    // + = need to increase, - = need to decrease. Uses passed in value v rather than IIA's solution value
    int compute_out_of_bounds(int c, int v) const;
    
    // assign a quality to using row m to increment x[tc]+=dx
    //   criteria
    //     primary: keep row variables in bounds (or moving them closer to in bounds)
    //       return false if a row variable gets more out-of-bounds
    //     secondary: moving row variables towards their goals
    //   store an aggregate score for the row in rows_by_quality
    // return true if m was added to rows_by_quality
    // for choosing the best row, not just the first one that works.
    //   be sure to sort rows_by_quality after calling this on all candidate rows
    bool queue_rows_by_quality( int tc, int dx, const MatrixSparseInt &MN, int m,
                               vector< pair<int /*row*/, double /*quality*/> > &rows_by_quality);

    bool verify_full_solution(bool print_unsatisfied_constraints, bool pretend_solved = false) const;
    
    // increase the solution by dx times R. Returns true if the solution is (still) within bounds
    bool increment_solution(const RowSparseInt &R, int dx);

    // accept_strict_improvement
    //   if the solution strictly improves
    //      solution+=dx*mrow
    //      update Q
    //      return true
    //   else
    //      blocking_cols = vars that would get worse (accept self), return false
    
    // used by both satisy_bounds and improve_solution
    //   bounds: define "improvement" by the priority decreasing
    //   improve: define "improvement" by the quality function, and satisfying the constraint function
    bool accept_strict_improvement(QWithReplacement &Q,
                                   SetValuesFn &quality_fn,
                                   double quality_threshold,
                                   SetValuesFn *constraint_fn, // satisfying_bounds if this is null, "improving" otherwise
                                   double constraint_threshold,
                                   int tc,
                                   const RowSparseInt &R,
                                   int dx,
                                   bool &self_block,
                                   int &self_block_coeff,
                                   BlockingCols &blocking_cols,
                                   int &num_block,
                                   map<int,int> &improved_cols);

    bool is_improvement(QElement &s,
                        QElement &t,
                        SetValuesFn *constraint_fn,
                        double constraint_threshold ) const;

    // x > g ? x/g : g/x, but in a safe computation. c is column, the other values are retrieved
    double get_R(int c, int &x, double &glo, double &ghi) const;
    void compute_quality_ratio(vector<QElement> &q, const vector<int> &cols) const;
    
    // compare by lexicographic min max
    // return
    //   1 if A<B
    //  -1 if B<A
    //   0 if A==B
    // this IIA is only used for printing
    int is_better( vector<QElement> &qA, vector<QElement> &qB ) const;

    // get the top of the Q, after ensuring it is up to date
    // if true, Q is empty and the returned QElement is just the default, and shouldn't be processed
    bool tip_top( priority_queue<QElement> &Q, SetValuesFn &val_fn, double threshold, QElement &t );
    
    // build a priority_queue, containing qcol and/or qrow, with priority set by set_val_fn
    // queue elements with a value not less than threshold are not put on the queue
    void build_Q(priority_queue< QElement > &Q,
                 SetValuesFn &set_val_fn,
                 double threshold,
                 const vector<int> &qcol) const;
    
  protected:
    // priorities for selecting variables to change.
    friend class SetValuesFn;
    friend class SetValuesBounds;
    friend class SetValuesRatioR;
    friend class SetValuesRatioG;
    friend class SetValuesRatio;
    friend class SetValuesOneCoeff;
    friend class SetValuesNumCoeff;
    friend class SetValuesNumCoeffV2;
    friend class SetValuesCoeffRowsGoal;
    friend class SetValuesCoeffRGoal;
    friend class SetValuesTiny;
  };
  
//  inline
//  IncrementalIntervalAssignment::IncrementalIntervalAssignment( IAResultImplementation *result_ptr ) : MatrixSparseInt(result_ptr), Nullspace(result_ptr)
//  {
//  }
  
  inline
  int IncrementalIntervalAssignment::get_B( int row ) const
  {
    return rhs[row];
  }
  
  inline
  void IncrementalIntervalAssignment::set_B( int row, int val )
  {
    if (rhs[row] != val)
    {
      hasBeenSolved=false;
      if ( row <= solved_used_row )
      {
        solved_used_row = -1;
        solved_used_col = -1;
      }
    }

    rhs[row] = val;
  }
  
  inline
  double IncrementalIntervalAssignment::get_goal(int c) const
  {
    return goals[c];
  }
  
  inline
  void IncrementalIntervalAssignment::set_no_goal(int c)
  {
    if (goals[c] != 0.)
    {
      hasBeenSolved=false;
      if ( c <= solved_used_col )
      {
        solved_used_row = -1;
        solved_used_col = -1;
      }
    }

    goals[c]=0.; // flag that the goal should be ignored
  }
  
  inline
  int IncrementalIntervalAssignment::get_M  ( int row, int col ) const
  {
    return rows[row].get_val(col);
  }
  
  inline
  void IncrementalIntervalAssignment::get_M( int row, const vector<int> *&cols, const vector<int> *&vals ) const
  {
    cols = &rows[row].cols;
    vals = &rows[row].vals;
  }

  inline
  void IncrementalIntervalAssignment::clear_M(int row)
  {
    rows[row].cols.clear();
    rows[row].vals.clear();
  }
  
  inline
  int IncrementalIntervalAssignment::row_sum(const vector<int>&cols,const vector<int>&coeff, int rhs)
  {
    return row_sum(cols,coeff,col_solution,rhs);
  }
  
  inline
  int IncrementalIntervalAssignment::row_sum(int r)
  {
    return row_sum(rows[r].cols,rows[r].vals,rhs[r]);
  }

  inline
  int IncrementalIntervalAssignment::get_is_solved() const
  {
    return hasBeenSolved;
    
  }
  
  inline
  void IncrementalIntervalAssignment::set_is_unsolved()
  {
    hasBeenSolved=false;
    solved_used_row = -1;
    solved_used_col = -1;
  }

  inline
  void IncrementalIntervalAssignment::set_constraint(int row, ConstraintType constraint_type)
  {
    if (constraint[row]!=constraint_type)
    {
      hasBeenSolved=false;
      if (row <= solved_used_row)
      {
        solved_used_row=-1;
        solved_used_col=-1;
      }
    }
    
    constraint[row]=constraint_type;
  }

  inline
  void IncrementalIntervalAssignment::set_lower ( int col, int bound )
  {
    col_lower_bounds[col]=bound;
    
    if (col_solution[col]<bound)
    {
      hasBeenSolved=false;
      if ( col <= solved_used_col )
      {
        solved_used_row = -1;
        solved_used_col = -1;
      }
    }
  }
  
  inline
  void IncrementalIntervalAssignment::set_upper ( int col, int bound )
  {
    col_upper_bounds[col]=bound;

    if (col_solution[col]>bound)
    {
      hasBeenSolved=false;
      if ( col <= solved_used_col )
      {
        solved_used_row = -1;
        solved_used_col = -1;
      }
    }
  }

  inline
  void RowSparseInt::multiply(int k)
  {
    for (auto &v : vals)
      v*=k;
  }

inline
void MatrixSparseInt::row_negate(int r)
{
  for (auto &v : rows[r].vals)
  {
    v = -v;
  }
}

inline
void IncrementalIntervalAssignment::row_negate(int r)
{
  MatrixSparseInt::row_negate(r);
  rhs[r]=-rhs[r];
  if (constraint[r] == IIA::LE)
    constraint[r] = IIA::GE;
  else if (constraint[r] == IIA::GE)
    constraint[r] = IIA::LE;
}

inline
void MatrixSparseInt::allrows_eliminate_column(int r, int c)
{
  matrix_allrows_eliminate_column(r, c, nullptr);
}
inline
void IncrementalIntervalAssignment::allrows_eliminate_column(int r, int c)
{
  matrix_allrows_eliminate_column(r, c, &rhs);
}


} // namespace

#endif  // INCREMENTALINTERVALASSIGNMENT_HPP
