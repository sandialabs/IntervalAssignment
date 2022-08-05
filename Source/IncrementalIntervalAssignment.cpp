//- file IncrementalIntervalAssignment.cpp
#include "IncrementalIntervalAssignment.h"

#include "IAQueue.h"

#include <numeric>
#include <cmath>
#include <cassert>

#include "IAResultImplementation.h"

#include "CpuTimer.h"

// math utilities
namespace IIA_Internal
{
  
  using IIA::ConstraintType;
  using IIA::EQ;
  using IIA::LE;
  using IIA::GE;
  using IIA::EVEN;
  
  using std::max;
  using std::min;
  using std::lround;
  using std::copy;
  using std::swap;
  using std::sort;
  using std::make_pair;
  using std::to_string;
  using std::iota;
  using std::lower_bound;
  using std::back_inserter;

  // forward declarations of internal methods, for reference
  int round_goal(double g); // return the value of x that minimizes R
  int gcd(int u, int v);
  int lcm(int u, int v);
  void row_simplify_vals(vector<int> &v); // divide row vals by gcd
  bool row_equal(const RowSparseInt &r, const RowSparseInt &s); // true if s = r or s = -r
  // may want a method that determines if s = m*r for integer m != 1

  // implementations
  int round_goal(double g)
  {
    const int lo = (int) g;
    if (lo != g)
    {
      auto hi = lo+1;
      auto gsame = sqrt(lo*hi);
      return (g>gsame) ? hi : lo;
    }
    else
    {
      return lo;
    }
  }

  // greatest common divisor of u and v
  // always returns a positive number, as we take abs of u and v
  int gcd(int u, int v)
  {
    // if C++17, return gcd
    
    // "The Euclidean Algorithm"
    auto p = max(abs(u),abs(v));
    auto q = min(abs(u),abs(v));
    while (true)
    {
      auto r = p % q;
      if ( r < 1 )
        break;
      p=q;
      q=r;
    }
    return q;
  }
  
  // least common multiple
  // u * v = GCD * LCM
  int lcm(int u, int v)
  {
    // ordered this way to avoid overflow, keep numbers small
    // (u/gcd)*v == (u*v)/gcduv, where we know (u % gcduv) == 0
    const auto gcduv = gcd(u, v);
    return gcduv ? (u / gcduv * v) : 0;
  }
  
  // lcm of submatrix
  // lcm_dep = "dep_coeff[0];"
  // for ( all other dep_coeff )
  // lcm_dep = lcm( lcm_dep, dep_coeff[i])
  
  int IncrementalIntervalAssignment::gcd_col(int r, int c)
  {
    auto &cr=col_rows[c];
    // assert(!cr.empty()); it can be empty, return -1 in that case
    int g=-1;
    for (auto rr : cr)
    {
      if (rr>=r)
      {
        if (g==-1)
          g =     get_M(rr,c);
        else
          g = gcd(get_M(rr,c),g);
        if (g==1)
          break;
      }
    }
    return g;
  }
  
  // return the divisor
  void row_simplify_vals(vector<int> &v)
  {
    if (v.empty())
      return;
    
    // get gcd of the entire row
    int g = abs(v[0]); // leading coefficient
    for (size_t c = 1; c < v.size(); ++c)
    {
      assert(v[c]!=0);
      g = gcd(v[c],g);
      if (g==1)
        break;
    }
    
    if (g==1)
    {
      // negate?
      if (v[0] < 0)
      {
        for (auto &vv : v)
        {
          vv = -vv;
        }
      }
    }
    // row division
    else
    {
      for (auto &vv : v)
      {
        assert(vv % g == 0);
        vv /= g;
      }
    }
    // return m;
  }
  
  bool row_equal(const RowSparseInt &r, const RowSparseInt &s)
  {
    // different number of non-zeros?
    if (r.cols.size() != s.cols.size())
      return false;
    // both empty
    if (r.cols.empty())
      return true;
    const bool is_negated = (r.vals[0] == -s.vals[0]);
    for (int i = 0; i < (int) r.cols.size(); ++i)
    {
      // different variable column?
      if (r.cols[i] != s.cols[i])
        return false;
      // different coeff?
      if (( is_negated && r.vals[0] != -s.vals[0]) ||
          (!is_negated && r.vals[0] !=  s.vals[0]))
      {
        return false;
      }
    }
    return true;
  }

  MatrixSparseInt::MatrixSparseInt(const MatrixSparseInt&M, int MaxMrow) : result(M.result)
  {
    rows.reserve(MaxMrow);
    copy(M.rows.begin(),M.rows.begin()+MaxMrow,back_inserter(rows));
    col_rows.resize(M.col_rows.size());
    fill_in_cols_from_rows();
  }
  
  MatrixSparseInt::MatrixSparseInt(const MatrixSparseInt&M, int MaxMrow, const MatrixSparseInt&N, bool no_duplicate_rows) : result(M.result)
  {
    rows.reserve(MaxMrow+N.rows.size());
    copy(N.rows.begin(),N.rows.end(),back_inserter(rows));
    if (no_duplicate_rows && !N.rows.empty())
    {
      for (int r=0; r<MaxMrow; ++r)
      {
        // does row r of M exist in N?
        auto &mrow = M.rows[r];
        assert(!mrow.cols.empty());
        int lead_col = mrow.cols[0];
        // check rows in N with same leading col
        bool found_duplicate_row = false;
        for (int s : N.col_rows[lead_col])
        {
          if (row_equal(mrow, N.rows[s]))
          {
            found_duplicate_row = true;
            break;
          }
        }
        if (!found_duplicate_row)
        {
          rows.push_back(mrow);
        }
      }
    }
    else
    {
      copy(M.rows.begin(),M.rows.begin()+MaxMrow,back_inserter(rows));
    }
    col_rows.resize(M.col_rows.size());
    fill_in_cols_from_rows();
  }

  int IncrementalIntervalAssignment::row_sum(const vector<int>&cols, const vector<int>&coeff, const vector<int> &sol, int rhs)
  {
    int sum=0;
    assert(cols.size()==coeff.size());
    for (size_t i = 0; i < cols.size(); ++i)
    {
      assert(cols[i]<(int)sol.size());
      sum += coeff[i]*sol[cols[i]];
    }
    sum -= rhs;
    return sum;
  }
  
  // assign independent vars to their closest integer goal
  void IncrementalIntervalAssignment::assign_vars_goals(bool overwrite)
  {
    for ( int i = 0; i <= used_col; ++i)
    {
      if (overwrite || col_solution[i]==0) // zzyk need some other flag than "0" here to indicate unassigned variables, as a sum-even var in a small loop might have been assigned 0 on purpose
      {
        int v = round_goal( goals[i] );
        v = max(v,col_lower_bounds[i]);
        v = min(v,col_upper_bounds[i]);
        col_solution[i] = v;
      }
    }
  }
  
  bool IncrementalIntervalAssignment::decide_roundup(int r, int c) const
  {
    // bounds
    const auto oob_lo = compute_out_of_bounds(c, col_solution[c]);
    const auto oob_hi = compute_out_of_bounds(c, col_solution[c]+1);
    if (oob_lo>0)
      return true;
    else if (oob_hi<0)
      return false;
    else
    {
      // can round up or down and still be in bounds
      
      const auto coeff = rows[r].get_val(c);
      
      // secondary criteria:
      //   goals vs. solution of vars in row[c]
      auto signed_diff_sum = 0.; // sum of solution-goal
      for (int j = 0; j < (int) rows[r].cols.size(); ++j)
      {
        const auto cc = rows[r].cols[j];
        const auto vv = rows[r].vals[j];
        const auto s = get_solution(cc);

        const auto opposite_actions = (cc != c && coeff>0 == vv>0);
        
        const auto oob_cc_up = compute_out_of_bounds(cc, s+1);
        const auto oob_cc_dn = compute_out_of_bounds(cc, s-1);
        
        if (oob_cc_dn>0) // we want cc to increase (or not decrease) to stay in bounds
        {
          if (c!=cc)
          {
            if (opposite_actions)
              signed_diff_sum -= oob_cc_dn; // subtract positive value -> decreases desire of roundup
            else
              signed_diff_sum += oob_cc_dn; // add positive value -> increases
          }
        }
        else if (oob_cc_up<0) // we want cc to decrease (or not increase) to stay in bounds
        {
          if (c!=cc)
          {
            if (opposite_actions)
              signed_diff_sum -= oob_cc_dn; // subtract negative value -> increases
            else
              signed_diff_sum += oob_cc_dn; // add negative value -> decreases
          }
        }
        else // if (oob_cc_up == 0 && oob_cc_dn == 0)
        {
          double glo, ghi;
          get_tied_goals(cc, glo, ghi);
          const auto g = (glo==ghi ? glo : sqrt(glo*ghi));
          if (g!=0.)
          {
            const auto signed_diff = s-g;
            
            // does rounding up encourage cc to increase or decrease?
            if (opposite_actions)
            {
              // rounding up causes cc to decrease.
              // If signed_diff>0, then rounding up is good
              signed_diff_sum += signed_diff;
            }
            else
            {
              // rounding up causes cc to increase.
              // If signed_diff<0, then rounding up is good
              signed_diff_sum -= signed_diff;
            }
            
          }
        }
      }
      return (signed_diff_sum >= -3.1); // empirical threshold
    }
  }
  
  // d = a/b, return true if a%b != 0
  // always round d down for fractional d
  // e.g. -7/2 = -4, and 7/2 = 3
  //   platform-independent results for negative values by avoiding negative values during division and mod
  bool a_div_b( int &d, int a, int b)
  {
    assert(b!=0);
    d = abs(a) / abs(b); // always positive, rounded towards zero
    const bool is_neg = ((a>=0) != (b>0)); // true if (a/b)<0 <--> a*b<0 <--> sign(a) != sign(b)
    const bool is_frac = (abs(a) % abs(b) != 0);
    if (is_neg)
    {
      d = -d;
      if (is_frac)
      {
        --d;
      }
    }
    return is_frac;
  }
  
  void test_a_div_b()
  {
    int d;
    for (int a=-9; a<10; ++a)
      for( int b=-9; b<10; ++b)
      {
        if (b==0)
          continue;
        const bool is_frac = a_div_b( d, a, b);
        printf("%d = %d / %d is %s.\n", d, a, b, (is_frac ? "fractional" : "whole"));
      }
  }
  
  bool IncrementalIntervalAssignment::assign_dependent_vars(const vector<int>&cols_dep,
                                                            const vector<int>&rows_dep,
                                                            vector<int>&cols_fail)
  {
    // test_a_div_b();
    
    assert(cols_dep.size()==rows_dep.size());
    bool OK=true;
    for (size_t i = 0; i < rows_dep.size(); ++i)
    {
      const auto r = rows_dep[i];
      const auto c = cols_dep[i];
      col_solution[c]=0; // to avoid it contributing to the row_sum
      const int sum = row_sum(r);
      const int coeff = rows[r].get_val(c);
      //
      //   col_solution[c] = -sum / coeff;
      if (sum == 0)
      {
        col_solution[c] = 0;
      }
      else
      {
        assert(coeff!=0);
        const bool is_frac = a_div_b( col_solution[c], -sum, coeff ); // solution = -sum / coeff
        if (is_frac)
        {
          // we have a problem, the dependent var value is not integer (its coeff_dep is not 1)
          OK=false;
          cols_fail.push_back(c);
          
          // decide whether to round up or down
          bool roundup = true;
          
          if (use_smart_adjust_for_fractions)
          {
            roundup = decide_roundup(r,c);
            result->debug_message("IIA dependent var x%d = %d/%d != integer. Rounding %s.\n",
                                 c, -sum, coeff, roundup ? "up" : "down");
          }
          
          // round solution value up or down
          if (roundup)
            ++col_solution[c];
          // else already rounded down, no need for --sol
          
          result->debug_message("IIA dependent var x%d = %d/%d != integer. Rounding up.\n",
                                c, -sum, coeff, roundup ? "up" : "down");
        }
      }
      result->debug_message("IIA Assigning dependent var x%d = %d/%d = %d%s\n",
                            c, -sum, coeff, col_solution[c],
                            (col_solution[c]<col_lower_bounds[c] || col_solution[c]>col_solution[c]) ? " out of bounds!" : "");
    }
    if (result->log_debug)
    {
      print_solution("after assigning dependent variables");
      if (!cols_fail.empty())
      {
        result->info_message("failed-to-assign variables, cols_fail A");
        print_vec(result, cols_fail);
      }
    }
    
    return OK;
  }
  
  void IncrementalIntervalAssignment::assign_dummies_feasible()
  {
    // solution_init has already been called
    // already assigned intVars to lround(goals), or to some other value from a prior mapping subproblem
    // get dummies close to what makes sense for the goals
    // i.e. initilize sum-even variables to row-sum / 2
    
    for (int c = 0; c <= used_col; ++c)
    {
      if (col_type[c]!=INT_VAR)
      {
        
        auto &dummy_rows = col_rows[c];
        if (dummy_rows.empty())
        {
          result->debug_message("dummy x%d adjustment skipped because it isn't in any row\n", c);
          continue;
        }

        // we only assign it based on the first row its in
        // if there is more than one dummy variable with different coefficients, we could take a second pass to try to satisfy it
        //   fill in if it helps a customer problem

        int old_val = col_solution[c];
        
        // assign col_solution[c]
        int r = dummy_rows[0];
        col_solution[c]=0;
        int sum = row_sum(r);
        int dummy_coeff = get_M(r,c);
        
        const bool is_frac = a_div_b( col_solution[c], -sum, dummy_coeff ); // solution = -sum / coeff
        if (is_frac)
        {
          bool roundup = true;
          if (use_smart_adjust_for_fractions)
          {
            roundup = decide_roundup(r, c);
            result->debug_message("constraint row %d not satisfied, dummy var %d rounded %s.\n",
                                 r, c, roundup ? "up" : "down");
          }
          if (roundup)
            ++col_solution[c];
          
          // could try adding in a mapping-phase nullspace vec here
        }
        if (result->log_debug)
        {
          print_row_iia(r);
          result->debug_message("dummy x%d changed from %d to %d based on row %d sum %d\n", c, old_val, col_solution[c], r, sum);
        }
      }
    }
  }
  

  bool IncrementalIntervalAssignment::find_tied_variables()
  {
    // after calling
    //  variables c that are unaffected have empty tied_variables[c]
    //            c          to be replaced by c2   tied_variables[c]={c2}
    //            c2         that are the replacement of a,b,...d  have   tied_variables[c]={c2,a,b,...d}

    assert(tied_variables.empty()); // otherwise need to clear out the data
    tied_variables.resize(number_of_cols);
    
    // to disable tied variables, uncomment the following return statement
    // return false;
    
    bool found=false;
    for (int r = 0; r <= used_row; ++r)
    {
      // here is the check as to whether variables are tied or not
      //   if we support non-equality constraints, then we should check that the constraint is an equality
      //   but currently we've already converted all inequalities to equalities by the time we've gotten here
      if (rows[r].cols.size()==2 && rhs[r]==0.) // if tied but not exactly equal, we don't handle that yet
      {
        // check coefficients, currently only val and -val are detected.
        //   Could extend to val and val, with sign bookkeeping, but that doesn't happen AFAIK.
        if (rows[r].vals[0] != -rows[r].vals[1])
          continue;
        found=true;
        
        // which col to use
        int c1 = rows[r].cols[0];
        int c2 = rows[r].cols[1];
        
        // default to c1 being the smaller index
        if (c1>c2)
        {
          swap(c1,c2);
        }
        // replace c1 and c2 with what they're tied to, if anything
        if (tied_variables[c1].size()==1)
        {
          c1 = tied_variables[c1].front();
        }
        if (tied_variables[c2].size()==1)
        {
          c2 = tied_variables[c2].front();
        }
        // if they are already tied to the same variable, do nothing
        if (c1==c2)
          continue;
        // make c1 the bigger set, for efficiency
        if (tied_variables[c1].size() < tied_variables[c2].size())
        {
          swap(c1,c2);
        }
        // c2 is not currently tied
        //   tie c2 to c1
        if (tied_variables[c2].empty())
        {
          tied_variables[c2].push_back(c1);
          // if c1 is not tied to anything yet, put itself as the first entry
          if (tied_variables[c1].empty())
            tied_variables[c1].push_back(c1);
          tied_variables[c1].push_back(c2);
        }
        // both c1 and c2 are tied, need to merge
        else
        {
          assert(c1!=c2);
          assert( tied_variables[c1].size()>1 );
          assert( tied_variables[c2].size()>1 );
          // vars are tied to c1, and tied to c2, merge
          for (auto c : tied_variables[c2])
          {
            tied_variables[c1].push_back(c);
            assert( tied_variables[c].size()==1 || c==c2);
            tied_variables[c].front()=c1;
          }
          tied_variables[c2].resize(1);
          assert(tied_variables[c2].front()==c1);
        }
      }
    }
    has_tied_variables=found;
    
    // debug
    // print which variables are tied to which
    if (result->log_debug)
    {
      int num_tied=0, num_tied_to=0;
      for (int c=0; (size_t) c<tied_variables.size(); ++c)
      {
        if (tied_variables[c].size()==1)
          ++num_tied;
        else if (tied_variables[c].size()>1)
        {
          ++num_tied_to;
        }
      }
      result->info_message("%d variables were removed and tied to %d variables\n", num_tied, num_tied_to);
      for (int c=0; (size_t) c<tied_variables.size(); ++c)
      {
        if (!tied_variables[c].empty())
        {
          result->info_message("%d tied :",c);
          for (auto c2 : tied_variables[c])
          {
            result->info_message(" %d",c2);
          }
          result->info_message("\n");
        }
      }
    }
    // debug
    
    return has_tied_variables;
  }

  // return in tied_var_rows all the rows of M that contain a tied variable column
  bool tied_variable_rows(const MatrixSparseInt &M, const vector< vector<int> > &tied_variables, set<int> &tied_var_rows)
  {
    for (int c=0; (size_t) c<tied_variables.size(); ++c)
    {
      if (tied_variables[c].size()==1)
      {
        // we need to replace c with tied_variables[c].front() in the matrix, in all the rows that c appears in
        assert(tied_variables[c].front()!=c);
        tied_var_rows.insert(M.col_rows[c].begin(),M.col_rows[c].end());
      }
    }
    return (!tied_var_rows.empty());
  }

  void IncrementalIntervalAssignment::replace_tied_variables()
  {
    // gather affected rows
    set<int> row_needs_sub;
    if (!tied_variable_rows(*this,tied_variables,row_needs_sub))
      return;

    // replace the coefficient of a tied variable with a coefficient of what it is tied to
    for (auto r : row_needs_sub)
    {
      auto &row = rows[r];
      map<int,int> rowents; // new row values
      for (size_t i=0; i<row.cols.size(); ++i)
      {
        int c = row.cols[i];
        if (tied_variables[c].size()==1)
        {
          c = tied_variables[c].front();
        }
        rowents[c]+=row.vals[i]; // initializes to 0 if doesn't exist yet
      }
      row.cols.clear();
      row.vals.clear();
      for (auto m : rowents)
      {
        if (m.second!=0)
        {
          row.cols.push_back(m.first);
          row.vals.push_back(m.second);
        }
      }
      // row is sorted because rowents was sorted
    }
    // efficiency improvement: could update the col_rows. For now, just rebuild the col_rows from scratch.
    fill_in_cols_from_rows();
    // this can result in an empty matrix, e.g. from the opposite sides of a single mapping face adjacent to a different mapping face
  }
  
  void IncrementalIntervalAssignment::cull_tied_variables(MatrixSparseInt &M)
  {
    
    // gather affected rows
    set<int> row_needs_sub;
    if (!tied_variable_rows(M,tied_variables,row_needs_sub))
      return;
    
    // remove tied variables from rows
    for (auto r : row_needs_sub)
    {
      RowSparseInt newrow;
      auto &row = M.rows[r];

      for (size_t i=0; i<row.cols.size(); ++i)
      {
        int c = row.cols[i];
        if (tied_variables[c].size()!=1)
        {
          newrow.cols.push_back(c);
          newrow.vals.push_back(row.vals[i]);
        }
      }
      swap(row,newrow);
    }
    // remove rows from tied variables
    for (int c=0; (size_t) c<tied_variables.size(); ++c)
    {
      if (tied_variables[c].size()==1)
        M.col_rows[c].clear();
    }
  }

  void IncrementalIntervalAssignment::adjust_solution_tied_variables()
  {
    for (int c = 0; (size_t) c < tied_variables.size(); ++c)
    {
      // adjust the solution of the tied-to variables only
      if (tied_variables[c].size()>1)
      {
        auto old_val = col_solution[c];
        
        auto &td = tied_data.at(c); // should already exist
        
        // assign new solution to be the geometric mean
        auto &lo = td.goal_lo; // shorthand
        auto &hi = td.goal_hi;
        if (lo!=hi)
        {
          double x = sqrt( (double) lo * hi );
          // figure out whether rounding x up or down minimizes the maximum of x/lo and hi/x
          int s = floor(x);
          int t = s+1;
          double xx = sqrt(s*t);
          col_solution[c] = (x>xx) ? t : s;
        }
        if (col_solution[c]<col_lower_bounds[c])
          col_solution[c]=col_lower_bounds[c];
        else if (col_solution[c]>col_upper_bounds[c])
          col_solution[c] = col_upper_bounds[c];
        
        // also assign solution to the tied variables
        
        
        // debug
        if (result->log_debug)
        {
          if (old_val != col_solution[c])
          {
            result->debug_message("MasterTie x%d adjusted from %d to %d\n", c, old_val, col_solution[c]);
          }
        }
        // debug
        
      }
    }
  }
  
  void IncrementalIntervalAssignment::cull_nullspace_tied_variable_rows(MatrixSparseInt&M, int &MaxMrow)
  {
    bool found=false;
    for (int r=0; r<MaxMrow; )
    {
      if (M.rows[r].cols.size()==1 &&
          tied_variables[ M.rows[r].cols.front() ].size()==1) // c = M.rows[r].cols.front()
      {
        // clear spurious row, and move to last place
        assert(MaxMrow);
        M.rows[r].vals.clear();
        M.rows[r].cols.clear();
        --MaxMrow;
        M.rows[r].cols.swap(M.rows[MaxMrow].cols);
        M.rows[r].vals.swap(M.rows[MaxMrow].vals);
        // col_rows are now invalid, because we swapped rows, and delete a row
        //   clear col_rows for tied variable is done below in fill_in_cols_from_rows
        //   alternative is to call M.swap_rows which updates col_rows as we go
        
        found=true;
        // need to check the new row we swapped into r, don't increment r
      }
      else
      {
        ++r;
      }
    }
    // rebuild col_rows
    if (found)
    {
      // assert(MaxMrow); // it's possible that the nullspace is reduced to nothing now
      M.rows.resize(MaxMrow);
      M.fill_in_cols_from_rows();
    }
  }
  
  void IncrementalIntervalAssignment::assign_tied_variables()
  {
    for (int c=0; (size_t) c<tied_variables.size(); ++c)
    {
      if (tied_variables[c].size()==1)
      {
        int c2 = tied_variables[c].front();
        assert(c2!=c);
        col_solution[c]=col_solution[c2];
      }
    }
  }
  
  int IncrementalIntervalAssignment::compute_out_of_bounds(int c, int v) const
  {
    // int v = col_solution[c];
    if (v<col_lower_bounds[c])
    {
      return col_lower_bounds[c]-v; // >0
    }
    if (v>col_upper_bounds[c])
    {
      return col_upper_bounds[c]-v; // <0
    }
    return 0;
  }
  
  bool IncrementalIntervalAssignment::get_tied_goals(int c, double &goal_lowest, double &goal_highest) const
  {
    if (has_tied_variables && tied_variables[c].size()>1)
    {
      auto &dt = tied_data.at(c);
      goal_lowest=dt.goal_lo;
      goal_highest=dt.goal_hi;
      return true; // (goal_highest != goal_lowest);
    }
    else
    {
      goal_lowest = goals[c];
      goal_highest = goal_lowest;
      return false;
    }
  }
  
  double IncrementalIntervalAssignment::set_goal(int c, double goal)
  {
    
    // scale goals to range [t2,infinity) to avoid blow-up elsewhere
    //   monotonic, continuous
    const double t1=1.e-3;
    const double t2=1.e-5;
    if (goal<t1)
    {
      // scale 0 <= goals < t1 to [2*t2, t1]
      //   linear function with g(0)=2*t2, and g(t1)=t1
      // scale goal < 0 to [t2, 2*t2]
      goal = (goal<0. ? t2*(1. + 1./(1.-goal))
              : 2.*t2 + ((t1 - 2.*t2)/t1) *goal);
    }
    
    if (goals[c] != goal)
    {
      hasBeenSolved=false;
      if ( c <= solved_used_col )
      {
        solved_used_row = -1;
        solved_used_col = -1;
      }
    }
    
    goals[c]=goal;
    
    return goal;
  }
  
  bool IncrementalIntervalAssignment::infeasible_constraints() const
  {
    bool rc=false;
    for (int r = 0; r <= used_row; ++r)
    {
      auto &row = rows[r];
      const int b = rhs[r];
      const auto con = constraint[r];
      if (row.cols.empty() && b!=0)
      {
        // 0 = constant
        bool bad = (con==EQ || (con==GE && b>0) || (con==LE && b<0));
        if (bad)
        {
          if (result->log_info || result->log_debug)
          {
            rc=true;
            result->info_message("Interval matching row %d reduced to 0 == %d, problem is overconstrained and infeasible.\n", r, b);
            if (result->log_debug)
              print_row_iia(r);
          }
          else
            return true;
        }
      }
      // only one variable in the row
      else if (row.cols.size()==1 && con==EQ)
      {
        int v = row.vals[0];
        // solution not integer.  x = b/v, but v does not divide b
        if (b % v)
        {
          if (result->log_info || result->log_debug)
          {
            rc=true;
            result->info_message("Interval matching row %d has no integer solution, problem is infeasible.\n",r);
            if (result->log_debug)
              print_row_iia(r);
          }
          else
            return true;
        }
        else
        {
          int c = row.cols[0];
          int s = b / v;
          // col_solution[c]=s;
          // assume constraint is equality, we're inside IncrementalIntervalAssignment after all
          
          // only solution is below lower bound
          if (s < col_lower_bounds[c])
          {
            if (result->log_info || result->log_debug)
            {
              rc=true;
              result->info_message("Interval matching row %d, its only solution is below the variable lower bound, problem is infeasible.\n",r);
              result->info_message("x%d only has solution %d < %d lower bound.\n", c, s, col_lower_bounds[c]);
            }
            else
              return true;
          }
          // only solution is above upper bound
          else if (s > col_upper_bounds[c])
          {
            if (result->log_info || result->log_debug)
            {
              rc=true;
              result->info_message("Interval matching row %d, its only solution is above the variable upper bound, problem is infeasible.\n",r);
              result->info_message("x%d only has solution %d > %d upper bound", c, s, col_upper_bounds[c]);
            }
            else
              return true;
          }
        }
      }
      else if (constraint[r]==LE || constraint[r]==GE)
      {
        // check if constraint is feasible first, return failure if not
        if (con==LE)
        {
          int sum=0;
          // what's the smallest the row sum could possibly be
          for (size_t i = 0; i < row.cols.size(); ++i)
          {
            auto c = row.cols[i];
            auto v = row.vals[i];
            int m;
            if (v<0)
            {
              m = col_upper_bounds[c];
              if (m==numeric_limits<int>::max())
              {
                sum = numeric_limits<int>::lowest(); // sum is unbounded below
                break;
              }
            }
            else
            {
              assert(v>0);
              m = col_lower_bounds[c];
              if (m==numeric_limits<int>::lowest())
              {
                sum = numeric_limits<int>::lowest(); // sum is unbounded below
                break;
              }
            }
            sum += v*m;
          }
          if (sum>b)
          {
            if (result->log_info || result->log_debug)
            {
              rc=true;
              result->info_message("Interval matching row %d sum must be <= %d, but smallest possible sum is %d; problem is infeasible.\n", r, b, sum);
              if (result->log_debug)
                print_row_iia(r);
            }
            else
              return true;
          }
        }
        else if (con==GE)
        {
          int sum=0;
          // what's the largest the row sum could possibly be
          for (size_t i = 0; i < row.cols.size(); ++i)
          {
            auto c = row.cols[i];
            auto v = row.vals[i];
            int m;
            if (v>0)
            {
              m = col_upper_bounds[c];
              if (m==numeric_limits<int>::max())
              {
                sum = numeric_limits<int>::max(); // sum is unbounded above
                break;
              }
            }
            else
            {
              assert(v<0);
              m = col_lower_bounds[c];
              if (m==numeric_limits<int>::lowest())
              {
                sum = numeric_limits<int>::max(); // sum is unbounded above
                break;
              }
            }
            sum += v*m;
          }
          if (sum<b)
          {
            if (result->log_info || result->log_debug)
            {
              rc=true;
              result->info_message("Interval matching row %d sum must be >= %d, but largest possible sum is %d; problem is infeasible.\n", r, b, sum);
              if (result->log_debug)
                print_row_iia(r);
            }
            else
              return true;
            
          }
        }
      }
    }
    return rc;
  }
  
  
  bool IncrementalIntervalAssignment::satisfy_constraints(const vector<int>&cols_dep,
                                                         const vector<int>&cols_ind,
                                                         const vector<int>&rows_dep)
  {
    // debug
    CpuTimer *timer(nullptr);
    if (result->log_debug_time)
    {
      timer = new CpuTimer;
    }
    result->debug_message("\nIIA: finding a solution that satisfies the constraints\n");
    // debug
    
    if (infeasible_constraints())
    {
      if (result->log_info)  result->error_message(    "Match Intervals problem is overconstrained and infeasible; there is no solution\n");
      else            result->debug_message("Match Intervals problem is overconstrained and infeasible; there is no solution\n");
      return false;
    }
    
    // try to satisfy constraints using direct back-substitution using the rref
    //   if it works, gets a solution close to the goals
    // backup is HNF, which always finds a solution if it exists, but the solution is often far from the goals
    
    bool rc=false;
    bool did_hnf=false;
    
    //rref
    // independent vars have already been assigned something
    if (!use_HNF_always)
    {
      // back-substitute values from indepedent vars to dependent vars
      rc=true; // will set it to false if it fails
      vector<int>cols_fail;
      if (!assign_dependent_vars(cols_dep,rows_dep,cols_fail))
      {
        // result->info_message("cols_fail B"); print_vec(result, cols_fail);
        result->debug_message("First pass at assigning dependent vars integer values failed. Now trying to adjust independent vars.\n");
        // work harder for a feasible assignment
        
        // adjust independent vars that only appear in a single row, see if that provides an integer solution to the dependent vars
        if (!adjust_independent_vars(cols_fail))
        {
          result->debug_message("Second pass at assigning dependent vars integer values failed.\n");
          rc=false;
        }
      }
    }
    //HNF
    if (!rc)
    {
      // fallback is HNF
      //   It always finds an integer solution if there is one
      result->debug_message("Now trying Hermite normal form. If that fails, the problem is infeasible and there is no solution.\n");
      
      did_hnf=true;
      MatrixSparseInt B(result), U(result);
      vector<int> hnf_col_order;
      // init g to col_solution, because col_solution was previously assigned based on goals and dummy-variable's row-sums
      vector<int> g(col_solution);
      auto OK = HNF(B,U,hnf_col_order,g);
      if (!OK)
      {
        result->error_message("HNF algorithm failed; implementation error.\n");
      }
      else
      {
        OK = HNF_satisfy_constraints(B,U,hnf_col_order,g);
        if (OK)
        {
          rc=true;
        }
        // else
        //   problem is infeasible, it's not an error if this happens
      }
    }
    
    // verify solution
    if (rc && !did_hnf)
    {
      verify_constraints_satisfied(result->log_info || result->log_debug); // hnf_satisfy_constraints calls this already itself
    }
    
    // debug
    if (result->log_debug_time)
    {
      double time_feas = timer->cpu_secs();
      result->info_message("Satisfy constraint time %f\n",time_feas);
      delete timer;
    }
    else if (result->log_debug)
    {
      result->debug_message("Satisfy constraint done\n");
    }
    //debug
    
    return rc;
  }
  
  void IncrementalIntervalAssignment::build_Q(priority_queue< QElement > &Q,
                                              SetValuesFn &set_val_fn,
                                              double threshold,
                                              const vector<int> &qcol) const
  {
    for (auto c : qcol )
    {
      QElement qe;
      qe.c = c;
      set_val_fn.set_values(*this,qe);
      if (qe.valueA>threshold)
        Q.push( qe );
    }
  }
  
  
  bool IncrementalIntervalAssignment::tip_top( priority_queue<QElement> &Q, SetValuesFn &val_fn, double threshold, QElement &t)
  {
    result->debug_message("Q size %d ", (int) Q.size());
    if (Q.empty())
    {
      result->debug_message("Q empty, all done!!!\n");
      t = QElement();
      return true;
    }
    
    t = Q.top(); // copy
    Q.pop();
    
    // Check that the element value is up to date.
    // If not, then put it back on the queue and get the new top
    while ( val_fn.update_values(*this,t) )
    {
      result->debug_message("Q stale.\n");
      if (t.valueA>threshold)
      {
        //      // alternatively
        //      //   if t is not smaller than the current top, so return it as if it were the top
        //      if (! (t<Q.top()) )
        //      {
        //        break;
        //      }
        result->debug_message("Q reque.\n");
        Q.push(t);
      }
      else
      {
        result->debug_message("Q tossed.\n");
      }
      
      if (Q.empty())
      {
        result->debug_message("Q empty, all done!!!\n");
        t = QElement();
        return true;
      }
      t = Q.top();
      Q.pop();
    }
    
    // clear any remaining values below the threshold
    if (t.valueA<=threshold)
    {
      result->debug_message("Q being emptied.\n");
      Q = priority_queue<QElement>();
    }
    result->debug_message("Q tip top %d values A:%g B:%g C:, remaining Q size %d\n", t.c, t.valueA, t.valueB, t.valueC, (int) Q.size());
    
    return false;
  }

  // desire to improve variables that are far from their goal, towards their goal
  // positive means variable would get better
  // negative means variable would get worse
  // 0 means no goal "don't care"
  //   ersatz g:x ratio that will also allow zero and negative x's
  double increment_desire_goal(int x, int dx, double g)
  {
    if (g==0.)  // no goal
      return 0.;

    const auto dist_gx   = (g- x);      const auto abs_gx   = fabs(dist_gx);
    const auto dist_gxdx = (g-(x+dx));  const auto abs_gxdx = fabs(dist_gxdx);
    const auto denom = 0.5 + fabs(g);
                        
    // incrementing would make it farther away from its goal
    if (abs_gxdx>abs_gx)
      return -abs_gxdx / denom;
    // incrementing would make it not farther, but flip the sign of g-x
    const bool pos_gxdx( dist_gxdx>0. ), pos_gx( dist_gx>0. );
    if ( pos_gxdx != pos_gx )
    {
      // return a positive value, but of the smaller one
      return min(abs_gxdx, abs_gx) / denom;
    }
    // incrementing would make it not farther, and keep it on the same side
    return abs_gx / denom;
  }

  bool IncrementalIntervalAssignment::unbounded_bounds_improvement(int tc, int dx, map<int,int> &improved_cols,
                                                                   QWithReplacement &Q, SetValuesFn &val_bounds_Q, const double threshold,
                                                                   const MatrixSparseInt &MN,
                                                                   BlockingCols &blocking_cols, vector<int> &deferred_rows)
  {
    if (tc >= (int) MN.col_rows.size())
      return false;
    
    // use rows that don't make *any* variable closer to its bounds
    //   unique for bounds, which are often one-sided, in contrast to "improve" where we must make tradeoffs
    //   this doesn't seem to help much, but it doesn't hurt
    
    // pick *best* unbounded row, not just the first we find
    //   no need for a queue, since we just use the best one then quit
    //
    // if row's variable is below it's bound, increment by that row until the bound is met
    //   pick the best such row, the variable farthest from its lower bound
    // quality: if a variable is below its goal, increment by that row until <= goal.
    //   pick the best such row, the one farthest below its goal
    // quality: pick the variable least over its goal, increment it once, then recheck.
  
    auto &tc_rows=MN.col_rows[tc];
    int biggest_bounds_deficit = 0;
    double biggest_quality(0.), best_degrade(numeric_limits<double>::max()); // prefer smaller degrade
    int best_unbounded_row = -1;
    for (int m : tc_rows )
    {
      auto &row = MN.rows[m];
    
      // a row is bounded if when we add it forever to the solution some variable bound is violated
      // get sign of coefficient for tc in this row.
      const auto dxx = row.get_val(tc) > 0 ? dx : -dx;
      bool unbounded = true;
      int row_deficit = 0;
      double row_quality(0.), row_degrade(0);
      for (size_t i = 0; i < row.cols.size(); ++i)
      {
        auto c = row.cols[i];
        auto v = row.vals[i];
        auto dc = dxx * v; // how much x[c] will change, signed
        if (( (dc < 0) && (col_lower_bounds[c] != numeric_limits<int>::lowest()) ) ||
            ( (dc > 0) && (col_upper_bounds[c] != numeric_limits<int>::max()   ) ))
        {
          unbounded = false;
          deferred_rows.push_back(m);
          break;
        }
        const int x = col_solution[c];

        // bounds
        // we already checked that it was unbounded in the right way, so any deficit is one that would be improved by an increment
        const int bounds_deficit = abs(compute_out_of_bounds(c,x)); // maybe we should divide by v, too.
        row_deficit = max(row_deficit,bounds_deficit);
        
        // goals
        auto dg = increment_desire_goal(x, dc, goals[c]);
        if (dg>0)
        {
          row_quality = max(row_quality, dg);
        }
        else
        {
          row_degrade = max(row_degrade, fabs(dg));
        }
      }
      // criteria:
      // unbounded
      // bigger deficit
      //   tied deficit and smaller degrade
      //   tied deficit and tied degrade and larger quality
      if (unbounded &&
          (biggest_bounds_deficit<row_deficit ||
           (biggest_bounds_deficit==row_deficit && row_degrade<best_degrade) ||
           (biggest_bounds_deficit==row_deficit && row_degrade==best_degrade && biggest_quality<=row_quality) ) )
      {
        best_unbounded_row = m;
        biggest_bounds_deficit = row_deficit;
        biggest_quality = row_quality;
        best_degrade = row_degrade;
      }
    }
        
    // if there was an unbounded row, use the best one
    if (best_unbounded_row>=0)
    {
      const int m = best_unbounded_row;

      // debug
      if (result->log_debug)
      {
        // result->error_message("Found a case where unbounded works!\n");
        result->debug_message("row %d unbounded ",(int)m);
        MN.rows[m].print_row(result);
      }
      // debug
      
      // accept_strict_improvement is overkill and inefficient in this case, but it works
      bool self_block = false;
      int self_block_coeff = 0;
      int num_block;
      if (accept_strict_improvement(Q, val_bounds_Q, threshold, nullptr, 0., tc, MN.rows[m], dx,
                                    self_block, self_block_coeff, blocking_cols, num_block, improved_cols))
      {
        return true;
      }
      else
      {
        result->error_message("Incremental Interval Assignment algorithm error. The row was unbounded but was somehow blocked, a contradiction.\n");
      }
    }
    return false;
  }
  
  bool bigger_quality(pair<int,double>& a, 
                      pair<int,double>& b)
  {
    return a.second > b.second;
  }

  bool IncrementalIntervalAssignment::queue_rows_by_quality( int tc, int dx, const MatrixSparseInt &MN, int m,
                                                            vector< pair<int,double> > &rows_by_quality)
  {
    auto &row = MN.rows[m];
    const auto dxx = row.get_val(tc) > 0 ? dx : -dx;
    double qbest = 0., qworst = 0.;
    for (size_t i = 0; i < row.cols.size(); ++i)
    {
      const auto c = row.cols[i];
      const auto v = row.vals[i];
      const int x = col_solution[c]; // current value
      const auto dc = dxx * v; // how much x[c] will change, signed
      const int xdc = x+dc; // new value
      
      // bounds
      // will incrementing make the bounds worse?
      if ( abs(compute_out_of_bounds(c,xdc)) > abs(compute_out_of_bounds(c,x)) )
      {
        return false;
      }
      
      // goals
      auto dg = increment_desire_goal(x, dc, goals[c]);
      qbest  = max(qbest,dg);
      qworst = min(qworst,dg);
    }
    const auto q = qbest + qworst;
    rows_by_quality.emplace_back( m, q );
    return true;
  }

  // look for any row that works
  bool IncrementalIntervalAssignment::any_bounds_improvement(int tc, int dx, map<int,int> &improved_cols,
                                                             QWithReplacement &Q, SetValuesFn &val_bounds_Q, const double threshold,
                                                             const MatrixSparseInt &MN, BlockingCols &blocking_cols, vector<int> &deferred_rows,
                                                             int &smallest_self_block_coeff, bool &give_up)
  {
    // pick *best* row, not just the first we find that works
    vector< pair<int,double> > rows_by_quality;
    for (auto m : deferred_rows )
    {
      queue_rows_by_quality(tc, dx, MN, m, rows_by_quality);
    }
    if (use_best_improvement)
    {
      sort( rows_by_quality.begin(), rows_by_quality.end(), bigger_quality);
    }

    for (auto rq : rows_by_quality)
    {
      int m = rq.first;
      // double q = rq.second
      bool self_block = false;
      int self_block_coeff = 0;
      int num_block = 0;
      // satisfy bounds
      if (accept_strict_improvement(Q, val_bounds_Q, threshold, nullptr, 0., tc, MN.rows[m], dx,
                                    self_block, self_block_coeff, blocking_cols, num_block, improved_cols))
      {
        return true;
      }
      else if (self_block)
      {
        result->debug_message("row %d selfblocks col %d\n",(int)m,tc);
        // give up later on this variable if it's not possible to get a combination with a smaller coefficient
        smallest_self_block_coeff=min(smallest_self_block_coeff,self_block_coeff);
        
        // give up immediately if the increment is the smallest it could be
        if (smallest_self_block_coeff==1)
        {
          give_up=true;
          return false;
        }
      }
      else
      {
        result->debug_message("row %d blocked\n",(int)m);
      }
    }
    return false;
  }
  
  
  bool IncrementalIntervalAssignment::satisfy_bounds(const MatrixSparseInt&M, MatrixSparseInt&MN)
  {
    CpuTimer *timer=nullptr;
    if (result->log_debug_time)
    {
      timer = new CpuTimer;
    }
    
    bool success(true); // return true if the bounds are satisfied
    
    vector<int>cols_fail;
    if (!bounds_satisfied(cols_fail))
    {
      if (MN.col_rows.empty())
      {
        result->debug_message("\nIIA: unable to find a solution that satisfies the variable bounds because the nullspace is empty.\n");
        return false;
      }

      result->debug_message("\nIIA: finding a solution that satisfies the variable bounds\n");
      bool rc(true); // set to false if we give up on any variable
                     // to do
                     // find sensitivity of dependent vars to independent vars
      
      // use priority queue, variable farthest from its bound
      //   find a null-space-basis vector, or a linear sum of them, to move vars closer to their bound.
      //   Strict improvement required at each step.
      
      // priority queue for variable bounds
      // build Q
      SetValuesBounds val_bounds_Q;
      const double threshold=0.;
      QWithReplacement Q(this,val_bounds_Q,threshold);
      Q.build( cols_fail );
      
      // store the blocking column that had the fewest blocking variables
      BlockingCols blocking_cols, all_blocking_cols;
      vector<int> deferred_rows;
      map<int,int> improved_cols;
      
      // We use MB2 for doing the Gaussian elimination, then concatenating MN with any rows we found
      //   continuing the partial elimination of MB2 on subsequent interations is a big speedup, 3x
      // implicit MB2 order of columns, can be deduced by rows [0,r)
      MatrixSparseInt MB2(M);
      int r=0; // next row of MB2 to reduce
      
      // debug
      //    int max_iter=100*Q.size(); // debug quit flag, just comment out otherwise
      int iter=0;
      // debug
      
      QElement t;
      while (!Q.empty())
      {
        // debug
        result->debug_message("iter %d",iter++);
        //      if (0 && ++iter>max_iter)
        //      {
        //        result->error_message("IIA: max_iter %d reached when trying to satisfy variable bounds.\n",max_iter);
        //        break;
        //      }
        // debug
        
        if (Q.tip_top(t))
          break;
        
        // get sign of needed change
        int dx = compute_out_of_bounds(t.c,t.solution) > 0 ? 1 : -1;
        
        // shortcut
        // if incrementing t.c doesn't improve quality, then skip all the rest
        QElement u = t;
        val_bounds_Q.set_values(*this,u,col_solution[t.c]+dx);
        if (u<t)
        {
          
          // search the nullspace for a vector that makes t better,
          // check that it doesn't degrade the min-max
          // if none, try a combination of nullspace vectors...
          int smallest_self_block_coeff=numeric_limits<int>::max();
          bool fixed=false;
          bool give_up = (MN.col_rows.size() <= t.c) || MN.col_rows[t.c].empty();
          if (give_up)
          {
            result->debug_message("IIA nullspace: It's not possible to improve x%d to be within bounds because it isn't in the nullspace at all. Giving up on this variable.\n", t.c);
            rc=false;
          }
          else
          {
            blocking_cols.clear();
            deferred_rows.clear();

            const bool debug_pave_research_code(false);
            
            // look for unbounded MN
            fixed = false;
            if (!fixed)
            {
              fixed = unbounded_bounds_improvement(t.c,dx,improved_cols, Q,val_bounds_Q,threshold, MN, blocking_cols,deferred_rows);
              if (debug_pave_research_code && fixed) result->info_message("bound fixed: by unbounded row of M+Nullspace\n");
            }
            // look for any MN
            if (!fixed)
            {
              fixed = any_bounds_improvement(t.c,dx,improved_cols, Q,val_bounds_Q,threshold, MN, blocking_cols,deferred_rows,smallest_self_block_coeff, give_up);
              if (debug_pave_research_code && fixed) result->info_message("bound fixed: by any row of Nullspace\n");
            }
            if (debug_pave_research_code && !fixed)
            {
              if (give_up)
                result->info_message("bound not fixed yet, giving up.\n");
              else
                result->info_message("bound not fixed yet, doing Gaussian elimination.\n");
            }
  
            // have to consider a combination of rows
            if (!fixed)
            {
              result->debug_message("\nIIA nullspace: No single nullspace vector gave strict improvement towards the bounds of x%d\n", t.c);
              // Find a nullspace combination that doesn't have the blocking vars in the forbidden directions
              
              // if we fail, then don't put t back on the queue, but keep trying to improve the other vars
              
              if (give_up)
              {
                result->debug_message("IIA nullspace: It's not possible to improve x%d to be within bounds. Giving up on this variable.\n", t.c);
                rc=false;
              }
              else if (M.col_rows[t.c].empty())
              {
                result->debug_message("Variable x%d isn't in any M nullspace row, so its impossible to improve!\n",t.c);
                // this isn't an error if we are doing autoscheme and the surface doesn't admit that scheme.
                result->debug_message("Unable to assign intervals within bounds, check mesh scheme, vertex types and sweep directions.\n");
                rc=false;
              }
              else
              {
                // Find a nullspace combination that doesn't have the blocking vars in them (either plus or minus)
                auto old_num_rows=MN.rows.size();
                
                do
                {
                  int unblocked_row=-1;
                  const bool worked = MB2.gaussian_elimination(t.c, blocking_cols, smallest_self_block_coeff, improved_cols, all_blocking_cols, r, unblocked_row);
                  
                  if (worked)
                  {
                    const RowSparseInt &R=MB2.rows[unblocked_row];
                    
                    // debug
                    if (result->log_debug)
                    {
                      result->debug_message("Gaussian success, new row ");
                      R.print_row(result);
                      // M.print_matrix("new nullspace");
                    }
                    
                    
                    int num_block=0;
                    bool self_block=false;
                    int self_block_coeff=0;
                    const size_t old_size = blocking_size(blocking_cols);
                    
                    if (accept_strict_improvement(Q, val_bounds_Q, threshold, nullptr, 0., t.c, R, dx,
                                                  self_block, self_block_coeff, blocking_cols, num_block, improved_cols))
                    {
                      MN.push_row(R);
                      fixed=true;
                    }
                    else
                    {
                      // try again, unless there is no hope of progress
                      
                      // give up if we can't get the self_block coefficient smaller
                      if (self_block)
                      {
                        smallest_self_block_coeff=min(smallest_self_block_coeff,self_block_coeff);
                        if (smallest_self_block_coeff==1)
                        {
                          result->debug_message("  Giving up! self_blocking coefficient is 1. latest self_block_coeff=%d.\n", smallest_self_block_coeff);
                          give_up=true;
                        }
                      }
                      // give up if no progress
                      else
                      {
                        auto new_size = blocking_size(blocking_cols);
                        if (old_size==new_size)
                        {
                          result->debug_message("  Giving up! blocking_cols size=%d didn't increase.\n", (int)new_size);
                          give_up=true;
                        }
                      }
                    }
                  }
                  else // (!worked)
                  {
                    give_up=true;
                    if (result->log_debug)
                    {
                      result->debug_message("  Giving up! Gaussian elimination failed.\n");
                      {
                        result->info_message("MB2 r=%d reduced rows\n",r);
                        MB2.print_matrix("MB2");
                        result->info_message("all_blocking_cols ");
                        print_map(all_blocking_cols);
                      }
                    }
                  }
                  
                  if ( (size_t)r==MB2.rows.size())
                  {
                    result->debug_message("  Giving up! there are no more rows to reduce\n");
                    give_up=true;
                  }
                  
                  
                } while (!fixed && !give_up);
                if (!fixed)
                  rc=false;
                
                // debug
                if (result->log_debug)
                {
                  size_t num_new_rows = MN.rows.size() - old_num_rows;
                  result->debug_message("IIA: nullspace: %lu rows added; was %lu, now %lu.\n",
                                        (unsigned long) num_new_rows,
                                        (unsigned long) old_num_rows,
                                        (unsigned long) MN.rows.size());
                  if (num_new_rows>0)
                  {
                    result->info_message("bounds old last row + new rows\n");
                    const size_t i_start = (old_num_rows>0 ? old_num_rows-1 : 0);
                    
                    for (size_t i = i_start; i<MN.rows.size(); ++i)
                    {
                      result->debug_message("row %d ",(int)i);
                      MN.rows[i].print_row(result);
                    }
                  }
                  result->debug_message("Blocking x columns = ");
                  print_map(blocking_cols);
                  if (fixed)
                  {
                    result->debug_message("IIA nullspace: fixed x%d . We made a row with a non-zero coeff for it, but zero for the blocking cols.\n", t.c);
                  }
                  else
                  {
                    result->debug_message("\n\nIIA nullspace: FAILED to fix x%d. Blocking columns won. I'll keep trying to satisfy_bounds on the other variables.\n\n", t.c);
                    MN.print_matrix("MN nullspace");
                    MB2.print_matrix("MB2 Gaussian-elim working nullspace");
                  }
                }
                // debug
                
              }
            } // !fixed
          }
        } // if (u<t)
        else
        {
          result->debug_message("IIA bounds Q: x%d can't seem to be any better than %d. Giving up on it.\n", t.c, col_solution[t.c]);
        }
      } // while !Q.empty
      
      success = bounds_satisfied(cols_fail);
      
      // debug
      if (result->log_debug)
      {
        result->debug_message("IIA bounds Q empty.  ");
        if (cols_fail.empty())
        {
          result->debug_message("IIA bounds are satisfied, ");
          if (rc==false)
          {
            result->debug_message("despite giving up on some variables along the way.\n");
          }
          else
            result->debug_message("as expected.\n");
        }
        else
        {
          result->debug_message("IIA bounds are UNSATISFIED for %lu variables: ", (unsigned long) cols_fail.size());
          print_vec(result, cols_fail);
          result->debug_message("\n");
        }
      }
      
      // bounds should be satisfied as long as rc was false
      if (rc && !success)
        result->error_message("IIA bounds algorithm failure. The Q never got stuck, but the bounds were still not satisfied at the end.\n");
      
    }
    else
    {
      result->debug_message("\nIIA: initial solution already satisfies the variable bounds\n");
//      result->warning_message("\nIIA: initial solution already satisfies the variable bounds\n");
    }

    // debug
    if (result->log_debug_time)
    {
      double time_bounds = timer->cpu_secs();
      result->info_message("Satisfy variable bounds time %f\n",time_bounds);
      delete timer;
    }
    else if (result->log_debug)
    {
      result->debug_message("Satisfy variable bounds done\n");
    }
    
    return success;
  }
  
  // utility method for accept_strict_improvement
  //
  // tell blocking_cols that col c increment by v is a blocker
  // v is the min blocking coefficient encountered,
  //   except if we have both a positive and negative block, then we set it to 0 (could be an error if the values are high...)
  void add_blocker( BlockingCols &blocking_cols, int c, int v )
  {
    assert(c>=0);
    assert(v!=0);
    auto old = blocking_cols.find(c);
    if (old==blocking_cols.end())
    {
      if (v<0)
        blocking_cols.emplace(c, make_pair(v,block_max));
      else
        blocking_cols.emplace(c, make_pair(block_min,v));
    }
    else
    {
      if (v<0)
      {
        auto &lo = old->second.first;
        lo = max(lo, v);
      }
      else
      {
        auto &hi = old->second.second;
        hi = min(hi, v);
      }
    }
  }
  
  // for improve
  // notionally true if !(s<t)
  bool IncrementalIntervalAssignment::is_improvement(QElement &s,
                                                     QElement &t,
                                                     SetValuesFn *constraint_fn,
                                                     double constraint_threshold ) const
  {
    if (constraint_fn)
    {
      if (!(s<t))
      {
        result->debug_message("Increment would make x%d have unacceptable quality wrt its goal.\n", s.c);
        return false;
      }
      else
      {
        QElement qc;
        qc.c=s.c;
        constraint_fn->set_values(*this, qc);
        if (qc.valueA>constraint_threshold)
        {
          result->debug_message("Increment would make x%d not satisfy its variable bounds.\n", qc.c);
          return false;
        }
      }
    }
    else
    {
      if (!(s.valueA<t.valueA))
      {
        result->debug_message("Increment would make x%d have unacceptable quality wrt its bounds.\n", s.c);
        return false;
      }
    }
    return true;
  }
  
  bool IncrementalIntervalAssignment::accept_strict_improvement(QWithReplacement &Q,
                                                                SetValuesFn &quality_fn,
                                                                double quality_threshold,
                                                                SetValuesFn *constraint_fn,
                                                                double constraint_threshold,
                                                                int tc,
                                                                const RowSparseInt &mrow,
                                                                int dx,
                                                                bool &self_block,
                                                                int &self_block_coeff,
                                                                BlockingCols &blocking_cols,
                                                                int &num_block,
                                                                map<int,int> &improved_cols)
  {
    num_block=0;
    const size_t n = mrow.cols.size(); // number of variables changed by mrow
    
    // flip sign of dx so dx*row moves t.c in the right direction
    if (mrow.get_val(tc)<0)
      dx=-dx;
    
    // compare t to all the QElements changed by solution+=dx*mrow
    vector<QElement> qold(n), qnew(n);
    int ti = -1; // index of t in the above
    
    // update solution, get new element quality values
    for (size_t i = 0; i < n; ++i)
    {
      auto c = mrow.cols[i];
      
      if (c == tc)
      {
        ti=(int)i;
      }
      
      // old q value
      qold[i].c = c;
      quality_fn.set_values(*this,qold[i]); // sets qold[i].solution // optimization, just get this from QWithReplacement !!!
      
      // update solution
      auto v = mrow.vals[i]; // coefficient
      col_solution[c] += dx*v;
      
      // new q value
      qnew[i].c = c;
      quality_fn.set_values(*this,qnew[i]); // sets qnew[i].solution
      
    }
    assert(ti>-1);
    
    // check new quality and constraints
    // check self-block first, because if it self-blocks, we don't want to accumulate any other blockers
    {
      if (!is_improvement(qnew[ti],qold[ti],constraint_fn,constraint_threshold))
      {
        self_block=true;
        auto c = mrow.cols[ti]; // column
        auto v = mrow.vals[ti]; // coefficient
        self_block_coeff=abs(v);
        
        result->debug_message("x%d self-blocks with coeff %d.\n", c,v);
      }
    }
    
    // c
    if (!self_block)
    {
      for (size_t i = 0; i < n; ++i)
      {
        if ((int)i==ti)
          continue;
        
        if (!is_improvement(qnew[i],qold[ti],constraint_fn,constraint_threshold))
        {
          ++num_block;
          auto c = mrow.cols[i]; // column
          auto v = dx * mrow.vals[i]; // coefficient, with corrected sign and magnitude
          add_blocker(blocking_cols,c,v);
          
          // debug
          // auto qn = qnew[i];
          // auto qo = qold[ti];
          result->debug_message("x%d increment by %d blocks desired increment of x%d; adding it to blockers.\n", c, v, tc);
          // debug
        }
      }
    }
    
    // blocked? restore solution
    if (num_block || self_block)
    {
      for (size_t i = 0; i < n; ++i)
      {
        auto c = mrow.cols[i];
        col_solution[c] = qold[i].solution;
      }
      
      // debug
      if (result->log_debug)
      {
        result->debug_message("Rejecting increment of x%d=%d by %d * ", tc, col_solution[tc], dx);
        mrow.print_row(result);
      }
      // debug
      
      return false;
    }
    
    // OK, we can increment by at least 1, try incrementing more!
    // keep improving as long as everyone is better than the prior value for t. In particular, stop if t stops improving.
    //   we could do a binary search, but slowly stepping by dx until it doesn't get better is probably fast enough, and maybe even better.
    // to do: experiment with what increment helps best
    auto total_dx=dx;
    bool good=true;
    while (good)
    {
      QElement worst = qnew[ti];
      total_dx+=dx;
      for (size_t i = 0; i < n; ++i)
      {
        auto c = mrow.cols[i];
        auto v = mrow.vals[i];
        col_solution[c] += dx*v;
        quality_fn.set_values(*this,qnew[i]);
        
        // check quality
        if (!is_improvement(qnew[i],worst,constraint_fn,constraint_threshold))
        {
          // undo partial increment
          for (size_t ii=0; ii<=i; ++ii)
          {
            auto cc = mrow.cols[ii];
            auto vv = mrow.vals[ii]; // coefficient
            col_solution[cc] -= dx*vv;
            quality_fn.set_values(*this,qnew[ii]);
          }
          total_dx-=dx;
          good=false;
          break; // out of loop i=0:n
        }
      }
      // if good, we will try incrementing again
    }
    
    // it improves min max t, accept it!
    if (result->log_debug)
    {
      result->debug_message("Accepted increment of x%d from %d to %d via %d * ",
                            tc, col_solution[tc] - total_dx * mrow.vals[ti], col_solution[tc], total_dx);
      mrow.print_row(result);
      print_solution("after accepted increment");
      // print_solution_summary("after accepted increment");
    }
    
    // update improved_cols
    for (size_t i = 0; i < n; ++i)
    {
      auto &qo=qold[i];
      auto &qn=qnew[i];
      const auto d = qn.solution-qo.solution;
      assert(d);
      improved_cols[qn.c]+=d;
    }
    
    // update Q
    Q.update(mrow.cols);
    
    return true;
  }
  
  void IncrementalIntervalAssignment::improve_solution(const MatrixSparseInt&M, MatrixSparseInt&MN)
  {
    if (MN.col_rows.empty())
    {
      result->debug_message("IIA: unable to improve_solution because nullspace is empty\n.");
      return;
    }
    
    // debug data
    CpuTimer *timer(nullptr), *timer_0(nullptr), *timer_1(nullptr), *accept_timer(nullptr);
        
    if (result->log_debug_time)
    {
      timer = new CpuTimer;
      timer_0 = new CpuTimer;
      timer_1 = new CpuTimer;
      accept_timer = new CpuTimer;
    }
    double time_wasted=0., time_improve=0.;
    double time_first_wasted=-1.;
    int num_iter=0;
    int num_useful_after_wasted=0;
    int num_useful=0;
    double time_0=0., time_0f=0., time_0not=0.,time_1=0.;
    size_t accept_calls(0);
    double accept_time=0.;
    size_t MN_rows_init=MN.rows.size();
    // debug
    
    // bounds are already satisfied
    // any step must stay within bounds
    // we improve towards the (floating point) goal intervals
    
    result->debug_message("\nIIA: improving feasible solution towards goals.\n");
    // use priority queue, variable farthest from its goal
    //   find a null-space-basis vector, or a linear sum of them, to move vars closer to their goals
    //   Strict improvement required at each step.
    
    int rc=true;
    
    // COMMIT test suite fails/passes invariant to which of thsese two priority fn we use
    // The first one worked better for the mesh scaling IIA
    SetValuesRatioR priority_fn_RatioR;  // based on current value+/-1
    const double priority_threshold_RatioR=0.01;
    SetValuesRatioG priority_fn_RatioG;  // based on current value, except if optimal (close to goal G) don't select it
    const double priority_threshold_RatioG=1.;

    SetValuesRatio quality_fn;
    const double quality_threshold=0;
    
    SetValuesBounds constraint_fn;
    const double constraint_threshold=0;
    
    // Instead of a priority_queue, we use QWithReplacement
    // That way, when we get a new priority, we can just replace it,
    // More runtime to maintain than the priority_queue, but time is dominated by actually processing and trying to improve the Q elements in the expensive cases.
    // With a priority_queue, we have duplicates in the Q and spend too much time once we are stuck.
    
    // build Q
    SetValuesFn *priority_fn =  use_Ratio_for_improvement_queue ? (SetValuesFn*) &priority_fn_RatioG : (SetValuesFn*) &priority_fn_RatioR;
    double priority_threshold = use_Ratio_for_improvement_queue ?          priority_threshold_RatioG :          priority_threshold_RatioR;
    QWithReplacement Q(this,*priority_fn,priority_threshold);
    // get all the columns, except the tied_variables
    vector<int>untied_int_vars;
    untied_int_vars.reserve(used_col+1);
    for (int c = 0; (size_t) c < tied_variables.size(); ++c)
    {
      if (tied_variables[c].size()!=1)
        untied_int_vars.push_back(c);
    }
    Q.build( untied_int_vars );
    // Q.print();
    
    // blocking_cols are just the ones that we discovered that block the current variable from incrementing
    // all_blocking_cols are the union blocking_cols from all prior increment attempts
    // improved_cols are the variables that were incremented (since they were added to all_blocking_cols, i.e. since the The  last gaussian elimination call)
    BlockingCols blocking_cols, all_blocking_cols;
    map<int,int> improved_cols;
    // implicit MV2 order of columns, can be deduced by rows [0,r)
    MatrixSparseInt MV2(M);
    int r=0; // next row of MV2 to reduce
    
    // debug
    int max_iter=10; // 100*Q.size(); // debug quit flag, just comment out otherwise
    int iter=0;
    // debug
    
    QElement t;
    while (!Q.empty())
    {
      
      // debug
      if (result->log_debug_time)
      {
        timer->cpu_secs(); // iter start
      }
      if (result->log_debug)
      {
        ++iter;
        result->debug_message("iter %d ",iter);
        if (/* DISABLES CODE */ (0) && iter>max_iter)
        {
          result->error_message("IIA: max_iter %d reached when trying to improve solution towards goals.\n",max_iter);
          break;
        }
      }
      // debug
      
      if (Q.tip_top(t))
        break;
      
      int dx = t.dx;
      assert(dx==1 || dx==-1);
      
      // debug
      if (result->log_debug)
      {
        result->debug_message("\nIIA improve_solution: x%d priority %g, sol %d, attempting %+d towards ", t.c, t.valueA, t.solution, dx);
        if (tied_variables[t.c].size()>1)
        {
          double glo, ghi;
          get_tied_goals(t.c, glo, ghi);
          result->debug_message("tied goals [%g, %g]\n", glo, ghi);
        }
        else
        {
          result->debug_message("goal %g\n", goals[t.c]);
        }
      }
      // debug
      
      if (t.valueA>priority_threshold)
      {
        // debug
        if (result->log_debug_time)
        {
          timer_0->cpu_secs();
        }
        
        bool fixed=false;
        bool give_up=false;
        int smallest_self_block_coeff=numeric_limits<int>::max();
        // does an existing row work?
        
        // see if an existing row of MN works?
        {
          blocking_cols.clear();
          int num_block=0;
          auto &tc_rows=MN.col_rows[t.c];
          if(!tc_rows.empty())
          {
            // pick *best* row, not just the first we find that works
            vector< pair<int,double> > rows_by_quality;
            for (auto m : tc_rows )
            {
              queue_rows_by_quality(t.c, dx, MN, m, rows_by_quality);
            }
            if (use_best_improvement)
            {
              sort( rows_by_quality.begin(), rows_by_quality.end(), bigger_quality);
            }

            for (auto rq : rows_by_quality)
            {
              int m = rq.first;
              // double q = rq.second

              // debug
              if (result->log_debug_time)
              {
                accept_timer->cpu_secs();
              }
              if (result->log_debug)
              {
                ++accept_calls;
              }

              //improve_solution
              bool self_block = false;
              int self_block_coeff = 0;
              const bool worked = accept_strict_improvement(Q,
                                                            quality_fn, quality_threshold,
                                                            &constraint_fn, constraint_threshold,
                                                            t.c, MN.rows[m], dx,
                                                            self_block, self_block_coeff, blocking_cols, num_block, improved_cols);
              
              // debug time
              if (result->log_debug_time)
              {
                accept_time+=accept_timer->cpu_secs();
              }
              
              if (worked)
              {
                fixed=true;
                break;
              }
              else if (self_block)
              {
                // give up later on this variable if it's not possible to get a combination with a smaller coefficient
                smallest_self_block_coeff=min(smallest_self_block_coeff,self_block_coeff);
                
                // give up immediately if the increment is the smallest it could be
                if (smallest_self_block_coeff==1)
                {
                  give_up=true;
                  break;
                }
              }
            } // for
          }
        }
        
        // debug
        if (result->log_debug_time)
        {
          auto t0 = timer_0->cpu_secs();
          time_0+=t0;
          if (fixed)
            time_0f+=t0;
          else
            time_0not+=t0;
          timer_1->cpu_secs();
        }
        // debug
        
        if (!fixed)
        {
          // if  (size_t)r+1<MV2.rows.size() then we can't reduce any more, but we can use the reduced MV2 to try to find an improvement vector
          result->debug_message("\nIIA nullspace: No single nullspace vector gave strict improvement towards the goal of x%d\n", t.c);
          if (give_up)
          {
            result->debug_message("IIA nullspace: It's not possible to improve x%d. Giving up on this variable.\n", t.c);
            rc=false;
          }
          else if (M.col_rows[t.c].empty()) // should we be checking MV2 instead of M?
          {
            result->debug_message("Variable x%d isn't in any M nullspace row, so its impossible to improve!\n",t.c);
            rc=false;
          }
          else
          {
            auto old_num_rows=MN.rows.size();

            // Add Nullspace to flip sign of blockers
            // Find a nullspace combination that has the opposite sign for the blocking vars
            if (use_Nullspace_sumeven_to_flip_sign_of_blocker && !Nullspace.rows.empty())
            {
              
              if (debug_Nullspace_sumeven_to_flip_sign_of_blocker)
              {
                result->info_message("\n\nTrying to %screment x%d\n", (t.dx>0 ? "in" : "de"), t.c);
                result->info_message("By using a sum-even nullspace row to flip the sign of blocking coeff, for variables in blocking_cols\n");
              }


              int num_block=0;
              // Should we just step over the rows of M instead of MN?
              auto &tc_rows=MN.col_rows[t.c];
              // figure out how to do the best one later
              for (auto r : tc_rows)
              {
                const auto r_val = MN.rows[r].get_val(t.c);
                const auto r_sign = (r_val > 0 ? 1 : -1) * t.dx; // sign after row is negated if needed to get sign == dx.
                
                if (smallest_self_block_coeff <= abs(r_val))
                {
                  // this row has a self-blocking coefficient, try some other row
                  continue;
                }
                // check sign of tc_rows.coeff t.c
                // find blockers
                // typedef map<int, pair<int,int> > BlockingCols;
                // copy row
                RowSparseInt row(MN.rows[r]); // copy
                const bool negate_row_r = ((r_val>0) != t.dx>0);
                if (negate_row_r)
                {
                  row.multiply(-1);
                }
                // row has the right sign for t.c so that we would add the row to improve solution[t.c]
                
                if (debug_Nullspace_sumeven_to_flip_sign_of_blocker)
                {
                  result->info_message("Initial MN row %d, with %s sign, ", r, (negate_row_r ? "flipped" : "original"));
                  row.print_row(result);
                  result->info_message("Current blocking_cols ");
                  print_map(blocking_cols);
                }

                bool try_row = true;
                bool progress = false;
                for (auto &blocking_col : blocking_cols)
                {
                  
                  // 1. see if row has blocking col with coefficient +/-1, compare to sign of blocking_cols i.e. first and second
                  // 2. if so search Nullspace for a vector with a "2" coeff for that blocker
                  // to do: figure out what to do if it changes the coeff of other blockers
                  //   add it to copy
                  //   see if it works
                  //   if not, does it also block with coeff +1? then we need to eliminate it.
                  //     otherwise, some other blocker has a bad coeff. That's the one we need to eliminate, not this one.

                  // 1.
                  int bc = blocking_col.first;
                  
                  auto bc_val = row.get_val(bc);
                  if (bc_val == 0)
                  {
                    continue; // blocker not in row, we're already good wrt this variable
                  }
                  if (abs(bc_val) % 2==0) // coefficient is even
                  {
                    // to do, indicate we want to add to make the coefficient zero
                    // or just deal with this during the Gaussian elimination phase
                    try_row = false;
                    continue;
                  }
                  // to do, generalize to look for any combo of 1 and 2 coefficients
                  // for now, just work on the coeff == 1 case
                  // e.g. we could also try to fix bc_val = -3, etc
                  if ( abs(bc_val)!=1 )
                  {
                    if (debug_Nullspace_sumeven_to_flip_sign_of_blocker)
                    {
                      result->info_message("Only bc_val == 1 or -1 is implemented currently. bc_val %d unsupported, skipping\n",bc_val);
                    }
                    try_row = false;
                    continue;
                  }

                  int bc_r_sign = bc_val > 0 ? 1 : -1;
                  auto block_neg = blocking_col.second.first;
                  auto block_pos = blocking_col.second.second;
                  // blocked in both directions? give up on r
                  if (block_neg == -1 && block_pos == 1)
                  {
                    // row has a blocker that we must eliminate for any chance of success
                    try_row = false;
                    // but keep trying to flip the sign of other blockers, so elimination will work on row
                    continue;
                  }

                  // do we want to add (true) or subtract (false) 2*nullspacerow to row?
                  try_row = true;
                  bool want_add_two = true;
                  if (block_neg == bc_r_sign)
                  {
                    want_add_two = true;
                  }
                  else if (block_pos == bc_r_sign)
                  {
                    want_add_two = false; // want to subtract 2
                  }
                  else
                  {
                    // block_pos or neg magnitude is greater than 1
                    // OK to try row without doing anything to this bc in row
                    continue;
                  }
                  
                  // 2.
                  // check rows in nullspace containing variable bc
                  assert(try_row);
                  try_row = false;
                  const auto &bc_rows = Nullspace.col_rows[bc];
                  for (int rN : bc_rows)
                  {

                    auto &Nrow = Nullspace.rows[rN];
                    int rN_bc_val = Nrow.get_val(bc);
                    // want a 2 for the blocking col, and a zero for t.c
                    if ( abs(rN_bc_val) != 2 || Nrow.get_val(t.c)!=0)
                    {
                      continue; // try the next row
                    }
                    bool rN_bc_pos = (rN_bc_val>0);
                    int mult = (want_add_two == rN_bc_pos) ? 1 : -1;

                    // debug
                    if (debug_Nullspace_sumeven_to_flip_sign_of_blocker)
                    {
                      result->info_message("Flipping sign of blocking col %dx%d in ", bc_val, bc);
                      row.print_row(result);
                      result->info_message(" by adding %d times ", mult);
                      Nrow.print_row(result);
                    }
                    
                    Nullspace.row_r_add_ms(row, mult, Nrow);
                    // to do
                    // Research
                    //   adding Nrow can make row worse, by adding more blocking cols, or making their coefficients worse
                    //   need a strategy for evaluating if it is worse and going back to original row if things are bad,
                    //   try various combinations of rows, etc.
                    //   could be complicated
                    
                    
                    // debug
                    if (debug_Nullspace_sumeven_to_flip_sign_of_blocker)
                    {
                      result->info_message(" to get ");
                      row.print_row(result);
                      result->info_message("\n");
                    }

                    // did it work?
                    const auto new_bc_val = row.get_val(bc);
                    if (new_bc_val>0)
                    {
                      assert(new_bc_val < block_pos);
                    }
                    else
                    {
                      assert(new_bc_val > block_neg);
                    }
                    assert( row.get_val(bc) == -bc_r_sign );
                    try_row = true;
                    break;
                  }
                  if (!try_row)
                  {
                    break;
                  }
                } // for blocking_col

                if (try_row)
                {
                  if (debug_Nullspace_sumeven_to_flip_sign_of_blocker)
                  {
                    result->info_message("Trying ");
                    row.print_row(result);
                  }

                  // try to improve solution using row
                  bool self_block = false;
                  int self_block_coeff = 0;
                  const bool worked = accept_strict_improvement(Q,
                                                                quality_fn, quality_threshold,
                                                                &constraint_fn, constraint_threshold,
                                                                t.c, row, dx,
                                                                self_block, self_block_coeff, blocking_cols, num_block, improved_cols);
                  if (worked)
                  {
                    if (debug_Nullspace_sumeven_to_flip_sign_of_blocker)
                    {
                      result->info_message("debug_Nullspace_sumeven_to_flip_sign_of_blocker worked :-O\n\n");
                    }
                    fixed = true;
                    break;
                  }
                  else
                  {
                    if (debug_Nullspace_sumeven_to_flip_sign_of_blocker)
                    {
                      result->info_message("failed, I'll try another row if any :-(\n\n");
                    }
                  }
                }
                
              } // for r
              if (debug_Nullspace_sumeven_to_flip_sign_of_blocker)
              {
                result->info_message("Moving on.\n\n");
              }

            }
            
            
            
            // Gaussian elim approach
            if (!fixed)
            {
            // Find a nullspace combination that doesn't have the blocking vars in them (either plus or minus)
            // try up to M.rows.size() times, since that's as many variables as we can eliminate
            // int max_tries = M.rows.size();
            do
            {
              int unblocked_row=-1;
              const bool worked = MV2.gaussian_elimination(t.c, blocking_cols, smallest_self_block_coeff, improved_cols, all_blocking_cols, r, unblocked_row);
              
              if (worked)
              {
                const RowSparseInt &R = MV2.rows[unblocked_row];
                
                // debug
                if (result->log_debug)
                {
                  result->debug_message("Gaussian success, new row ");
                  R.print_row(result);
                  // M.print_matrix("new nullspace");
                }
                
                
                int num_block=0;
                bool self_block=false;
                int self_block_coeff=0;
                const size_t old_size = blocking_size(blocking_cols);
                if (accept_strict_improvement(Q,
                                              quality_fn, quality_threshold,
                                              &constraint_fn, constraint_threshold,
                                              t.c, R, dx,
                                              self_block, self_block_coeff, blocking_cols, num_block, improved_cols))
                {
                  MN.push_row(R);
                  fixed=true;
                }
                else
                {
                  // try again, unless there is no hope of progress
                  
                  // give up if we can't get the self_block coefficient smaller
                  if (self_block)
                  {
                    smallest_self_block_coeff=min(smallest_self_block_coeff,self_block_coeff);
                    if (smallest_self_block_coeff==1)
                    {
                      result->debug_message("  Giving up! self_blocking coefficient is 1. latest self_block_coeff=%d.\n", smallest_self_block_coeff);
                      give_up=true;
                    }
                  }
                  // give up if no progress
                  else
                  {
                    auto new_size = blocking_size(blocking_cols);
                    if (old_size==new_size)
                    {
                      result->debug_message("  Giving up! blocking_cols size=%d didn't increase.\n", (int)new_size);
                      give_up=true;
                      
                      // save the row for fine_tune, and future iterations
                      MN.push_row(R); // zzyk is this actually helpful?
                    }
                  }
                }
              }
              else // (!worked)
              {
                give_up=true;
                if (result->log_debug)
                {
                  result->debug_message("  Giving up! Gaussian elimination failed.\n");
                  {
                    result->info_message("MV2 r=%d reduced rows\n",r);
                    MV2.print_matrix("MV2");
                    result->info_message("all_blocking_cols ");
                    print_map(all_blocking_cols);
                  }
                }
              }
              
              if ( (size_t)r==MV2.rows.size())
              {
                result->debug_message("  Giving up! there are no more rows to reduce\n");
                give_up=true;
              }
            } while (!fixed && !give_up);
            }
            
            if (!fixed)
            {
              rc = false;
              
              // debug
              if (/* DISABLES CODE */ (0))
              {
                result->info_message("var x%d stuck at %d, not fixed ", t.c, t.solution);
                if (give_up)
                  result->info_message("giving up ");
                // result->info_message("Blocking self_block=%s num_block=%d cols ", (self_block ? "yes" : "no"), num_block);
                // print_map(blocking_cols);
              }
              // debug
            }
            
            // debug
            if (result->log_debug)
            {
              size_t num_new_rows = MN.rows.size() - old_num_rows;
              result->debug_message("IIA: nullspace: %lu rows added; was %lu, now %lu.\n",
                                    (unsigned long) num_new_rows,
                                    (unsigned long) old_num_rows,
                                    (unsigned long) M.rows.size());
              if (num_new_rows>0)
              {
                result->info_message("improve old last row + new rows\n");
                const size_t i_start = (old_num_rows>0 ? old_num_rows-1 : 0);
                
                for (size_t i = i_start; i<MN.rows.size(); ++i)
                {
                  result->debug_message("row %d ",(int)i);
                  MN.rows[i].print_row(result);
                }
              }
              result->debug_message("Blocking x columns = ");
              print_map(blocking_cols);
              if (fixed)
              {
                result->debug_message("IIA nullspace: fixed x%d . We made a row with a non-zero coeff for it, but zero for the blocking cols.\n", t.c);
              }
              else
              {
                result->debug_message("IIA nullspace: FAILED to fix x%d . Blocking columns won. I'll keep trying to improve the other variables.\n\n", t.c);
                // M.print_matrix("failed nullspaceD");
              }
            }
            // debug
            
          } // else (no existing row was unblocked)
        } // !fixed
        
        // debug
        if (result->log_debug_time)
        {
          time_1 += timer_1->cpu_secs();
          double iter_time = timer->cpu_secs(); // iter end
          if (fixed)
          {
            time_improve+=iter_time;
          }
          else
          {
            if (time_first_wasted==-1.)
            {
              time_first_wasted=time_improve+time_wasted;
            }
            time_wasted+=iter_time;
          }
        }
        if (result->log_debug)
        {
          ++num_iter;
          if (fixed)
          {
            ++num_useful;
            if (num_useful<num_iter)
            {
              ++num_useful_after_wasted;
            }
          }
        }
        // debug
        
      } // if (u<t)
      else
      {
        result->debug_message("IIA goals Q: x%d can't seem to be any better than %d; its goal is %g. Giving up on it.\n",
                              t.c, col_solution[t.c],goals[t.c]);
      }
    } // while !Q.empty
    
    
    // debug
    // time on wasted iterations?
    result->debug_message( "Return code rc %d\n",rc);
    if (result->log_debug && result->log_debug_time)
    {
      auto total_time = (time_wasted+time_improve);
      auto time_after_first_waste = total_time - time_first_wasted;
      if (time_first_wasted==-1.)
      {
        result->debug_message( "Improve solution, %d useful iteration time = %f (100 percent)\n",num_iter, total_time);
      }
      else
      {
        int num_wasted = num_iter-num_useful;
        result->debug_message("Improve solution %d wasted iterations = %f (%f percent), %d useful iterations = %f (%f percent). %d time  %f (%f percent) on useful iterations before first failure, %d useful iterations after. total time %f (%f percent) after.\n",
                              num_wasted, time_wasted, time_wasted*100. / total_time ,
                              num_useful,
                              time_improve, time_improve*100. / total_time,
                              num_useful-num_useful_after_wasted,
                              time_first_wasted, time_first_wasted*100. / total_time,
                              num_useful_after_wasted,
                              time_after_first_waste, time_after_first_waste*100. / total_time
                              );
      }
      result->debug_message("section 0 time %f (%f fixed, %f not), section 1 time %f\n",time_0, time_0f, time_0not, time_1);
      result->debug_message("section 0 accept_strict calls %lu, time %f\n",accept_calls,accept_time);
      result->debug_message("number of rows: start %lu, end %lu\n", MN_rows_init, MN.rows.size() );
      delete timer;
      delete timer_0;
      delete timer_1;
      delete accept_timer;
    }
    if (result->log_debug)
    {
      const int num_wasted = num_iter-num_useful;
      if (num_wasted==0)
      {
        result->debug_message( "Improve solution, %d useful iterations (100 percent)\n",num_iter);
      }
      else
      {
        result->debug_message("Improve solution %d wasted iterations, %d useful iterations. %d on useful iterations before first failure, %d useful iterations after.\n",
                              num_wasted,
                              num_useful,
                              num_useful-num_useful_after_wasted,
                              num_useful_after_wasted
                              );
      }
      result->debug_message("section 0 accept_strict calls %lu\n",accept_calls);
      result->debug_message("number of rows: start %lu, end %lu\n", MN_rows_init, MN.rows.size() );
    }
    else if (result->log_debug_time)
    {
      auto total_time = (time_wasted+time_improve);
      auto time_after_first_waste = total_time - time_first_wasted;
      if (time_first_wasted==-1.)
      {
        result->debug_message( "Improve solution, time = %f (100 percent)\n", total_time);
      }
      else
      {
        result->debug_message("Improve solution wasted iterations = %f (%f percent), useful iterations = %f (%f percent). time  %f (%f percent) on useful iterations before first failure, useful iterations after. total time %f (%f percent) after.\n",
                              time_wasted, time_wasted*100. / total_time ,
                              time_improve, time_improve*100. / total_time,
                              time_first_wasted, time_first_wasted*100. / total_time,
                              time_after_first_waste, time_after_first_waste*100. / total_time
                              );
      }
      result->debug_message("section 0 time %f (%f fixed, %f not), section 1 time %f\n",time_0, time_0f, time_0not, time_1);
      result->debug_message("section 0 accept_strict calls, time %f\n", accept_time);
      delete timer;
      delete timer_0;
      delete timer_1;
      delete accept_timer;
    }
    // debug
    
    // debug
    //  if (0)
    //  {
    //    result->info_message("Nullspace after done adding rows to improve solution: \n");
    //    M.print_matrix();
    //    print_solution("improved solution");
    //  }
    // debug
    
    // debug
    if (result->log_debug)
    {
      result->debug_message("IIA goal Q empty.\n");
      print_solution("improved solution");
      // print_solution_summary("improved solution");
    }
    // debug
    
  }
  
  
  double IncrementalIntervalAssignment::get_R(int c, int &x, double &glo, double &ghi) const
  {
    get_tied_goals(c,glo,ghi);
    x = get_solution(c);

    if (glo==0)
    {
      return 0.;
    }
    else
    {
      if (x<=0)
      {
        const auto R1lo = glo<1 ? 1/glo : glo;
        const auto R1hi = ghi<1 ? 1/ghi : ghi;
        const auto R1 = std::max(R1lo,R1hi);
        return 1000.*(2-x)*R1;
      }
      if (x<=glo)
      {
        return ghi/x; // >= glo/x
      }
      else if (x>=ghi)
      {
        return x/glo; // >= x/glo
      }
      // glo < x < ghi
      else
      {
        return std::max( x/glo, ghi/x );
      }
    }
  }

  void IncrementalIntervalAssignment::compute_quality_ratio(vector<QElement> &q, const vector<int> &cols) const
  {
    SetValuesRatio fn;
    q.clear();
    q.reserve( cols.size()*2 );
    double glo, ghi;
    for (auto c : cols)
    {
      // tied, set to worse of ghi and glo
      if (get_tied_goals(c,glo,ghi))
      {
        // hi
        QElement qhi;
        qhi.c = c;
        qhi.solution = col_solution[c];
        qhi.solution_set = true;
        // lo
        QElement qlo(qhi);

        fn.set_values_goal(*this,ghi,qhi);
        fn.set_values_goal(*this,glo,qlo);
        
        // assign the larger value
        q.push_back( qlo<qhi ? qhi : qlo );
      }
      // not tied, lo==hi
      else
      {
        assert(glo==ghi);
        q.emplace_back();
        auto &qe = q.back();
        qe.c = c;
        qe.solution = col_solution[c];
        qe.solution_set = true;
        fn.set_values_goal(*this,glo,qe);
      }
    }
    sort(q.begin(),q.end());
  }
  
  bool IncrementalIntervalAssignment::increment_solution(const RowSparseInt &R, int dx)
  {
    bool rc=true;
    for (size_t i = 0; i < R.cols.size(); ++i)
    {
      auto c = R.cols[i];
      auto v = R.vals[i];
      col_solution[c] += dx*v;
      if (col_solution[c]<col_lower_bounds[c] || col_solution[c]>col_upper_bounds[c])
      {
        rc=false;
      }
    }
    return rc;
  }
  
  void IncrementalIntervalAssignment::fine_tune(const MatrixSparseInt&MN)
  {
    // see if any row of M will improve the lexicographic min-max of the deviations
    
    // while there was improvement
    //   iterate over rows of MN
    //     try updating the solution by +/- the row
    //       compare lex quality of the row columns before and after
    //       if better, even if not strictly better, then accept the increment
    
    bool improved = false;
    
    vector<QElement> qold, qnew;
    
    //  print_problem("before fine tune");
    //  M.print_matrix("M");
    //  N.print_matrix("N");
    
    do
    {
      improved = false;
      const auto maxr = MN.rows.size();
      for (size_t i = 0; i < maxr; ++i)
      {
        const RowSparseInt &R = MN.rows[i];
        
        // get vector of current quality
        compute_quality_ratio(qold,R.cols);
        
        // try +1
        if (increment_solution(R,1))
        {
          compute_quality_ratio(qnew,R.cols);
          if (1==is_better(qnew,qold))
          {
            improved = true;
            //          result->info_message("fine tune made an improvement using row %lu = ", (unsigned int) i);
            //          R.print_row(result);
            continue;
          }
        }
        
        // try -1
        if (increment_solution(R,-2))
        {
          compute_quality_ratio(qnew,R.cols);
          if (1==is_better(qnew,qold))
          {
            improved = true;
            //          result->info_message("fine tune made an improvement using row %lu = ", (unsigned int) i);
            //          R.print_row(result);
            continue;
          }
        }
        
        // restore solution to +0
        increment_solution(R,1);
      }
    } while (improved);
  }
  
  
  void IncrementalIntervalAssignment::categorize_vars(vector<int>&col_order,vector<int>&cols_dep, vector<int>&cols_ind, vector<int>&rows_dep)
  {
    cols_dep.reserve(used_row+1); // can't be more *dependent* variables than rows.
    cols_ind.reserve(col_order.size()); // all vars could be independent
    rows_dep.reserve(used_row+1);
    
    vector<int>col_is_dependent(used_col+1,0);
    int first_empty_row=used_row+2; // big value
    for (int r = 0; r <= used_row; ++r)
    {
      // this is the only useful thing this method does!
      // otherwise, the first used_row entries of col_order are the dependent variables
      auto &rc=rows[r].cols;
      if (rc.empty())
      {
        if (r<first_empty_row)
          first_empty_row=r;
        continue;
      }
      
      rows_dep.push_back(r);
      int c = col_order[r];
      cols_dep.push_back(c);
      col_is_dependent[c]=1;
      
      // sanity check that c is correct
      assert(c<=used_col);
      assert(col_rows[c].size()==1);
      assert(col_rows[c][0]==r);
      assert(r<first_empty_row); // index c is off if we have some blank rows and then some non-blank ones. If this happens, need to row-swap earlier to fix.
    }
    // all other columns are dependent
    for (int c = 0; c < (int) col_is_dependent.size(); ++c)
    {
      if (col_is_dependent[c]==0)
        cols_ind.push_back(c);
    }
    // sanity check
    assert(cols_dep.size()==rows_dep.size());
    assert(cols_dep.size()+cols_ind.size()==col_is_dependent.size());
    
    if (result->log_debug)
    {
      result->info_message("Categorize_vars\n");
      result->info_message("Dependent rows: ");   print_vec(result, rows_dep);
      result->info_message("Dependent vars: ");   print_vec(result, cols_dep);
      result->info_message("Independent vars: "); print_vec(result, cols_ind);
      result->info_message("\n");
    }
  }
  
  bool IncrementalIntervalAssignment::adjust_independent_vars(vector<int>&cols_fail)
  {
    // heuristic: search for another var that only appears in this row, and increment it until we can make the sum divisible by it
    // won't work if there is no such row
    vector<int>still_fail;
    for (auto c : cols_fail)
    {
      bool fixed=false;
      if (col_rows[c].size()==1)
      {
        int r = col_rows[c][0];
        for (size_t i = 0; i < rows[r].cols.size(); ++i)
        {
          auto cc = rows[r].cols[i];
          if (cc!=c && col_rows[cc].size()==1)
          {
            int &cv  = col_solution[c];
            int &ccv = col_solution[cc];
            cv=0;  // the independent var
            ccv=0; // the dependent var, we will re-assign this later regardless
            int sum=row_sum(r);
            int c_coeff = rows[r].vals[0];
            int cc_coeff = rows[r].vals[i];
            // increase cc var until sum is divisible by c's coefficient
            for (int j=0; j < c_coeff && !fixed; ++j)
            {
              if ((sum + cc_coeff*ccv) % c_coeff == 0)
              {
                cv=(sum + cc_coeff*ccv)/c_coeff;
                assert(row_sum(r)==0);
                fixed=true;
                break;
              }
              else
                ++ccv;
            }
          }
        }
      }
      if (!fixed)
      {
        still_fail.push_back(c);
      }
    }
    
    cols_fail.swap(still_fail);
    return cols_fail.empty() ? true : false;
  }
  
  bool IncrementalIntervalAssignment::generate_tied_data()
  {
    if (!has_tied_variables)
      return true;
    
    for ( int c=0; c < (int) tied_variables.size(); ++c)
    {
      const auto sz = tied_variables[c].size();
      if (sz<2)
        continue;
      
      auto &dt = tied_data[c];
      dt.orig_bound_lo=col_lower_bounds[c];
      dt.orig_bound_hi=col_upper_bounds[c];
      dt.goal_lo=goals[c];
      dt.goal_hi=dt.goal_lo;
      
      int hi=dt.orig_bound_hi;
      int lo=dt.orig_bound_lo;
      for (int i=1; (size_t) i<tied_variables[c].size(); ++i)
      {
        int c2 = tied_variables[c][i];
        lo = max(lo,col_lower_bounds[c2]);
        hi = min(hi,col_upper_bounds[c2]);
        
        auto g = goals[c2];
        if (g==0.)
          continue;
        if (g < dt.goal_lo || dt.goal_lo==0.)
          dt.goal_lo = g;
        if (g > dt.goal_hi || dt.goal_hi==0.)
          dt.goal_hi = g;
      }
      if (lo>hi)
      {
        result->warning_message("Infeasible interval bounds: interval bound ranges do not overlap for intervals that are constrained to be identical.\n");
        // to do, print something about the ids of the entities for the user
        if (result->log_debug)
        {
          result->info_message("Tied variables = ");
          print_vec(result, tied_variables[c]);
        }
        return false;
      }
      else
      {
        col_lower_bounds[c]=lo;
        col_upper_bounds[c]=hi;
      }
    }
    return true;
  }
  
  bool IncrementalIntervalAssignment::bounds_unsatisfied(int c)
  {
    // for tied variables, checks are on the tied-to variable instead
    if (c < ((int) tied_variables.size()) && tied_variables[c].size()==1)
      return false;
    return (col_solution[c]<col_lower_bounds[c] || col_solution[c]>col_upper_bounds[c]);
  }
  
  bool IncrementalIntervalAssignment::bounds_satisfied(vector<int>&cols_fail)
  {
    cols_fail.clear();
    for ( int c=0; c <= used_col; ++c)
    {
      if (bounds_unsatisfied(c))
        cols_fail.push_back(c);
    }
    return cols_fail.empty();
  }

  void MatrixSparseInt::row_simplify(int r)
  {
    auto &v = rows[r].vals;
    row_simplify_vals(v);
  }

  void IncrementalIntervalAssignment::row_simplify(int r)
  {
    auto &v = rows[r].vals;
    if (rhs[r]!=0)
    {
      v.push_back(rhs[r]);
      row_simplify_vals(v);
      rhs[r] = v.back();
      v.pop_back();
    }
    else
    {
      row_simplify_vals(v);
    }
  }
  
  bool IncrementalIntervalAssignment::create_nullspace(const vector<int>&cols_dep,
                                                       const vector<int>&cols_ind,
                                                       const vector<int>&rows_dep,
                                                       MatrixSparseInt&M,
                                                       int &MaxMrow)
  {
    // M is #independent_vars rows, by #variables columns
    // these two lines are only used if we are trying to append to an existing nullspace
    // auto Mrow = M.rows.size();
    // const auto MrowStart = Mrow;
    size_t Mrow = 0;
    M.rows.resize(Mrow + cols_ind.size());
    M.col_rows.resize(cols_dep.size()+cols_ind.size());
    
    // calculate least-common-multiple of the diagonal coefficients
    // Implicitly, as we assemble the nullspace, each row of A is multiplied by lcm / leading value of that row
    // A *= lcm / leading value of each row
    int L=1;
    for (size_t row_i = 0; row_i <rows_dep.size(); ++row_i)
    {
      const auto r = rows_dep[row_i];
      assert(r==(int)row_i);
      const auto c = cols_dep[row_i];
      const auto v = abs(rows[r].get_val(c));
      assert(v!=0);
      L=lcm(L,v);
    }
    
    result->debug_message("creating nullspace starting in row %d\n",(int)Mrow);
    result->debug_message("lcm=%d\n",L);
    
    if (L==0)
      return false;
    
    for (size_t i = 0; i < cols_ind.size(); ++i, ++Mrow)
    {
      // i is the index into cols_ind
      // Mrow is the row of M we are filling in
      
      // c = column of this we are turning into a row vector of the nullspace
      auto c = cols_ind[i];
      
      auto &row = M.rows[Mrow];
      
      // result->info_message("Filling in row %d from matrix column %d\n",(int)i,c);
      // grab the entries for the independent vars
      // transform from row index to cols_dep index
      // row.cols=col_rows[c];
      row.cols.resize(col_rows[c].size());
      for (size_t j=0; j<col_rows[c].size(); ++j)
      {
        row.cols[j]=cols_dep[ col_rows[c][j] ];
      }
      row.vals=column_coeffs(c);
      // result->info_message("intial "); row.print_row(result);
      
      // multiply by -1 * L / leading value of row
      for (size_t j = 0; j < row.vals.size(); ++j)
      {
        // a non-zero entry
        // row corresponding to row.vals[j] non-zero entry == column it came from
        const auto cc = col_rows[c][j];
        const auto rr=rows_dep[cc];
        const auto lead=rows[rr].get_val( cols_dep[rr] );
        auto &vv=row.vals[j];
        // result->info_message("  old value = %d\n",vv);
        assert( abs(L) % abs(lead) == 0 );
        vv*= -L / lead;
        
        //      result->info_message("j=%d, row.vals[%d]=%d",(int)j,(int)j, vv);
        //      result->info_message("  permuted col=%d, original matrix row index=%d=%d ", c, cc, rr); rows[rr].print_row(result);
        //      result->info_message("  lead=rows[%d].get_val(%d)=%d\n",rr,cols_dep[rr],lead);
        //      result->info_message("  new value = %d\n",vv);
      }
      
      // result->info_message("after multiply by row_dep / leading value of "); row.print_row(result);
      
      // grap the entries for the dependent vars
      // append column of identity matrix, times L
      // row.cols.push_back(cols_dep.size() + i);
      row.cols.push_back(cols_ind[i]);
      row.vals.push_back(L);
      row.sort();
      
      row_simplify_vals(row.vals);
      // result->info_message("after simplify "); row.print_row(result);
      // result->info_message("\n");
      
//      // if it is identical to a passed-in vector, then erase it
//      for (int m = 0; m < MrowStart; ++m)
//      {
//        if (row.cols == M.rows[m].cols && row.vals == M.rows[m].vals)
//        {
//          result->debug_message("Clearing identical row %d because it is the same as passed-in row %d\n",Mrow,m);
//          row.cols.clear();
//          row.vals.clear();
//          --Mrow;
//          M.rows.resize( M.rows.size()-1 );
//          break;
//        }
//      }
    }
    M.fill_in_cols_from_rows();
    MaxMrow=(int)M.rows.size();
    return true;
  }
  
  // Partion A and B: C = A cap B.  A = A \ B.  B = B \ A.
  void remove_intersection( set<int> &A, set<int> &B, set<int> &C)
  {
    C.clear();
    if (A.empty() || B.empty())
      return;
    set<int> AnotB, BnotA;
    auto a=A.begin();
    auto b=B.begin();
    while (true)
    {
      if (*a<*b)
      {
        AnotB.emplace_hint(AnotB.end(),*a);
        ++a;
      }
      else if (*a>*b)
      {
        BnotA.emplace_hint(BnotA.end(),*b);
        ++b;
      }
      else // *a==*b
      {
        C.emplace_hint(C.end(),*a);
        ++a;
        ++b;
      }
      if (a==A.end())
      {
        while (b!=B.end())
        {
          BnotA.emplace_hint(BnotA.end(),*b);
          ++b;
        }
        A.swap(AnotB);
        B.swap(BnotA);
        return;
      }
      if (b==B.end())
      {
        while (a!=A.end())
        {
          AnotB.emplace_hint(AnotB.end(),*a);
          ++a;
        }
        A.swap(AnotB);
        B.swap(BnotA);
        return;
      }
    } // while
  }
  
  bool MatrixSparseInt::coefficient_irreducible(int col, int max_coeff)
  {
    // are entries for col some multiple of max_coeff?
    // this will catch the sum-even variables that started with only a "2" coefficient
    // But, it will give up in some cases where it would be possible to reduce the coefficient
    // e.g.
    //    2  1  1
    //  + 2 -1  3
    // -----------
    //    4  0  4  --->  1 0 1
    
    auto &cr=col_rows[col];
    if (cr.empty())
      return true;
    int g = abs(rows[ cr[0] ].get_val(col));
    for (size_t i=1; i<cr.size(); ++i)
    {
      g = gcd(g, rows[ cr[i] ].get_val(col));
    }
    if (g>=max_coeff)
    {
      result->debug_message("The gcd of all of x%d's coefficients is too large: %d >= %d.\n", col, g, max_coeff);
      return true;
    }
    
    // Are we stuck?
    // Stuck if there is some variable k with gcd(coeff(k),c_coeff)==1 and it only appears in the same ratio everywhere
    set<int> checked_cols;
    for (auto r : cr )
    {
      for (auto c : rows[r].cols)
      {
        if (c == col)
          continue;
        auto p = checked_cols.insert(c);
        if (!p.second) // already checked
          continue;
        
        // check c against col
        // is c not a blocker in some row?
        bool is_ok = false;
        for (int r : col_rows[c])
        {
          int vcol = rows[r].get_val(col);
          int vc   = rows[r].get_val(c);
          if (vcol==0 || vc==0)
          {
            is_ok=true;
            break;
          }
          int g = gcd(vcol,vc);
          if (abs(vcol/g)<max_coeff)
          {
            // variable doesn't block reducing col coefficient in this row
            is_ok=true;
            break;
          }
        }
        if (!is_ok)
        {
          // variable blocks reducing col in every row
          // irreducible if it is a blocker in the same ratio in every row.
          int old_vcol=0, old_vc=0;
          for (int r : col_rows[c])
          {
            int vcol = rows[r].get_val(col);
            int vc   = rows[r].get_val(c);
            if (vcol<0)
            {
              vcol=-vcol;
              vc=-vc;
            }
            int g = gcd(vcol,vc);
            vcol /= g;
            vc /= g;
            if (old_vcol==0)
            {
              old_vcol=vcol;
              old_vc=vc;
            }
            else
            {
              if (old_vcol!=vcol || old_vc!=vc)
              {
                is_ok=true;
                break;
              }
            }
          }
          if (!is_ok)
          {
            result->debug_message("The ratio of x%d's coefficients is x%d is %d to %d in every row.\n", col, c, old_vcol, old_vc);
            return true;
          }
        }
      }
    }
    return false;
  }
  
  void find_new_blocking_col( const BlockingCols &blocking_cols, const BlockingCols &all_blocking_cols, vector<int> &new_blocking_cols)
  {
    for (auto bc : blocking_cols)
    {
      int b_col = bc.first;
      if (all_blocking_cols.find(b_col) == all_blocking_cols.end())
      {
        new_blocking_cols.push_back(b_col);
      }
    }
  }
  
  void MatrixSparseInt::unreduce(int &r, int c, BlockingCols &all_blocking_cols)
  {
    // if it is in all_blocking_cols, then remove it. Else done.
    auto it = all_blocking_cols.find(c);
    if (it == all_blocking_cols.end())
      return; // wasn't in all_blocking_cols
    all_blocking_cols.erase(it);
    
    // decide whether to unreduce its row
    if (r<1)
      return; // no rows were reduced to begin with
    
    // if it is in more than one row, then it was never successfully reduced, so do nothing
    if (col_rows[c].size()!=1)
      return;
    
    // if it isn't the only blocker in row rr, i.e. the row is rref but not diagonal, then remove it but don't mark the row as unreduced
    const int rr = col_rows[c][0];
    for (int k : rows[rr].cols)
    {
      if (k!=c && all_blocking_cols.find(k)!=all_blocking_cols.end())
      {
        return;
      }
    }
    
    // it is the only blocker in a row, i.e. its row is diagonal and not rref, so remove that row from the ones that are reduced
    --r;
    row_swap(r,rr);
  }
  
  bool MatrixSparseInt::gaussian_elimination(int col, const BlockingCols &blocking_cols, int max_coeff, map<int,int> &improved_cols,
                                             BlockingCols &all_blocking_cols, int &r, int &unblocked_row)
  {
    // if work_hard is true, then we switch which column we reduce when there is a choice and we can't reduce all of them
    //   it doesn't seem to add significant time, but keeps the algorithm from getting stuck sometimes
    //   if it is a time bottleneck, we could try only doing it only when solving constraints, or when other options failed
    bool work_harder = true;
    
    // debug
    if (result->log_debug)
    {
      print_matrix("Starting Gaussian elimination");
      // print_matrix_summary("Starting Gaussian elimination");
      
      result->debug_message("Trying to diagonalize %lu columns ", (long unsigned) blocking_cols.size());
      print_map(blocking_cols);
      
      result->debug_message("%d rows of %d already diagonalized\n", r, (int) rows.size());
      result->debug_message("Already diagonalize %lu columns ", (long unsigned) all_blocking_cols.size());
      print_map(all_blocking_cols);
      result->debug_message("Since last call we've improved %lu columns ", (long unsigned) improved_cols.size());
      ::IIA_Internal::print_map(result,improved_cols);
    }
    // debug
    
    // sanity check, col itself isn't a blocker; that is handled by self-block instead
    assert(blocking_cols.find(col)==blocking_cols.end());
    
    // erase improved_cols that aren't actually changed, or are a new blocker
    if (!all_blocking_cols.empty() && !improved_cols.empty())
    {
      for (auto it = improved_cols.begin(); it != improved_cols.end(); )
      {
        if (it->second==0  ||  blocking_cols.find(it->first)!=blocking_cols.end())
        {
          it = improved_cols.erase(it);
        }
        else
        {
          ++it;
        }
      }
      // remove remaining improved_cols from all_blocking_cols
      //   update all_blocking_cols and r
      // undiagonalize col (if it is diagonalized)
      if (!improved_cols.empty())
      {
        const auto old_all_block_size = all_blocking_cols.size();
        for (auto imp : improved_cols)
        {
          unreduce(r, imp.first, all_blocking_cols);
        }
        
        if (result->log_debug)
        {
          auto num_removed = old_all_block_size - all_blocking_cols.size();
          result->info_message("all_blocking_cols smaller by %lu\n", (unsigned long) num_removed);
          if (num_removed)
          {
            print_matrix("gaussian ");
            // print_matrix_summary("gaussian ");
            
            result->debug_message("%d rows of %d remain diagonalized\n", r, (int) rows.size());
            result->debug_message("remaining diagonalize %lu columns ", (long unsigned) all_blocking_cols.size());
            print_map(all_blocking_cols);
          }
        }
      }
    }
    improved_cols.clear();
    
    // find which blocking cols are new; these are the ones to reduce
    vector<int> new_blocking_cols;
    find_new_blocking_col( blocking_cols, all_blocking_cols, new_blocking_cols);
    if (!new_blocking_cols.empty())
    {
      
      // unreduce col, so we can reduce any blockers in its row instead
      unreduce(r, col, all_blocking_cols);
      
      // to research: should we add new_blocking_cols to all_blocking_cols regardless of whether we can reduce them?
      // for now, we only store the eliminated ones, but this might confuse the caller
      
      unblocked_row=-1; // row we can increment by
      for (int lead : new_blocking_cols) // col we're working on
      {
        result->debug_message("reducing row %d col %d\n",r,lead);
        
        // out of rows to reduce
        if ((size_t)r==rows.size())
        {
          result->debug_message( "Perdu! Stuck! We already reduced every row!"); // diagonalize a different row here, too
          break; // but still search for an OK row, as is we succeeded
        }
        
        // goal is this(r,lead)=1
        // find a row with lead coeff == 1 (or smallest) to pivot to r
        int rr=-1, rr0=-1, rr1=-1, rr2=-1; // the row we will use, and its backups in case there isn't anything ideal
        int backup_coeff = numeric_limits<int>::max();
        // search starting from the last rows, as they are the unreduced ones...
        auto &col_r=col_rows[lead];
        for (int i = ((int) col_r.size())-1; i>=0; --i)
        {
          int s = col_r[i];
          if (s < r)
          {
            // no good, we've already fixed this row, so it has other leading nonzeros
            // but are those other rows still relevant?
            break; // no more rows to consider, since we started with the last row and worked backwards
          }
          else
          {
            int lead_coeff = rows[s].get_val(lead);
            if (abs(lead_coeff)==1)
            {
              // prefer a row without col
              if (rows[s].get_val(col))
                rr0=s;
              else
              {
                rr=s; // perfect
                break;
              }
            }
            else
            {
              if (abs(lead_coeff)<backup_coeff)
              {
                backup_coeff=abs(lead_coeff);
                // prefer a row without col
                if (rows[s].get_val(col))
                  rr2=s;
                else
                  rr1=s;
              }
            }
          }
        }
        
        // pick the best row we found
        if (rr==-1)
        {
          if (rr0>-1)
          {
            rr=rr0; // coefficient==1, has col
          }
          else
          {
            if (rr1>-1)
            {
              rr=rr1; // coefficient!=1, no col
            }
            else
            {
              if (rr2>-1)
              {
                rr=rr2; // coefficient!=1, has col
              }
              else
              {
                // lead doesn't exist in any row rr >= r
                // if it is in only one row, then this is equivalent to it being diagonalized, so success!
                // remember ... RREF not diagonalized form
                if (col_rows[lead].size()==1)
                {
                  all_blocking_cols.insert(*blocking_cols.find(lead));
                  continue; // done with lead, move on to next new blocker
                }
                
                // if we got to here then we didn't find any good rows >r, and new blocker lead is unreduced
                // if lead is in two rows, then we could choose which 2 of the 3 variables to diagonalize...
                if (work_harder)
                {
                  // diagonalize lead instead of some col that isn't in blocking_cols, if any row has such a case
                  bool fixed = false;
                  for (int i = ((int) col_r.size())-1; i>=0; --i)
                  {
                    auto rb = col_r[i];
                    // check if lead is the only blocker in the row
                    bool lead_only_blocker=true;
                    for ( auto cb : rows[rb].cols )
                    {
                      if (cb!=lead)
                      {
                        if (blocking_cols.find(cb)!=blocking_cols.end())
                        {
                          lead_only_blocker=false;
                        }
                      }
                    }
                    if (lead_only_blocker)
                    {
                      // use this row to eliminate lead
                      matrix_allrows_eliminate_column(rb, lead);
                      all_blocking_cols.insert(*blocking_cols.find(lead));
                      // don't increment r, since we modified an already-reduced row instead
                      fixed = true;
                      break; // success
                    }
                  }
                  if (fixed) continue; // success
                                       // else we didn't find any row where lead was the only blocker, give up
                }
                
                // give up, we aren't going to reduce lead
                continue; // move on to the next new blocker
              }
            }
          }
        }
        
        // do the reduce on col lead and row rr
        result->debug_message("matrix_row_eliminate using row %d to kill col %d\n",rr, lead);
        
        // pivot
        row_swap(rr,r);
        
        matrix_allrows_eliminate_column(r, lead);
        
        //        if (result->log_debug)
        //        {
        //          result->debug_message("Done reducing row %d and col %d\n",r,lead);
        //          print_matrix("done row col");
        //          // print_matrix_summary("done row col");
        //        }
        
        // sanity check
        assert( col_rows[lead].size()==1 && col_rows[lead][0]==r);
        if (col_rows[lead].size()!=1 || col_rows[lead][0]!=r)
        {
          result->error_message("IIA: Failed Gaussian elimination\n");
          if (result->log_debug)
          {
            print_matrix("Failed Gaussian elimination");
            result->debug_message("  Gaussian elimination must have an implementation error. Failed to remove column %d from all rows except %d.\n",lead,r);
          }
          return false;
        }
        
        // success, next row/col
        // col_order.push_back(lead);
        // update all_blocking_cols
        all_blocking_cols.insert(*blocking_cols.find(lead));
        ++r;
      }
    }
    
    if (result->log_debug)
    {
      result->debug_message("Gaussian elimination finished\n");
      print_matrix("Gaussian elimination finished");
    }
    
    // look for a good row to return, regardless of whether we actually reduced anything
    set<int> rows_blocked, rows_large_coeff, rows_both;
    
    // is there a good row with 0 for the blocking coefficients?
    bool print_first=true;
    for (int i = ((int) col_rows[col].size()) - 1; i>=0; --i)
    {
      auto cr = col_rows[col][i];
      auto &rowcr=rows[cr];
      row_simplify_vals(rowcr.vals); // already done??
      
      // these are the rows with no blocking variable coefficient, which we check first
      if (cr>=r)
      {
        int c_coeff=abs(rowcr.get_val(col));
        if (c_coeff<max_coeff)
        {
          if (result->log_debug)
          {
            result->debug_message("Row %d has zero coefficients for all the blocking columns and small x%d coefficient %d ",cr, col, c_coeff);
            rows[cr].print_row(result);
          }
          
          // return first unblocked row
          unblocked_row=cr;
          return true; // success // option A
        }
        else
        {
          rows_large_coeff.insert(cr);
          if (result->log_debug)
          {
            result->debug_message("Row %d has a big x%d coefficient %d>=%d",cr,col,c_coeff,max_coeff);
            rowcr.print_row(result);
          }
        }
      }
      else
        // if (cr<r)  // these are the rows with a blocking variable coefficient; could be more than one if rows < block_cols.size()
      {
        // debug
        if (print_first)
        {
          result->debug_message("every row with small col %d coefficient had a non-zero coefficient for some blocking column\n",col);
          print_first=false;
        }
        // debug
        
        int c_coeff=rowcr.get_val(col);
        if (is_row_blocked(cr,col,blocking_cols))
        {
          rows_blocked.insert(cr);
        }
        else
        {
          // unblocked
          if (abs(c_coeff)<max_coeff)
          {
            if (result->log_debug)
            {
              result->debug_message("Row %d has acceptable (but not zero) coefficients for the blocking columns and small %d coefficient %d ",cr, col, c_coeff);
              rows[cr].print_row(result);
            }
            
            // save first unblocked row
            unblocked_row=cr;
            return true; // option A
          }
        }
      }
    }
    
    // if (unblocked_row!=-1)
    //  return true; // option B
    
    remove_intersection( rows_blocked, rows_large_coeff, rows_both );
    
    // are we stuck?
    // if all rows have a large coefficient, is the gcd too big to be able to reduce it small enough?
    if (rows_blocked.empty())
    {
      if (coefficient_irreducible(col,max_coeff))
      {
        result->debug_message("It's impossible to reduce the coefficient of x%d below %d.\n", col, max_coeff);
        return false;
      }
    }
    
    // ==============  below here is heuristics for reducing the leading coefficient, getting c_coeff < max_coeff ===========
    
    for (auto cr : rows_large_coeff)
    {
      if (cr<r)
      {
        if (reduce_coefficient(cr,col,r,max_coeff))
        {
          if (result->log_debug)
          {
            int c_coeff = abs(rows[cr].get_val(col));
            result->debug_message( "Final coefficient reduction success! row %d x%d coeff %d<%d ", cr, col, c_coeff, max_coeff);
            rows[cr].print_row(result);
          }
          unblocked_row=cr;
          return true;
        }
      }
      else
      {
        result->debug_message("IIA: Gaussian elimination to both reduce_coefficient and fix blocking cols is not implemented.\n");
      }
    }
    
    // debug
    if (result->log_debug)
    {
      if (max_coeff<numeric_limits<int>::max())
      {
        result->debug_message("IIA: Gaussian elimination unable to find an unblocked vector with a small enough coefficient (<%d) for x%d.\n",max_coeff, col);
      }
      else
      {
        result->debug_message("IIA: Gaussian elimination unable to find an unblocked vector for x%d.\n", col);
      }
    }
    // debug
    
    return false;
  }
  
  
  size_t MatrixSparseInt::blocking_size(BlockingCols &blocking_cols)
  {
    size_t count(0);
    for (auto bc : blocking_cols)
    {
      if (bc.second.first!=block_min)
        ++count;
      if (bc.second.second!=block_max)
        ++count;
    }
    return count;
  }
  
  bool MatrixSparseInt::is_blocked(int v, pair<int,int> block)
  {
    assert(block.first<0);
    assert(block.second>0);
    return (v<=block.first || v>=block.second);
  }
  
  bool MatrixSparseInt::is_row_blocked(int row, int col, const BlockingCols &blocking_cols )
  {
    // check that the non-zero entries are not blocked
    auto &rowcr=rows[row];
    int c_coeff=rowcr.get_val(col);
    
    for (size_t j=0; j<rowcr.cols.size(); ++j)
    {
      int cc = rowcr.cols[j];
      auto bcc=blocking_cols.find(cc);
      if (bcc!=blocking_cols.end())
      {
        int cc_coeff=rowcr.get_val(cc);
        // flip sign if c_coeff is negative
        if (c_coeff<0)
          cc_coeff = -cc_coeff;
        
        if (is_blocked(cc_coeff,bcc->second))
        {
          if (result->log_debug)
          {
            result->debug_message("Row %d is blocked by x%d coefficient %d ",row, cc, cc_coeff);
            const auto lo = bcc->second.first;
            const auto hi = bcc->second.second;
            if (lo==block_min)
            {
              result->debug_message(" >= %d\n", hi);
            }
            else if (hi==block_max)
            {
              result->debug_message(" <= %d\n", lo);
            }
            else
            {
              result->debug_message("outside (%d,%d)\n", lo, hi);
            }
            rowcr.print_row(result);
          }
          return true;
        }
      }
    }
    if (result->log_debug)
    {
      result->debug_message("Row %d has acceptable (but perhaps non-zero) coefficients for all the blocking columns", row);
      rowcr.print_row(result);
    }
    return false; // success, it isn't blocked
  }
  
  bool MatrixSparseInt::reduce_coefficient(int row, int col, int r, int max_coeff)
  {
    // there is a chance, but no guarantee, we can make the coefficient smaller, by combining integer combinations of rows
    
    if (result->log_debug)
    {
      auto &roww=rows[row];
      int c_coeff=roww.get_val(col);
      result->debug_message("IIA: Gaussian elimination removed the blocking variables from row %d, but the x%d coefficient is too big, %d>=%d.\n",row,col,abs(c_coeff),max_coeff);
      rows[row].print_row(result);
    }
    
    // If this line is hit, then it would make a good test case for gaussian_elimination to reduce coefficients.
    result->warning_message("IIA: Gaussian elimination, getting the leading coefficient smaller may be possible, but this is not implemented yet.\n");
    result->warning_message("Interval Matching: encountered a case that is theoretically possible but not encountered before in practice. Please share this example journal file to Scott A. Mitchell.\n");
    
    return false; // not implemented, because we have no test cases that need it
    
    // rest of method has been deleted - it wasn't finished.
    // See repo version of IncrementalIntervalAssignment prior to 13 Dec 2019 to recover it.
  }
  
  int MatrixSparseInt::push_row(const RowSparseInt &row)
  {
    const int new_r = (int)rows.size();
    rows.push_back( row );
    for (auto c : rows[new_r].cols)
      col_rows[c].push_back(new_r);
    return new_r;
  }
  
  vector<int> IncrementalIntervalAssignment::column_coeffs(int c) const
  {
    vector<int> vals;
    auto &col = col_rows[c];
    for (size_t i = 0; i < col.size(); ++i)
    {
      auto r = col[i];
      vals.push_back(get_M(r,c));
    }
    return vals;
  }
    
  // For an explanation of Reduced Row Echelon Form, see
  // https://people.sc.fsu.edu/~jburkardt/c_src/row_echelon_integer/row_echelon_integer.html
  //
  // For how this helps us compute the nullspace of a maxtrix (multiply), see
  // https://math.stackexchange.com/questions/88301/finding-the-basis-of-a-null-space
  //
  // For how that helps us solve the Integer optimization problem, see
  // https://cs.stackexchange.com/questions/82086/linear-programming-restricted-to-rational-coefficients
  // and
  // Slide 27, totally unimodular ILP solve (which we don't quite have)
  // https://www.isical.ac.in/~arijit/courses/autumn2016/ILP-Lecture-1.pdf
  
  bool IncrementalIntervalAssignment::rref_elim(int &r, int rr, int c, vector<int> &rref_col_order, vector<int> &rref_col_map)
  {
    // result->debug_message(" rref filling in row %d (diagonal is col %d) ", r, c );
    
    row_swap(r,rr); // now we're working on r, not rr
    
    // mark column as used
    rref_col_map[c]=(int)rref_col_order.size();
    rref_col_order.push_back(c);
    
    // eliminate c from all other rows
    allrows_eliminate_column(r,c); // plus also update the rhs
    
    // ensure c coefficient is positive, not negative (e.g. 1, not -1)
    if (rows[r].get_val(c)<0)
    {
      row_negate(r);
    }
    
    // next row
    ++r;
    
    // debug
    // sanity check, was the column eliminated?
    // For general matrics there is at most one non-zero, and for IA matrices exactly one non-zero
    if (col_rows[c].size()!=1)
    {
      result->error_message("IIA: RREF failure, unable to get a single non-zero column for dependent variable %d\n",c);
      result->info_message("col %d rows = ",c);
      print_vec(result, col_rows[c]);
      print_problem("RREF stuck case");
      return false;
    }
    assert( col_rows[c].size()==1 );
    // debug
    
    return true;
  }
  
  void RowSparseInt::sparse_to_full(vector<int> &full) const
  {
    if (!cols.empty() && cols.back() >= full.size())
    {
      full.resize(cols.back()+1);
    }
    int j = 0;
    for (int c=0; c<full.size(); ++c)
    {
      if (j<cols.size() && cols[j] == c)
      {
        full[c] = vals[j];
        ++j;
      }
      else
      {
        full[c] = 0;
      }
    }
  }

  int RowSparseInt::get_col_index(int c) const
  {
    for (int i=0; i<(int)cols.size(); ++i)
      if (cols[i]==c)
        return i;
    assert(0);
    return -1;
  };
  
  bool IncrementalIntervalAssignment::rref_step0(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map)
  {
    result->debug_message("rref_step0\n");
    // Collect vars that are already reduced
    // Use any var that has a "1" for a coefficient, and there is only one row that contains that variable
    priority_queue<QElement> Q;
    SetValuesOneCoeff val_OneCoeff_Q;
    vector<int>all_vars(used_col+1);
    iota(all_vars.begin(),all_vars.end(),0);
    const double threshold=0.1; // 0 is bad, 1 is good
    build_Q(Q, val_OneCoeff_Q, threshold, all_vars);
    
    while (!Q.empty())
    {
      // no need for tip_top or re-prioritizing, could just use a static sorted vector
      QElement t = Q.top();
      Q.pop();
      
      const int c = t.c;
      const int cr = col_rows[c][0];
      assert( col_rows[c].size()==1 );
      
      if (cr<rref_r) // if var is in a lower-ranked row, then its an independent var of that row already, and can't be a dependent var
        continue;
      
      assert( abs(rows[cr].get_val(c)) == 1);
      if (!rref_elim(rref_r, cr, c, rref_col_order, rref_col_map))
        return false;
    }
    return true;
  }

  bool IncrementalIntervalAssignment::rref_step_numrows(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map)
  {
    result->debug_message("rref_step_numrows\n");
    priority_queue<QElement> Q;
    SetValuesNumCoeff val_NumCoeff;
    vector<int>all_vars(used_col+1);
    iota(all_vars.begin(),all_vars.end(),0);
    const double threshold=numeric_limits<double>::lowest(); // anything is fine
    build_Q(Q, val_NumCoeff, threshold, all_vars);
    
    QElement t;
    while (!Q.empty())
    {
      if (tip_top(Q, val_NumCoeff, threshold, t))
        return true;
      
      const int c = t.c;
      // find an unreduced row with a smallest coefficient
      int cr = -1;
      int best_coeff = 0;
      for ( auto r : col_rows[c] )
      {
        if (r>=rref_r)
        {
          row_simplify(r);
          auto coeff = abs(rows[r].get_val(c));
          if (cr==-1 || coeff<best_coeff)
          {
            cr=r;
            best_coeff=coeff;
          }
        }
      }
      if (cr==-1) // e.g. -1 means no unreduced row had a non-zero coeff, can't use this variable
        continue;
      
      assert(cr>=rref_r);
      if (!rref_elim(rref_r, cr, c, rref_col_order, rref_col_map))
        return false;
      row_simplify(rref_r-1);
      
      if (rref_r>used_row)
        return true;
    }
    return true;
  }

  
  bool IncrementalIntervalAssignment::rref_step1(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map, bool map_only_phase)
  {
    if (rref_r>used_row)
      return true;
    
    result->debug_message("rref_step1\n");
    
    SetValuesCoeffRowsGoal val_CoeffRowsGoal_Q; // map_only_phase
    SetValuesCoeffRGoal val_CoeffRGoal_Q; // even phase
    QWithReplacement Q(this, val_CoeffRowsGoal_Q, numeric_limits<double>::lowest()/4);
    if (use_better_rref_step1_pivots_in_sumeven_phase && !map_only_phase) // zzyk, test doing this in the map_only_phase as well
    {
      Q.set_val_fn(val_CoeffRGoal_Q);
    }
    vector<int>all_vars;
    all_vars.reserve(used_col+1);
    for (int c=0; c<=used_col; ++c)
    {
      if (rref_col_map[c]==-1)
        all_vars.push_back(c);
    }

    // workspace
    set<int>   changed_var_set;
    vector<int> update_var;

    for (;;)
    {
      Q.build(all_vars);
      // Q.print(); zzyk
      
      if (Q.empty())
        return true;
      
      all_vars.clear();
      
      result->debug_message("Trying %d candidate columns\n",(int)Q.size());
      
      bool progress_since_pushback=false;
      QElement t;
      while (!Q.empty() && rref_r<=used_row)
      {
        if (Q.tip_top(t))
          break;
        
        const int c = t.c;

        // skip if already reduced, isn't in any row, or all the rows its in are reduced already
        if (rref_col_map[c]>0 || col_rows[c].empty() || col_rows[c].back() < rref_r)
          continue;

        // result->info_message("processing t, col %d, valA %g, valB %g, valC %g\n",t.c,t.valueA,t.valueB,t.valueC);
        
        // rr row to move into r
        // find the one with a lead coeff of 1 if any, and prefer the shortest row
        int rr=-1;
        for (int rrr : col_rows[c])
        {
          if (rrr<rref_r)
            continue;
          if (abs(rows[rrr].get_val(c))==1)
          {
            if (rr==-1 || rows[rrr].cols.size() < rows[rr].cols.size())
            {
              rr=rrr;
            }
          }
        }
        if (rr==-1)
        {
          // we may just need to row_simplify
          for (int rrr : col_rows[c])
          {
            if (rrr<rref_r)
              continue;
            row_simplify(rrr);
            if (abs(rows[rrr].get_val(c))==1)
            {
              if (rr==-1 || rows[rrr].cols.size() < rows[rr].cols.size())
              {
                rr=rrr;
              }
            }
          }
        }
        if (rr==-1)
        {
          // try again later
          
          // debug
          if (result->log_debug)
          {
            result->debug_message("skipping %d in %d rows for now, no 1 coefficient in a pivot row anymore\n",c,(int)col_rows[c].size());
            for (auto s : col_rows[c])
            {
              if (s<rref_r)
                result->info_message("*");
              else if (s==rr)
                result->info_message("+");
              else
                result->info_message("-");
              print_row_iia(s);
            }
          }
          //debug
          
          if (all_vars.empty())
            progress_since_pushback=false;
          all_vars.push_back(c);
        }
        else
        {
          // debug
          if (result->log_debug)
          {
            result->debug_message("elim %d from %d rows\n",c,(int)col_rows[c].size()-1);
            for (auto s : col_rows[c])
            {
              if (s<rref_r)
                result->info_message("*");
              else if (s==rr)
                result->info_message("+");
              else
                result->info_message("-");
              print_row_iia(s);
            }
          }
          
          // changed_vars that we need to update the priority for
          changed_var_set.clear();
          for (auto r : col_rows[c])
          {
            changed_var_set.insert(rows[r].cols.begin(),rows[r].cols.end());
          }
            
          if(!rref_elim(rref_r, rr, c, rref_col_order, rref_col_map))
            return false;
          
          if (rref_r>used_row)
            return true;
          
          progress_since_pushback=true;
          
          // update variables
          update_var.clear();
          for (auto c : changed_var_set )
          {
            if (rref_col_map[c]<0 && !col_rows[c].empty() && col_rows[c].back() >= rref_r)
              update_var.push_back(c);
          }
          Q.update( update_var, false );
        }
      }
      // try again with the left-over vars, if any
      if (all_vars.empty() || !progress_since_pushback)
        return true;
    }
    return true;
  }
  
  bool IncrementalIntervalAssignment::rref_step2(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map)
  {
    if (rref_r>used_row)
      return true;
    
    result->debug_message("rref_step2\n");
    // eliminate in order of increasing gcd (decreasing -gcd)
    // find the gcd of all columns, store in a priority queue, and add a new Q entry for the changed columns when we do a row-eliminate
    
    // row simplify once before finding gcd
    for (int rr=rref_r; rr<=used_row; ++rr)
      row_simplify(rr);
    
    // build priority Q
    priority_queue< pair<int,int> >  Q; // gcd,col
    for (int c=0; c<=used_col;++c)
    {
      if (rref_col_map[c]==-1)
      {
        int g = gcd_col(rref_r,c);
        if (g>0)
        {
          // result->info_message("gcd_to_cols.insert <%d, %d>\n", g, c);
          Q.emplace(-g,c);
        }
      }
    }
    
    // process in gcd order, re-inserting if gcd got worse
    while (rref_r<=used_row && !Q.empty())
    {
      auto gc=Q.top();
      Q.pop();
      // result->info_message("popped <%d, %d>, remaining size %d\n", -gc.first, gc.second, (int)Q.size());
      
      // tip top, update val
      int c = gc.second;
      assert(c<=used_col);
      int g=gcd_col(rref_r,c);
      if (g>0)
      {
        if (g>-gc.first)
        {
          // result->info_message("reinsert gcd_to_cols.insert <%d, %d>\n", g, c);
          Q.emplace(-g,c); // reinsert in queue
        }
        else
        {
          // update gc.first to avoid confusion
          if (gc.first!=-g) gc.first=-g;
          
          // swap a row into position
          // best one is the one with the smallest coefficient, and rr>=r
          auto &cr = col_rows[c];
          {
            int best_row=-1;
            int best_row_c=numeric_limits<int>::max();
            for (size_t i=0; i<cr.size(); ++i)
            {
              if (cr[i]>=rref_r)
              {
                if (best_row==-1)
                {
                  best_row=cr[i];
                  best_row_c=abs(rows[cr[i]].get_val(c));
                }
                else
                {
                  const int v=abs(rows[cr[i]].get_val(c));
                  if (v<best_row_c)
                  {
                    best_row=cr[i];
                    best_row_c=v;
                  }
                }
              }
            }
            assert(best_row>-1);
            assert(best_row<(int)rows.size());
            assert(best_row_c<numeric_limits<int>::max());
            row_swap(rref_r,best_row);
          }
          
          // ok, we're working with row r
          // keep multiplying it and subtracting the multiple of other rows until we get down to their gcd
          auto &row=rows[rref_r];
          // result->info_message("working with row %d ", r); row.print_row(result);
          auto &row_vals=rows[rref_r].vals;
          int ci = row.get_col_index(c);
          auto &c_coeff=row_vals[ci];
          // result->info_message("coeff of col %d is %d\n", c, c_coeff);
          
          // find subset of rows to work with
          //   ri points to row r
          int ri=0;
          while (cr[ri]<rref_r)
            ++ri;
          assert( cr[ri]==rref_r );
          
          // for all the remaining rows, until we get c_coeff minimal
          while (abs(c_coeff)>g && ri+1<(int)cr.size())
          {
            // result->info_message("Coeff not minimal, |coeff|>g, |%d|>%d\n", c_coeff, g);
            
            int rr=cr[++ri]; // next row to multiply and subtract
            assert(rr>rref_r);
            int cc_coeff=rows[rr].get_val(c); // coeff of c in that row
            if (abs(c_coeff)>=abs(cc_coeff))
            {
              result->error_message("why am I working with row %d coeff %d, when I could have been using row %d coeff %d?\n", rref_r, c_coeff, rr, cc_coeff);
              print_problem("bad rref");
              return false;
            }
            assert(abs(c_coeff)<abs(cc_coeff));
            int gg = gcd(c_coeff,cc_coeff); // this is as good as we can get using this row
            if (gg==abs(c_coeff))
              continue;
            int  m = abs( c_coeff)/gg;
            int mm = abs(cc_coeff)/gg;
            int n=1;
            while ( (n*m) % mm != 1 )
              ++n;
            int nn = (n*m) / mm;
            // row[r] = n * row[r] - nn * row[rr]
            const int new_m = n * m - nn * mm;
            assert(new_m==1);
            if (c_coeff<0)
              n = -n;
            if (cc_coeff<0)
              nn = -nn;
            row_us_minus_vr(rr, rref_r, n, nn);
            assert(c_coeff>0);
            assert(c_coeff==new_m*gg);
            if(c_coeff!=new_m*gg) // some coding mistake? or weird coefficients?
            {
              result->error_message("IIA: some coding mistake in rref_step2, unexpected value after matrix row algebra.\n");
              return false;
            }
          }
          // result->info_message("final coeff of col %d is %d\n", c, c_coeff);
          // result->info_message("final row r %d is",r); row.print_row(result);
          
          // do the elimination
          if (!rref_elim(rref_r, rref_r, c, rref_col_order, rref_col_map))
            return false;
          row_simplify(rref_r-1);
          
          // print_problem("rref2 progress");
          
        } // Q element work
      } // g is a good value >0
    } // while !Q.empty
    return true;
  }
  
  bool IncrementalIntervalAssignment::rref_stepZ(int &rref_r, vector<int> &rref_col_order, vector<int> &rref_col_map)
  {
    result->debug_message("rref_stepZ\n");
    for (int c=0; c<=used_col; ++c)
    {
      if (rref_r>used_row)
        break;
      
      // column to use = c = lead
      if (rref_col_map[c]!=-1 || col_rows[c].empty())
        continue;
      
      // row to use, with smallest coeff
      int rr=-1;
      for (auto rrr : col_rows[c])
      {
        if (rrr<rref_r)
          continue;
        
        row_simplify(rrr);
        if (rr==-1 || abs(rows[rr].get_val(c))>abs(rows[rrr].get_val(c)))
        {
          rr = rrr;
        }
      }
      if (rr==-1)
        continue;
      
      if (!rref_elim(rref_r, rr, c, rref_col_order, rref_col_map))
        return false;
    }
    
    // Assign order to remaining columns
    for (int c=0; c<=used_col; ++c)
    {
      if (rref_col_map[c]==-1)
      {
        assert(rref_col_map[c]==-1);
        rref_col_map[c]=(int)rref_col_order.size();
        rref_col_order.push_back(c);
      }
    }
    
    return true;
  }
  
  bool IncrementalIntervalAssignment::rref_constraints(vector<int> &rref_col_order, bool map_only_phase)
  {
    // rref data
    // col_order is implicit order of columns, the sequence in which rref eliminated them
    // col_map is -1 if we have not added them to the order, otherwise their index in col_order
    
    // order of columns in rref
    rref_col_order.clear();
    const size_t col_size=used_col+1+10; // 10 is a buffer
    rref_col_order.reserve(col_size);
    
    // map from true column c to the index in which c appears in col_order
    // last entry is a safety, always -1
    vector<int> rref_col_map(col_size,-1);
    
    int rref_r=0; // next row to reduce
    
    // Step 0
    // if any var is in only one col, and has a "1" for a coefficient, we're good immediately!
    if (!rref_step0(rref_r, rref_col_order, rref_col_map))
      return false;
    
    // Step 1
    // Use any var that has a "1" for a coefficient in a row, prefering those with fewer rows
    if (!rref_step1(rref_r, rref_col_order, rref_col_map, map_only_phase))
      return false;
    
    // Step 2
    // Select vars in order of increasing gcd of its coefficients (in the unused rows)
    // if the gcd is 1, then it is possible to combine rows to get a leading coefficient of 1
    if (!rref_step2(rref_r, rref_col_order, rref_col_map))
      return false;
    
    // Step (3) Last Resort
    // Use any var with a non-zero coefficient in a remaining row!
    if (!rref_stepZ(rref_r, rref_col_order, rref_col_map))
      return false;
    
    return true;
  }
   
  bool IncrementalIntervalAssignment::rref_improve(vector<int> &rref_col_order)
  {
    // rref data
    // col_order is implicit order of columns, the sequence in which rref eliminated them
    // col_map is -1 if we have not added them to the order, otherwise their index in col_order
    
    // order of columns in rref
    rref_col_order.clear();
    const size_t col_size=used_col+1+10; // 10 is a buffer
    rref_col_order.reserve(col_size);
    
    // map from true column c to the index in which c appears in col_order
    // last entry is a safety, always -1
    vector<int> rref_col_map(col_size,-1);
    
    int rref_r=0; // next row to reduce
    
    // Step 0
    // pick ones in fewer rows
    if (!rref_step_numrows(rref_r, rref_col_order, rref_col_map))
      return false;
    
    // Step (end) Last Resort, safety, left-overs
    // Use any var with a non-zero coefficient in a remaining row!
    if (!rref_stepZ(rref_r, rref_col_order, rref_col_map))
      return false;
    
    return true;
  }
  
  void MatrixSparseInt::row_union(int r, int s, vector<int> &rsc)
  {
    auto &rc=rows[r].cols;
    auto &sc=rows[s].cols;
    set_union(rc.begin(),rc.end(),
                   sc.begin(),sc.end(),
                   back_inserter(rsc));
  }
  
  void MatrixSparseInt::row_swap(int r, int s)
  {
    if (r==s)
      return;
    
    assert(r>-1);
    assert(s>-1);
    assert((size_t)r<rows.size());
    assert((size_t)s<rows.size());
    
    rows[r].cols.swap(rows[s].cols);
    rows[r].vals.swap(rows[s].vals);
    
    // find colums in r or s
    vector<int>rsc;
    row_union(r,s,rsc);
    // update columns that contain r or s
    for (auto c : rsc)
    {
      for (auto &t : col_rows[c])
      {
        if      (t==r) t=s;
        else if (t==s) t=r;
      }
      // to do: potential speedup, keep track of indices of t and s, and just sort between those indices
      sort(col_rows[c].begin(),col_rows[c].end());
    }
  }
  
  void IncrementalIntervalAssignment::row_swap(int r, int s)
  {
    // optimization
    //   this is a lot of data movement. performance hit?
    //   Consider instead keeping a vector of indices and just updating those.
    if (r==s)
      return;
    
    // swap rows
    MatrixSparseInt::row_swap(r,s);
    
    swap(rhs[r],rhs[s]);
    swap(constraint[r],constraint[s]);
    
    // swap data from IncrementalIntervalAssignment
    if (parentProblem)
    {
      swap(parentRows[r],parentRows[s]);
    }
  }

  void MatrixSparseInt::matrix_allrows_eliminate_column(int r, int c, vector<int> *rhs)
  {
    // eliminate c from other rows, except r
    
    // optimized for speed because this takes most of the time when improving
    // It does this:
    //    auto kill_rows = col_rows[c]; // copy, since matrix_row_elimination changes col_rows[c]
    //    for (auto s : kill_rows)
    //    {
    //      if (s==r)
    //        continue;
    //
    //      matrix_row_eliminate_column(r, s, c); // c is always removed from col_rows everywhere, so could update that last. todo optimization
    //      row_simplify_vals(rows[s].vals);
    //    }
    
    
    if (col_rows[c].size()>1)
    {
      const auto u = rows[r].get_val(c);
      assert(u!=0);
      auto &r_cols=rows[r].cols;
      auto &r_vals=rows[r].vals;
      const auto imax=r_cols.size();
      assert(imax>0);
      
      for (auto s : col_rows[c])
      {
        if (s==r)
          continue;
        
        // matrix_row_eliminate_column(r, s, c);
        
        auto v = rows[s].get_val(c);
        assert(v!=0);
        
        // matrix_row_us_minus_vr(r, s, u, v);
        // matrix_row_us_minus_vr(rows[r], rows[s], s, &col_rows, u, v);
        int min_val = numeric_limits<int>::max()/2;
        auto &s_cols=rows[s].cols;
        auto &s_vals=rows[s].vals;
        {
          size_t i=0,j=0; // index into r, s
          const auto jmax=s_cols.size();
          assert(jmax>0);
          
          // new row we are building
          // speedup by making these permanent working space in the matrix
          new_cols.clear();
          new_vals.clear();
          new_cols.reserve(imax+jmax);
          new_vals.reserve(imax+jmax);
          
          // go down both rows in sync
          while (i<imax || j<jmax)
          {
            // just r[i] is non-zero
            if (i<imax && (j==jmax || r_cols[i]<s_cols[j]))
            {
              new_cols.push_back(   r_cols[i]);
              const auto nv = -v*r_vals[i];
              new_vals.push_back(nv);
              min_val = min(min_val, abs(nv));
              
              // let variable r_cols[i] know that it has a non-zero in row s (and it was zero before)
              auto &crow = col_rows[r_cols[i]];
              // add s to crow, maintaining the sorted order
              crow.push_back(s);
              inplace_merge(crow.begin(), prev(crow.end()), crow.end()); // restore sorted order
              
              ++i;
            }
            // just s[j] is non-zero
            else if (j<jmax && (i==imax || r_cols[i]>s_cols[j]))
            {
              new_cols.push_back(  s_cols[j]);
              const auto nv =    u*s_vals[j];
              new_vals.push_back(nv);
              min_val = min(min_val, abs(nv));
              ++j;
            }
            // both r[i],s[j] are non-zero
            else
            {
              assert(i<imax);
              assert(j<jmax);
              assert(r_cols[i]==s_cols[j]);
              const auto nv = u*s_vals[j]-v*r_vals[i];
              if (nv != 0)
              {
                new_cols.push_back( r_cols[i]);
                new_vals.push_back( nv );
                min_val = min(min_val, abs(nv));
              }
              else
              {
                // we special case c because we know it's just going to be in one row in the end
                if (r_cols[i]!=c)
                {
                  // let variable r_cols[i] know that it no longer has a non-zero in row s
                  auto &crow = col_rows[r_cols[i]];
                  // remove s from crow, maintaining the sort
                  auto next = upper_bound(crow.begin(),crow.end(),s);
                  move(next,crow.end(),prev(next));
                  crow.pop_back();
                }
              }
              ++i;
              ++j;
            }
          }
          // save the new row as rows[s]
          s_cols.swap(new_cols);
          s_vals.swap(new_vals);
        }
        
        if (rhs)
        {
          auto b = u * (*rhs)[s] - v * (*rhs)[r];
          if (b!=0)
          {
            min_val = min(min_val, abs(b));
            if (min_val != 1)
            {
              s_vals.push_back(b);
              row_simplify_vals(s_vals);
              b = s_vals.back();
              s_vals.pop_back();
            }
          }
          (*rhs)[s] = b;
        }
        else if (min_val != 1)
          row_simplify_vals(s_vals);
        // else the row already has a 1, so it can't be simplified, so don't try
        
        assert(rows[s].get_val(c)==0);
      } // for all rows
      
      // update col_rows[c] to be just r
      col_rows[c].clear();
      col_rows[c].push_back(r);
    }
  }
    
  void MatrixSparseInt::matrix_row_us_minus_vr(const RowSparseInt &r_row, RowSparseInt &s_row, int s, vector< vector<int> > *col_rows_loc, int u, int v)
  {
    assert(u!=0);
    assert(v!=0);
    auto &r_cols=r_row.cols;
    auto &s_cols=s_row.cols;
    auto &r_vals=r_row.vals;
    auto &s_vals=s_row.vals;
    size_t i=0,j=0; // index into r, s
    const auto imax=r_cols.size();
    const auto jmax=s_cols.size();
    assert(imax>0);
    assert(jmax>0);
    
    // optimization, use permanent working space in the matrix
    new_cols.clear();
    new_vals.clear();
    new_cols.reserve(imax+jmax);
    new_vals.reserve(imax+jmax);
    
    // divide each of these into separate functions and run again
    
    while (i<imax || j<jmax)
    {
      // just r[i] is non-zero
      if (i<imax && (j==jmax || r_cols[i]<s_cols[j]))
      {
        new_cols.push_back(   r_cols[i]);
        new_vals.push_back(-v*r_vals[i]);
        
        if (col_rows_loc)
        {
          // let variable r_cols[i] know that it has a non-zero in row s (and it was zero before)
          auto &crow = (*col_rows_loc)[r_cols[i]];
          // add s to crow, maintaining the sorted order
          crow.push_back(s);
          inplace_merge(crow.begin(), prev(crow.end()), crow.end()); // restore sorted order
        }
        
        ++i;
      }
      // just s[j] is non-zero
      else if (j<jmax && (i==imax || r_cols[i]>s_cols[j]))
      {
        new_cols.push_back(  s_cols[j]);
        new_vals.push_back(u*s_vals[j]);
        ++j;
      }
      // both r[i],s[j] are non-zero
      else
      {
        assert(i<imax);
        assert(j<jmax);
        assert(r_cols[i]==s_cols[j]);
        const int new_val = u*s_vals[j]-v*r_vals[i];
        if (new_val != 0)
        {
          new_cols.push_back( r_cols[i]);
          new_vals.push_back( new_val);
        }
        else
        {
          if (col_rows_loc)
          {
            // let variable r_cols[i] know that it no longer has a non-zero in row s
            auto &crow = (*col_rows_loc)[r_cols[i]];
            // result->info_message("column %d's rows, before removing row %d\n",r_cols[i],s); print_vec(result, crow);
            // remove s from crow, maintaining the sort
            auto next = upper_bound(crow.begin(),crow.end(),s);
            move(next,crow.end(),prev(next));
            crow.pop_back();
            // result->info_message("column %d's rows, after removing row %d\n",r_cols[i],s); print_vec(result, crow);
          }
        }
        ++i;
        ++j;
      }
    }
    s_cols.swap(new_cols);
    s_vals.swap(new_vals);
  }
  
  void IncrementalIntervalAssignment::multiply_rhs(int r, int s, int u, int v)
  {
    if (u==1 && v==0)
      return;
    
    // rhs update
    auto  rb = rhs[r];
    auto &sb = rhs[s];
    if (rb!=0 || sb!=0)
    {
      sb = u*sb-v*rb;
    }
  }
  
  void IncrementalIntervalAssignment::row_us_minus_vr(int r, int s, int u, int v)
  {
    this->matrix_row_us_minus_vr(r,s,u,v);
    multiply_rhs(r,s,u,v);
  }
  
  void MatrixSparseInt::transpose(MatrixSparseInt &T,  size_t num_rows, size_t num_cols) const
  {
    T.col_rows.clear();
    T.rows.clear();
    if (num_rows==0)
      num_rows=rows.size();
    if (num_cols==0)
      num_cols=col_rows.size();
    T.col_rows.resize(num_rows);
    T.rows.resize(num_cols);
    // For all entries, flip row and column and assign to T
    for (int r = 0; r < (int)rows.size(); ++r)
    {
      auto &row = rows[r];
      for (size_t i = 0; i < row.cols.size(); ++i)
      {
        auto c = row.cols[i];
        auto v = row.vals[i];
        auto &Trow = T.rows[c];
        Trow.cols.push_back(r);
        Trow.vals.push_back(v);
      }
    }
    // Sorting is not needed, since we pushed in order of increasing r
    //   do it anyway for safety, its fast compared to everything else we do
    for (auto &Trow : T.rows)
    {
      Trow.sort();
    }
    // needed
    T.fill_in_cols_from_rows();
  }
  
  void MatrixSparseInt::identity(MatrixSparseInt &I, size_t sz)
  {
    I.rows.clear();
    I.col_rows.clear();
    I.rows.resize(sz);
    I.col_rows.resize(sz);
    for (size_t i=0; i < sz; ++i)
    {
      I.rows[i].cols.push_back((int)i);
      I.rows[i].vals.push_back(1);
      I.col_rows[i].push_back((int)i);
    }
  }

  bool IncrementalIntervalAssignment::HNF(MatrixSparseInt &B, MatrixSparseInt &U, vector<int> &hnf_col_order,
                                          vector<int> &g)
  {
    if (result->log_debug)
    {
      result->debug_message("IIA starting HNF\n");
      print_problem("Trying HNF on this");
    }
    
    if (result->log_debug)
    {
      result->debug_message("IIA starting HNF\n");
      print_problem("Trying HNF on this");
    }
    
    // Reduce to only the non-zero rows, as RREF may have removed some redundant constraints from the end
    int r(used_row);
    while (r>-1 && rows[r].cols.empty())
    {
      --r;
    }
    // no constraints, set any solution we want!
    if (r==-1)
    {
      // col_solution is already assigned
      assert( col_solution == g );
      return true;
    }
    const int this_rows = r+1; // m
    
    // Reduce to non-zero columns
    int c(used_col);
    while(c>-1 && col_rows[c].empty())
    {
      // col_solution[c] is arbitrary, as it's in no constraints
      // col_solution[c]=g[c]; // already assigned,
      --c;
    }
    const int this_cols = used_col+1; // n
    
    return MatrixSparseInt::HNF(this_rows, this_cols, B, U, hnf_col_order, g, use_HNF_goals);
  }

  void MatrixSparseInt::nonzero(int &r, int &c) const
  {
    r = (int) rows.size(); --r;
    while (r>-1 && rows[r].cols.empty())
    {
      --r;
    }

    c = (int) col_rows.size(); --c;
    while(c>-1 && col_rows[c].empty())
    {
      --c;
    }
  }

  bool check_inverse(const MatrixSparseInt &U, const MatrixSparseInt &Uinv)
  {
    // check U * Uinv == I
    MatrixSparseInt UinvT(Uinv.result);
    Uinv.transpose(UinvT);
    vector<int> UinvColr(U.rows.size(),0), er(U.rows.size(),0);
    for (int r=0; r<U.rows.size(); ++r)
    {
      // col r of Uinv == row r of UinvT
      UinvT.rows[r].sparse_to_full(UinvColr);
      U.multiply(UinvColr, er);
      // er should be elementary vector with 1 in position r
      if (er[r]!=1)
        return false;
      for (int rr=0; rr<er.size(); ++rr)
      {
        if (rr!=r && er[rr]!=0)
        {
          return false;
        }
      }
    }
    return true;
  }

  bool MatrixSparseInt::HNF(int this_rows, int this_cols, MatrixSparseInt &B, MatrixSparseInt &U, vector<int> &hnf_col_order, vector<int> &g, bool use_HNF_goals) const
  {

    // Description of HNF Hermite Normal Form and why it is what we want to do
    // http://sites.math.rutgers.edu/~sk1233/courses/ANT-F14/lec3.pdf
    // http://sites.math.rutgers.edu/~sk1233/courses/ANT-F14/lec4.pdf
    
    // Form column-form of HNF, the lower triangular matrix
    // Implemented via row operations on the transpose instead of column operations,
    //   because row-ops are what is implemented, and, more importantly, row-ops are what we've optimized the matrix data structure for.
    // Write all comments as if we are doing column operations, to make it easier to follow, but "IMP" describes what that translates to in our implementation
    
    // Copy this to B.  IMP copy this^T to B
    // B = This.  IMP B = This^T
    // U = I
    
    const bool debug_HNF(false);
    
    transpose(B, this_rows, this_cols);
    identity(U, this_cols);
    // debug
    MatrixSparseInt UinvT(U);
    identity(UinvT, this_cols);
    
    hnf_col_order.resize(this_cols);
    iota(hnf_col_order.begin(),hnf_col_order.end(),0);
    
    // debug
    if (result->log_debug)
    {
      B.print_matrix("HNF B initial");
      U.print_matrix("HNF U initial");
      UinvT.print_matrix("HNF UinvT initial");
      MatrixSparseInt Uinv(UinvT);
      UinvT.transpose(Uinv);
      assert( check_inverse(U,Uinv) );
      result->info_message("HNF column order initial");
      print_vec(result, hnf_col_order);
      // result->info_message("HNF U^{-1}g initial"); print_vec(result, g);
    }
    
    // zzyk to do, if needed
    // If columns are shorter than the diagonal, we should do row swaps to make all pivots on the diagonal
    // This might not happen because we did RREF and removed redundant constraints already?
    // i.e. if we get this
    //   1                        1
    //   1       => row swap      0 1
    //   0 1                      1
    
    // IMP hnf_r is row we're working on, and hnf_c the pivot for that row
    int hnf_c(0);
    for (int hnf_r = 0; hnf_r < (int)this_rows; ++hnf_r)
    {
      result->debug_message("hnf_r %d, hnf_c %d\n",hnf_r, hnf_c);
      
      // For each column, ensure top-row coefficient is positive, by multiplying columns by -1 as needed
      // IMP for each row, ensure first coefficient is positive, by multiplying rows by -1 as needed
      // to do optimization - only do this loop once, for the first variable. Otherwise update the sign of the row as needed after the multiplies
      for (int br : B.col_rows[hnf_r])
      {
        if (br<hnf_c)
          continue;
        auto &row = B.rows[br];
        // assert(row.get_val(hnf_r) != 0); could be zero, if already reduced, no more columns to reduce, etc.
        
        if (row.get_val(hnf_c)<0)
        {
          for (auto &vv : B.rows[br].vals)
          {
            vv = -vv;
          }
          // same operation on U
          for (auto &vv : U.rows[br].vals)
          {
            vv = -vv;
          }
          // negate UinvT row == negate Uinv column
          for (auto &vv : UinvT.rows[br].vals)
          {
            vv = -vv;
          }
          // same operation on g
          // g[br] = -g[br];
        }
        assert(B.rows[br].get_val(hnf_c) >= 0);
      }
      if (/* DISABLES CODE */ (debug_HNF))
      {
        B.print_matrix("HNF B after column hnf_r lower diagonal entries made positive");
        U.print_matrix("HNF U after column hnf_r lower diagonal entries made positive");
        UinvT.print_matrix("HNF UinvT after column hnf_r lower diagonal entries made positive");
        MatrixSparseInt Uinv(UinvT);
        UinvT.transpose(Uinv);
        assert( check_inverse(U,Uinv) );
        // result->info_message("HNF U^{-1}g after column hnf_r lower diagonal entries made positive"); print_vec(result, g);
      }
      
      
      // Get one non-zero in row hnf_r (ignore columns to the left) by adding multiples of columns together
      //   IMP. get one non-zero in col by adding multiples of other rows, ignoring higher rows
      
      // Dumb simple way to do this is to iteratively loop over all rows, reducing coefficients by the smallest coefficient we have so far,
      //   until no more reduction is possible.
      // There is some faster way using just one pass and using the gcd, but this works for us, and the slowdown is minor because our coefficients tend to be 1.
      size_t extra_nonzeros=0;
      int small_r=-1; // row with the smallest (or only) positive coefficient
      do
      {
        // find row with smallest coeff (they are all positive)
        small_r=-1;
        extra_nonzeros=0;
        int small_v=numeric_limits<int>::max(); // coefficient of the row with the smallest coefficient
        
        for (int br : B.col_rows[hnf_r])
        {
          if (br<hnf_c)
            continue;
          auto &row = B.rows[br];
          const auto v = row.get_val(hnf_c);
          assert(v != 0);
          
          if (v<small_v)
          {
            small_v = v;
            small_r = br;
          }
        }
        result->debug_message("row %d has smallest coefficient %d for col %d\n",small_r,small_v,hnf_c);
        
        // subtract small_r from all other rows
        // copy col_rows since matrix_row_us_minus_vr changes it
        auto cr = B.col_rows[hnf_c]; // copy
        for (int br : cr)
        {
          if (br<hnf_c)
            continue;
          if (br==small_r)
          {
            continue;
          }
          int v = B.rows[br].get_val(hnf_r);
          assert(v!=0);
          int m = v / small_v; // integer
          // row[br] := row[br] - m * row[small_r]
          B.matrix_row_us_minus_vr(small_r, br, 1, m);
          U.matrix_row_us_minus_vr(small_r, br, 1, m);
          // Uinv  row br += m row small_r
          // UinvT row small_r += m row br
          UinvT.matrix_row_us_minus_vr(br, small_r, 1, -m);
          // g := U^{-1}g, where for Uinv use +m
          // g[br] += m*g[small_r]; // ?? applied in wrong order, need to do UinvTT.multiply(g,e)
          if (v % small_v == 0)
          {
            assert( B.rows[br].get_val(hnf_r) == 0);
          }
          else
          {
            assert( B.rows[br].get_val(hnf_r) != 0);
            ++extra_nonzeros;
          }
          // zzyk to do, robustness
          // test for overflow
          // if we get overflow, the integers get too big, we need to implement the extra step in
          //   Section 4.2 in http://sites.math.rutgers.edu/~sk1233/courses/ANT-F14/lec3.pdf
          // B^hat = [ Z kI_{m by m} ] where Z are m independent columns of B, and k is det(Z)
          // Do HNF on B^hat instead of B, since they have the same integer span. Just throw away the extra m variables at the end.
          // Further, prevent overflow as follows:
          //   If entry (i,j) of B exceeds threshold k, we pick the ith column from the kIm×m sub-matrix and subtract a large enough
          //   multiple of this column from the jth column of the matrix so that the (i,j) entry becomes less than k
          //   (This amounts to subtracting a large enough integer multiple of k from the (i,j) entry without disturbing the other
          //    entries of the jth column of the matrix).
          // Can this be simplified? We already did RREF, so getting the determinant is just the product of the RREF diagonal.
          // Do we need to explicitly add kI_{m by m} to B^hat? Can we just subract k if an entry is too large (or add k if it is too small)?
          // I'm confused, lec3 implies the extra columns of B can be discarded, but don't we need those to reduce the diagonal coefficients?
        }
        
      } while (extra_nonzeros);
      
      if (/* DISABLES CODE */ (debug_HNF))
      {
        result->debug_message("done reducing hnf_r %d hnf_c %d\n",hnf_r,hnf_c);
        B.print_matrix("HNF B after reducing hnf_r");
        U.print_matrix("HNF U after reducing hnf_r");
        UinvT.print_matrix("HNF UinvT after reducing hnf_r");
        MatrixSparseInt Uinv(UinvT);
        UinvT.transpose(Uinv);
        assert( check_inverse(U,Uinv) );
        // result->info_message("HNF U^{-1}g after reducing hnf_r"); print_vec(result, g);
      }
      
      if (small_r==-1)
      {
        // move this col to the end and keep working? are we sure we are done?
        result->debug_message("HNF done early, there was a redundant constraint or too few constraints, no nonzeros in col %d\n",hnf_r);
        continue;
      }
      
      // swap small_r into hnf_r
      if (small_r != hnf_r)
      {
        B.row_swap(small_r,hnf_r);
        U.row_swap(small_r,hnf_r);
        
        // swap cols of Uinv == swap rows of UinvT == swap g entries
        UinvT.row_swap(small_r,hnf_r);
        // swap(g[small_r],g[hnf_r]);
        
        swap(hnf_col_order[small_r], hnf_col_order[hnf_r]);
      }
      
      if (/* DISABLES CODE */ (debug_HNF))
      {
        result->debug_message("swapped rows %d and %d after reducing hnf_r %d hnf_c %d\n", small_r, hnf_r, hnf_r, hnf_c);
        B.print_matrix("HNF B after reducing hnf_r");
        U.print_matrix("HNF U after reducing hnf_r");
        UinvT.print_matrix("HNF UinvT after reducing hnf_r");
        MatrixSparseInt Uinv(UinvT);
        UinvT.transpose(Uinv);
        assert( check_inverse(U,Uinv) );
        // result->info_message("HNF U^{-1}g after reducing hnf_r"); print_vec(result, g);
      }
      
      ++hnf_c;
    }
    if (result->log_debug)
    {
      result->debug_message("Done reducing column (row) heights\n");
      B.print_matrix("HNF B reduced");
      U.print_matrix("HNF U reduced");
      UinvT.print_matrix("HNF UinvT reduced");
      MatrixSparseInt Uinv(UinvT);
      UinvT.transpose(Uinv);
      assert( check_inverse(U,Uinv) );
      // result->info_message("HNF U^{-1}g reduced"); print_vec(result, g);
    }
    
    // Ensure off-diagonal entries are non-negative and smaller than the diagonal for each row
    // IMP ensure off-diagaonl entries are non-negative and smaller than the diagonal for each column
    // IMP nothing to do for first row, start at r=1
    for (int r = 1; r < (int) B.rows.size(); ++r)
    {
      auto &row = B.rows[r];
      // IMP diagonal entry is the first one of his row
      if (row.vals.empty())
        continue;
      const auto d = row.vals[0];
      const auto c = row.cols[0];
      assert(d>0);
      // IMP now ensure smaller row values are non-negative and smaller
      auto cr = B.col_rows[c]; // copy 'cause it changes
      for (int rr : cr)
      {
        if (rr==r)
          break; // IMP all subsequent rows are too big, should always happen on the last col_rows[c]
        
        const auto v = B.rows[rr].get_val(c);
        int m=0;
        // IMP extra minus to negate the "minus_vr" operation later
        if (v<0)
        {
          m = - (abs(v) / d + (abs(v) % d ? 1 : 0)); // - ceil(|v|/d)
        }
        else if (v>=d)
        {
          m = v / d;  // floor (v/d)
        }
        else
          continue;
        
        // row[rr] := row[rr] - m * row[r]
        B.matrix_row_us_minus_vr(r, rr, 1, m);
        U.matrix_row_us_minus_vr(r, rr, 1, m);
        // for Uinv and g, opposite sign
        // Uinv  row rr += m row r
        // UinvT row  r += m row rr
        UinvT.matrix_row_us_minus_vr(rr, r, 1, -m);
        // g[rr] += m*g[r]; // ?? wrong order

        if (/* DISABLES CODE */ (0 && debug_HNF))
        {
          result->debug_message("Done ensuring smaller off-diagonals for diagonal of row %d and other row %d\n", r, rr);
          B.print_matrix("HNF B off-diagonals");
          U.print_matrix("HNF U off-diagonals");
          UinvT.print_matrix("HNF UinvT off-diagonals");
          MatrixSparseInt Uinv(UinvT);
          UinvT.transpose(Uinv);
          assert( check_inverse(U,Uinv) );
          // result->info_message("HNF U^{-1}g reduced"); print_vec(result, g);
        }
      }
    }
    
    if (result->log_debug)
    {
      result->debug_message("Done ensuring smaller off-diagonals\n");
      B.print_matrix("HNF B off-diagonals");
      U.print_matrix("HNF U off-diagonals");
      UinvT.print_matrix("HNF UinvT off-diagonals");
      MatrixSparseInt Uinv(UinvT);
      UinvT.transpose(Uinv);
      assert( check_inverse(U,Uinv) );
      // result->info_message("HNF U^{-1}g reduced"); print_vec(result, g);
    }
    
    // transpose B and U at end, since we did row operations
    MatrixSparseInt BT(result), UT(result);
    B.transpose(BT);
    U.transpose(UT);
    swap(BT,B);
    swap(UT,U);
    
    // UinvT is now just Uinv
    // apply to g
    vector<int> gg(g.size(),1.);
    if (use_HNF_goals)
    {
      // use c = Uinv g for independent variables.
      UinvT.multiply(g, gg);
    }
    else
    {
      // use c = 1 for independent variables.

    }
    g.swap(gg);
    

    
    // debug
    if (result->log_debug)
    {
      B.print_matrix("HNF B untransposed final");
      U.print_matrix("HNF U untransposed final");
      result->info_message("HNF column order final");
      print_vec(result, hnf_col_order);
      UinvT.print_matrix("HNF Uinv untransposed final");
      assert( check_inverse(U,UinvT) );
      result->info_message("HNF U^{-1}g untransposed final"); print_vec(result, g);
    }
    
    // debug
    // sanity check
    //   B lower triangular
    //   positive diagonals
    //   row entries positive and smaller than the diagonal entry
    //   which non-conforming rows
    if (result->log_debug)
    {
      for (int r = 0; r < (int)B.rows.size(); ++r)
      {
        auto &row = B.rows[r];
        auto d = row.get_val(r);
        if (d==0)
          result->warning_message("HNF B row %d diagonal is zero, skipping other sanity checks.\n", d);
        else
        {
          if (d<0)
          {
            result->error_message("HNF B row %d diagonal is negative, %d<0.\n", r, d);
          }
          for (size_t i = 0; i < row.cols.size(); ++i)
          {
            int c = row.cols[i];
            int v = row.vals[i];
            if (c>r)
            {
              result->error_message("HNF B not lower triangular, B[%d,%d] = %d\n",r,c,v);
            }
            if (v<0)
            {
              result->error_message("HNF B lower triangular entry is negative, B[%d,%d] = %d\n",r,c,v);
            }
            if (v>d)
            {
              result->error_message("HNF B off-diagonal entry is not smaller than the diagaonal entry B[%d,%d] = %d !< %d = B[%d,%d]\n",r,c,v,d,r,r);
            }
          }
        }
      }
    }
    assert( check_inverse(U,UinvT) );
    
    return true;
  }

  bool IncrementalIntervalAssignment::HNF_satisfy_constraints(MatrixSparseInt &B, MatrixSparseInt &U, vector<int> &hnf_col_order,
                                                              const vector<int> &g)
  {
    bool OK = HNF_satisfy_rhs(B, U, hnf_col_order,
                              rhs, col_solution, used_col+1, g);
    if (!OK)
    {
      result->debug_message("HNF, no solution");
      return false;
    }
    
    if (result->log_debug)
    {
      print_solution("HNF constraint solution");
    }
    
    // verify constraints
    bool satisfied=verify_constraints_satisfied(result->log_info || result->log_debug);
    if (!satisfied && result->log_debug)
    {
      print_solution("bad solution does not satisfy constraints");
      result->debug_message("The HNF that resulted in this bad solutions is:\n");
      B.print_matrix("B matrix");
      U.print_matrix("U matrix");
      print_problem("from original problem");
    }
    
    return satisfied;    
  }

  
  bool IncrementalIntervalAssignment::HNF_satisfy_rhs(MatrixSparseInt &B, MatrixSparseInt &U, vector<int> &hnf_col_order,
                                                      const vector<int> &rhs, vector<int> &col_solution, int col_size,
                                                      const vector<int> &g)
  {
    // Let e = [B^{−1}b, 0, 0, . . . , 0] where e is a vector of length n with first m entries B^{−1}b
    // If e is integral, then output Uc, else output "no solution"
    
    // Solve e = B^{-1}b by back substitution, i.e. find c : Be = b
    // if there is no solution, return false
    
    // multiply e by U,  solution = Ue
    // verify solution, if good, return true, else return false;
    
    const bool debug_HNF(false);

    vector<int> e(col_size,0); // rhs == b
                               // cc are the unassigned variables, coming from the redundant constraints
    
    // backsub
    for (int r=0; r< (int)B.rows.size(); ++r)
    {
      const auto &row = B.rows[r];
      auto rs = row_sum(row.cols,row.vals,e,rhs[r]);
      auto coeff = row.get_val(r);
      if (coeff==0)
      {
        e[r]=1;
      }
      else if (rs % coeff)
      {
        if (B.result->log_warning)
          B.result->warning_message("HNF no integer solution because row_sum[%d] %% coeff[%d] != 0:  %d %% %d = %d\n", r, r, rs, coeff, rs % coeff);
        return false;
      }
      else
      {
        e[r] = -rs / coeff;
      }
    }
    // zzyk-hnf
    // set e[c] to be U^{-1}g, where g is the vector of goals
    
    // The non-diagonal columns of U are nullspace vectors we can add in any linear combination to get the values we want
    // we can assign these remaining variables to anything we want
    //   we'd like to pick something that get's us closer to the bounds or goals
    //   this starts to look like the optimization step, maybe just add the columns of U to the nullspace...
    for (size_t c=B.rows.size(); c< e.size(); ++c)
    {
      // e[c] can be anything and not affect satisfying the constraints
      // improvement since IMR paper: assign to goal
      //   Note whatever transformation of variables HNF did to e was already applied to the goals as well.
      e[c] = g[c];
    }
    
    if (/* DISABLES CODE */ (debug_HNF))
    {
      U.print_matrix("Final U");
      U.get_result()->info_message("e (y): ");
      print_vec( U.result, e);
    }
    
    // rhs = Uc
    // U is empty if there were no constraints. Treat it as identity
    if (U.rows.empty())
    {
      assert(e.size()==col_size);
      assert(col_solution.size()==col_size);
      col_solution = e;
    }
    else
    {
      U.multiply(e,col_solution);
    }

    if (/* DISABLES CODE */ (debug_HNF))
    {
      U.get_result()->info_message("col_solution (x): ");
      print_vec( U.result, col_solution);
    }

    return true;
  }
  
  bool IncrementalIntervalAssignment::verify_constraints_satisfied(bool print_me)
  {
    // verify constraints
    bool satisfied=true;
    for (int r=0; r<=used_row; ++r)
    {
      auto rs = row_sum(r);
      if (constraint[r] == EQ)
      {
        if (rs!=0)
        {
          if (print_me)
          {
            print_row_iia(r);
            print_solution_row(r);
            result->error_message("Interval matching equality constraint %d not satisfied, row_sum = %d != 0.\n",r,rs);
          }
          satisfied=false;
        }
      }
      else if (constraint[r] == LE)
      {
        if (rs>0)
        {
          if (print_me)
          {
            print_row_iia(r);
            print_solution_row(r);
            result->error_message("Interval matching less-than-or-equal constraint %d not satisfied, row_sum = %d > 0.\n",r,rs);
          }
          satisfied=false;
        }
      }
      else if (constraint[r] == GE)
      {
        if (rs<0)
        {
          if (print_me)
          {
            print_row_iia(r);
            print_solution_row(r);
            result->error_message("Interval matching greater-than-or-equal constraint %d not satisfied, row_sum = %d < 0.\n",r,rs);
          }
          satisfied=false;
        }
      }
    }
    return satisfied;
  }
  
  bool MatrixSparseInt::multiply(const vector<int> &x, vector<int> &y) const
  {
    if (col_rows.size()!=x.size())
    {
      result->error_message("Matrix-vector size mismatch. Can't multipy matrix with %lu columns with column vector of different length %lu.\n", (unsigned long) col_rows.size(), (unsigned long) x.size());
      return false;
    }
    // y.resize(rows.size());
    if (y.size()<rows.size())
    {
      result->error_message("Matrix-vector size mismatch. Matrix multiplication would produce a vector of length %lu, but result vector is shorter, length %lu.\n", (unsigned long) col_rows.size(), (unsigned long) y.size());
      return false;
    }
    for (size_t r=0; r < rows.size(); ++r)
    {
      auto &row = rows[r];
      auto &yr = y[r];
      yr=0;
      for (size_t i=0; i < row.cols.size(); ++i)
      {
        auto c = row.cols[i];
        auto v = row.vals[i];
        yr += v * x[c];
      }
    }
    
    return true;
  }
  
  void IncrementalIntervalAssignment::copy_bounds_to_sub( IncrementalIntervalAssignment *sub_problem ) const
  {
    // set lower/upper bounds, etc., in sub-problem.
    for (int c = 0; c < sub_problem->lastCopiedCol; c++ )
    {
      const auto pc = sub_problem->parentCols[c];
      sub_problem->col_lower_bounds[c] = col_lower_bounds[ pc ];
      sub_problem->col_upper_bounds[c] = col_upper_bounds[ pc ];
      sub_problem->col_type[c]         = col_type[pc];
      sub_problem->goals[c]            = goals[pc];
      sub_problem->col_solution[c]     = col_solution[pc];
    }
  }
  
  void IncrementalIntervalAssignment::copy_submatrix(const vector <int> &sub_rows, const vector <int> &sub_cols,
                                                     IncrementalIntervalAssignment *target ) const
  {
    // for each parent column, which sub_matrix column is it?
    map<int,int> col_map;
    for ( int c = 0; c < sub_cols.size(); c++ )
    {
      col_map[ sub_cols[c] ] = c;
    }

    // for ( r : rows )
    //   target_row = row_map[r]
    //   copy this row r into target row target_row
    //     covert this column c into target column col_map[c]
    target->used_row=((int)sub_rows.size())-1;
    target->used_col=((int)sub_cols.size())-1;
    
    // t is row of sub_problem
    for (int t = 0; t < (int) sub_rows.size(); ++t)
    {
      int r = sub_rows[t];
      // for all (col,val) pairs of the row
      for (size_t i = 0; i < rows[r].cols.size(); ++i)
      {
        const auto c = rows[r].cols[i];
        auto dd = col_map.find(c);
        if (dd != col_map.end()) // if column is part of the subproblem...
        {
          int d = dd->second; // col of sub_problem
          target->rows[t].cols.push_back(d);
          target->rows[t].vals.push_back(rows[r].vals[i]);
        }
      }
      assert((size_t)t<target->rhs.size());
      assert((size_t)t<target->constraint.size());

      target->rhs[t]=rhs[r];
      target->constraint[t]=constraint[r];
    }
    
    // d is col of sub_problem
    for (int d = 0; d < (int) sub_cols.size(); ++d)
    {
      int c = sub_cols[d];

      target->col_solution    [d]=col_solution    [c];
      target->col_lower_bounds[d]=col_lower_bounds[c];
      target->col_upper_bounds[d]=col_upper_bounds[c];
    }
    
    target->fill_in_cols_from_rows();
  }
  
  int IncrementalIntervalAssignment::next_available_row( ConstraintType constraint_type)
  {
    if (++used_row>=number_of_rows)
      add_more_rows();
    constraint[used_row]=constraint_type;
    return used_row;
  }
  
  int IncrementalIntervalAssignment::next_dummy()
  {
    if (++used_col >= number_of_cols)
    {
      add_more_columns();
    }
    // change if from an INT_VAR to a DUMMY_VAR
    col_type[used_col]=DUMMY_VAR;
    goals[used_col]=0;
    col_lower_bounds[used_col] = numeric_limits<int>::lowest();
    col_upper_bounds[used_col] = numeric_limits<int>::max();
    return used_col;
  }

  int IncrementalIntervalAssignment::next_intvar()
  {
    if (++used_col >= number_of_cols)
    {
      add_more_columns();
    }
    return used_col;
  }

  void IncrementalIntervalAssignment::count_used_rows_cols()
  {
    int max_r = -1, max_c = -1;
    for (int r = 0; r < number_of_rows; ++r)
    {
      auto &row_cols = rows[r].cols;
      if (!row_cols.empty())
      {
        max_r = r;
      }
      for (auto c : row_cols)
      {
        max_c = max(c, max_c);
      }
    }
    used_row = max_r;
    used_col = max_c;
  }
  
  
  void IncrementalIntervalAssignment::resize_rows()
  {
    // use number_of_xxx not how many are used so far
    rows.resize(number_of_rows);
    constraint.resize(number_of_rows,EQ);
    rhs.resize(number_of_rows,0);
  }
  
  void IncrementalIntervalAssignment::resize_cols()
  {
    // default to standard for an interval variable INT_VAR range [1,inf], goal 1.
    col_rows.resize(number_of_cols);
    col_lower_bounds.resize(number_of_cols,1);
    col_upper_bounds.resize(number_of_cols,numeric_limits<int>::max());
    col_solution.resize(number_of_cols,0);
    col_type.resize(number_of_cols,INT_VAR);
    goals.resize(number_of_cols,1.);
    
    // skip tied variables
    // we only resize_cols on parent problems, and only allocate tied_variables on sub_problems
  }
  
  void IncrementalIntervalAssignment::freeze_problem_size()
  {
    resize_rows();
    resize_cols();
  }
  
  void IncrementalIntervalAssignment::add_more_rows()
  {
    // const auto old_rows = number_of_rows;
    number_of_rows = (number_of_rows * 3) / 2 + 8;
    resize_rows();
  }
  
  void IncrementalIntervalAssignment::add_more_columns()
  {
    number_of_cols = (number_of_cols * 3) / 2 + 8;
    resize_cols();
  }
  
  int IncrementalIntervalAssignment::new_row(int num_rows)
  {
    // dynamically update in a way that we can continue solving where we left off
    if (num_rows<1)
      return -1;
    
    int return_row = used_row+1;
    used_row += num_rows;
    
    // add more rows if needed
    if (used_row>=(int)rows.size())
    {
      add_more_rows();
    }
    return return_row;
  }
  
  int IncrementalIntervalAssignment::new_col(int num_cols)
  {
    // dynamically update in a way that we can continue solving where we left off
    if (num_cols<1)
      return -1;
    
    int return_col = used_col+1;
    used_col += num_cols;
    
    // add more rows if needed
    if (used_col>=(int)col_rows.size())
    {
      add_more_columns();
    }
    return return_col;
  }
  
  void IncrementalIntervalAssignment::reserve_cols( int num_cols )
  {
    number_of_cols=max(number_of_cols,num_cols);
  }
  void IncrementalIntervalAssignment::reserve_rows( int num_rows )
  {
    number_of_rows=max(number_of_rows,num_rows);
  }
  
  // from https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
  template <typename T, typename Compare>
  vector<size_t> sort_permutation(const vector<T>& vec,
                                            Compare compare)
  {
    vector<size_t> p(vec.size());
    iota(p.begin(), p.end(), 0);
    sort(p.begin(), p.end(),
              [&](size_t i, size_t j){ return compare(vec[i], vec[j]); });
    return p;
  }
  template <typename T>
  vector<T> apply_permutation(const vector<T>& vec,
                                   const vector<size_t>& p)
  {
    vector<T> sorted_vec(vec.size());
    transform(p.begin(), p.end(), sorted_vec.begin(),
                   [&](size_t i){ return vec[i]; });
    return sorted_vec;
  }
  
  void RowSparseInt::sort()
  {
    auto p = sort_permutation(cols,
                              [](int const& a, int const& b)
                              { return a<b; });
    auto cols_vec = apply_permutation(cols,p);
    cols.swap(cols_vec);
    auto vals_vec = apply_permutation(vals,p);
    vals.swap(vals_vec);
  }
  
  void IncrementalIntervalAssignment::sort_rows(int row_start)
  {
    for (int r=row_start; r<=used_row;++r)
    {
      rows[r].sort();
    }
  }
  
  void IncrementalIntervalAssignment::fill_in_cols_from_rows()
  {
    // clear all the columns first, if not already
    //   needed if we previously solved this, and then u-submap constraints are added, and now we are resolving
    for (auto &cr : col_rows)
    {
      cr.clear();
    }
    
    // fill in col_row data
    for (int r=0; r<=used_row; ++r)
    {
      auto &cols = rows[r].cols;
      for (auto c : cols)
      {
        assert(c<(int)col_rows.size());
        col_rows[c].push_back(r);
        used_col=max(used_col,c);
      }
    }
  }
  
  void MatrixSparseInt::fill_in_cols_from_rows()
  {
    // clear current columns
    for (auto &cr : col_rows)
    {
      // assert(cr.empty()); // it is not empty when we are refilling from cull_nullspace_tied_variable_rows
      cr.clear();
    }
    
    // count how many columns there are
    int max_col=-1;
    for (size_t r=0; r<rows.size(); ++r)
    {
      auto &cols = rows[r].cols;
      for (auto c : cols)
      {
        max_col = max( max_col, int(c) );
      }
    }
    if (max_col >= (int)col_rows.size())
      col_rows.resize(max_col+1);
    
    // fill in data
    for (size_t r=0; r<rows.size(); ++r)
    {
      auto &cols = rows[r].cols;
      for (auto c : cols)
      {
        col_rows[c].push_back((int)r);
      }
    }
  }
  
  int RowSparseInt::get_val( int col ) const
  {
    auto found_col = lower_bound(cols.begin(),cols.end(),col);
    if (found_col!=cols.end() && *found_col==col)
      return vals[ found_col-cols.begin() ];
    return 0;
  }
  int IncrementalIntervalAssignment::get_M_unsorted  ( int row, int col ) const
  {
    auto &r = rows[row];
    for (size_t i = 0; i < r.cols.size(); ++i)
    {
      if (r.cols[i]==col)
        return r.vals[i];
    }
    return 0;
  }
  
  void IncrementalIntervalAssignment::set_M( int row, int col, int val )
  {
    // for efficiency
    //   rows columns may get out of order
    //   column rows are not set
    // caller may need to sort rows and fill in cols
    auto &cols = rows[row].cols;
    auto &vals = rows[row].vals;
    // overwrite old value, if any
    for (size_t i=0; i< cols.size(); ++i)
    {
      if (cols[i]==col)
      {
        if (val==0)
        {
          if (i+1<cols.size())
          {
            swap(cols[i],cols.back());
            swap(vals[i],vals.back());
          }
          cols.pop_back();
          vals.pop_back();
          // rows[row].sort();
        }
        else
          vals[i]=val;
        val = 0; // flag to not add the value
        break;
      }
    }
    if (val!=0)
    {
      // push new non-zero
      cols.push_back(col);
      vals.push_back(val);
      // rows[row].sort();
    }
    
    hasBeenSolved=false;
    if (row <= solved_used_row)
    {
      solved_used_row = -1;
      solved_used_col = -1;
    }
  }
  
  void IncrementalIntervalAssignment::set_M( int row,  vector<int> &cols, vector<int> &vals )
  {
    rows[row].cols.swap(cols);
    rows[row].vals.swap(vals);
    // not sorted yet
    // columns not filled in yet

    hasBeenSolved=false;
    if (row <= solved_used_row)
    {
      solved_used_row = -1;
      solved_used_col = -1;
    }
  }

  // verify Ax=b
  bool MatrixSparseInt::verify_Ax_b(const vector<int> &x, const vector<int> &b) const
  {
    // bb = Ax
    vector<int> bb(rows.size(),0);
    bool OK = multiply(x, bb);
    if (!OK)
    {
      result->error_message("verify_Ax=b fail, computation failure.");
      return false;
    }
    // b.size == bb.size
    if (bb.size() != b.size())
    {
      result->error_message("verify_Ax=b fail, Ax not the same size as b.");
      return false;
    }
    // b == bb
    for (int i=0; i<bb.size(); i++)
    {
      if (b[i]!=bb[i])
      {
        result->error_message("verify_Ax=b fail, b[%d]==%d but Ax[%d]==%d.",i,b[i],i,bb[i]);
        return false;
      }
    }
    return true;
  }

  int MatrixSparseInt::lcm_entries() const
  {
    int L=1;
    for (auto &row : rows)
    {
      for (auto v : row.vals)
      {
        assert(v!=0);
        L=lcm(L,v);
      }
    }
    return L;
  }

  void remove_duplicate_rows(const MatrixSparseInt &M, int MaxMrow, MatrixSparseInt &N)
  {
    set<int> delete_row;
    for (int r=0; r<MaxMrow; ++r)
    {
      // does row r of M exist in N?
      auto &mrow = M.rows[r];
      assert(!mrow.cols.empty());
      int lead_col = mrow.cols[0];
      // for rows in N with same leading col
      for (int s : N.col_rows[lead_col])
      {
        if (row_equal(mrow, N.rows[s]))
        {
          delete_row.insert(s);
        }
      }
    }
    if (!delete_row.empty())
    {
      // overwrite rows of N
      auto si = delete_row.begin();
      int ss = 0; // assign row s to be ss
      for (int s = 0; s<N.rows.size(); ++s)
      {
        if (si!=delete_row.end() && s==*si)
        {
          si++;
          // don't increment ss
        }
        else
        {
          if (s!=ss)
          {
            N.rows[ss]=N.rows[s];
          }
          ss++; // ss increments with s
        }
      }
      N.rows.resize(ss);
      N.fill_in_cols_from_rows();
    }
    
  }

  int IncrementalIntervalAssignment::solve_sub(bool do_improve, bool map_only_phase)
  {
    // solve a subproblem
    //   get row echelon form
    //   find nullspace from row echelon form
    //   get initial feasible solution
    //     and satisfy variable bounds
    //   optimize solution quality
    //     increment by (linear combinations of) vectors from the nullspace
    
    // debug
    result->debug_message( "\n---- solving solve_sub ----\n");
    CpuTimer timer;
    double time_rref(0.), time_constrain(0.), time_bound(0.), time_opt(0.), time_tune(0.);
    // debug
    
    if (find_tied_variables())
    {
      if (result->log_debug)
      {
        print_problem("before removing tied variables\n");
      }
      
      // remove tied variables from the matrix and replace them with the tied-to variables
      replace_tied_variables();
      cull_tied_variables(Nullspace);
      
      // generate the upper/lower bounds/goals for the tied-to variables
      if (!generate_tied_data())
      {
        result->constraints_satisfied=false;
        result->bounds_satisfied=false;
        return false;
      }
      
      // adjust the initial solution for the tied-to variables
      if (should_adjust_solution_tied_variables)
        adjust_solution_tied_variables();
      
      if (result->log_debug)
      {
        print_problem("after removing tied variables and adjusting solution\n");
      }
      
      // to debug just the tied variables
      // return false;
    }
    
    // initial values, something close to feasible,
    /* feasible= */
    assign_dummies_feasible(); // if we already have a feasible solution, this does nothing
    
    if (result->log_debug)
    {
      print_problem("solve_sub before RREF");
      // print_problem_summary("solve_sub before RREF");
    }
    
    // it saves time and gets better output to compute a different rref for solving the constraints vs. improving the bounds and solution
    // if (do_restore) is false, then we just reuse the constraint rref for the bounds/improve step
    const bool do_restore=true;
    
    // satisfy constraints
    vector<int> rref_col_order;
    vector<int>cols_dep, cols_ind, rows_dep;
    {
      // save original matrix A,b
      MatrixSparseInt *this_matrix = this;
      EqualsB *this_B = this;
      MatrixSparseInt Asave(*this_matrix);
      EqualsB Bsave(*this_B);
      
      // row echelon for constraints
      bool rref_OK = rref_constraints(rref_col_order,map_only_phase);
      
      // debug
      if (result->log_debug)
      {
        print_problem("solve_sub after RREF constraints",rref_col_order);
        // print_problem_summary("solve_sub after RREF constraints");
      }
      if (result->log_debug_time)
      {
        time_rref = timer.cpu_secs();
        result->info_message("rref time %f\n",time_rref);
      }
      // debug
      
      if (!rref_OK)
      {
        result->error=true;
        result->constraints_satisfied=false;
        result->bounds_satisfied=false;
        return false;
      }
      
      // satisfy constraints using direct back-substitution using the rref, with HNF as backup
      categorize_vars(rref_col_order, cols_dep, cols_ind, rows_dep);
      const bool is_feasible = satisfy_constraints(cols_dep, cols_ind, rows_dep);
      
      // debug
      // print_problem("solve_sub initial_feasible_solution ");
      if (result->log_debug)
      {
        print_solution("solve_sub initial constraints solution ");
      }
      if (result->log_debug_time)
      {
        time_constrain = timer.cpu_secs();
        result->info_message("feasible constraints solution time %f\n",time_constrain);
      }
      // debug
      
      if (!is_feasible)
      {
        result->constraints_satisfied=false;
        result->bounds_satisfied=false;
        return false;
      }
      
      // restore original matrix, before rref
      if (do_restore)
      {
        swap(*this_matrix,Asave);
        swap(*this_B,Bsave);
        
        if (result->log_debug)
        {
          print_problem("solve_sub after restoring to before RREF state");
          // print_problem_summary("solve_sub after restoring to before RREF state");
        }
      }
      
    }
    
    // research, skip bounds and improve
    if (satisfy_constraints_only)
    {
      if (has_tied_variables)
      {
        assign_tied_variables();
      }
      if (result->constraints_satisfied)
      {
        hasBeenSolved = true;
      }
      return result->constraints_satisfied;
    }
    // zzyk to do, for large problems, attempt to satisfy bounds and improve solution just using Nullspace, without any rref
    
    // setup nullspace for improving bounds and quality
    // M is nullspace
    // rows 0..MaxMrow span the nullspace, and beyond that are extra rows that were useful candidate increment vectors
    MatrixSparseInt M(result);
    int MaxMrow=0;
    // MaxMrow is last row of nullspace M before we augmented it
    // to prepend Nullspace to M, use the following two lines
    //   MatrixSparseInt &M = Nullspace;
    //   int MaxMrow=(int) M.rows.size(); // last row of nullspace M before we augmented it

    {
      // save original matrix A,b
      MatrixSparseInt *this_matrix = this;
      EqualsB *this_B = this;
      MatrixSparseInt Asave(*this_matrix);
      EqualsB Bsave(*this_B);
      
      if (do_restore)
      {
        rref_col_order.clear();
        rref_improve(rref_col_order);
        
        // debug
        if (result->log_debug)
        {
          print_problem("solve_sub after RREF improve", rref_col_order);
          // print_problem_summary("solve_sub after RREF improve");
        }
        //debug
        
        cols_dep.clear(); cols_ind.clear(); rows_dep.clear();
        categorize_vars(rref_col_order, cols_dep, cols_ind, rows_dep);        
      }
      
      // M nullspace
      create_nullspace(cols_dep, cols_ind, rows_dep, M, MaxMrow);
      cull_nullspace_tied_variable_rows(M, MaxMrow);
      
      // debug
      if (result->log_debug) // || use_map_nullspace)
      {
        M.print_matrix("M nullspace");
        // M.print_matrix_summary("M nullspace ");
        if (use_map_nullspace)
        {
          Nullspace.print_matrix("N nullspace");
          // Nullspace.print_matrix_summary("N nullspace ");
        }
      }
      
      if (do_restore)
      {
        // restore original matrix, before rref
        swap(*this_matrix,Asave);
        swap(*this_B,Bsave);
        // print_problem("solve_sub after restore");
      }
    }
    
    // MN = concatenate M+Nullspace
    //   if use_map_nullspace is false, then Nullspace is an empty matrix here and MN==M
    MatrixSparseInt MN(M,MaxMrow,Nullspace,true);
    if (result->log_debug && use_map_nullspace)
    {
      MN.print_matrix("MN nullspace");
    }
    
    bool in_bounds = false;
    
    // normal method
    in_bounds = satisfy_bounds(M, MN);
    
    // satisfy_bounds with tiny subspace
    //   for all the variables that are out of bounds
    //   find tiny subspace
    //   do the above on it instead
    if (!in_bounds)
    {
      MatrixSparseInt Mtiny(M);
      bool more_nullspace_rows = false;
      for ( int c = 0; c <= used_col; ++c )
      {
        if (bounds_unsatisfied(c))
        {
          if (result->log_debug)
          {
            result->debug_message("bounds not satisfied for %d\n",c);
            print_problem("stuck problem");
          }
          
          if (tiny_subspace(c, Mtiny, MN))
          {
            more_nullspace_rows = true;
          }
        }
      }
      if (more_nullspace_rows)
      {
        in_bounds = satisfy_bounds(Mtiny, MN);
      }
    }
    
    // debug
    if (result->log_debug)
    {
      print_solution("solve_sub constraints + bounds solution ");
    }
    if (result->log_debug_time)
    {
      time_bound = timer.cpu_secs();
      result->info_message("feasible bounds solution time %f\n",time_bound);
    }
    // debug
    
    if (!in_bounds)
    {
      result->bounds_satisfied=false;
      return false;
    }
    
    hasBeenSolved=true; // feasible solution
    
    if (do_improve)
    {
      // redo rref and nullspace if helpful
      //   variables with non-1 coefficients, such as sum-even variables with a coeff of 2, were chosen as independent vars above.
      //   however, for optimization, we want them to be *dependent* variables
      
      if (/* DISABLES CODE */ (0) && use_map_nullspace)
      {
        result->info_message("do_improve phase:\n");
        M.print_matrix("M nullspace");
        // M.print_matrix_summary("M nullspace ");
        Nullspace.print_matrix("N nullspace");
        // Nullspace.print_matrix_summary("N nullspace ");
      }
      
      improve_solution(M, MN);
      
      // debug
      if (result->log_debug_time)
      {
        time_opt = timer.cpu_secs();
        result->info_message("optimize quality time %f\n",time_opt);
      }
      // debug
      
      fine_tune(MN);
      
      // debug
      if (result->log_debug_time)
      {
        time_tune = timer.cpu_secs();
        result->info_message("fine tune quality time %f\n",time_tune);
      }
      // debug
    }
    
    if (has_tied_variables)
    {
      assign_tied_variables();
      
      // debug
      if (result->log_debug)
      {
        print_solution("tied variable solution");
        // print_solution_summary("tied variable solution");
      }
    }
    
    if (result->log_debug_time)
    {
      double time_assign_tied = timer.cpu_secs();
      result->info_message("assign tied variable time %f\n",time_tune);
      result->info_message("total solve_sub time %f\n", time_rref + time_constrain + time_bound + time_opt + time_tune + time_assign_tied);
    }
    result->debug_message( "\n---- done solve_sub ----\n\n");

    // save M in nullspace to pass to parent
    Nullspace = MatrixSparseInt(M,MaxMrow);
    // re-insert tied-to variables
    if (!tied_variables.empty())
    {
      for (int c=0; (size_t) c<tied_variables.size(); ++c)
      {
        if (tied_variables[c].size()>1)
        {
          for ( int r : Nullspace.col_rows[c] )
          {
            const int v = Nullspace.rows[r].get_val(c);
            for ( int c2 : tied_variables[c] )
            {
              if (c2==c) continue;
              Nullspace.rows[r].cols.push_back(c2);
              Nullspace.rows[r].vals.push_back(v);
              Nullspace.col_rows[c2].push_back(r);
            }
            Nullspace.rows[r].sort();
          }
        }
      }
      // Nullspace.fill_in_cols_from_rows(); // not needed
    }
    
    
    return true;
  }

  void IncrementalIntervalAssignment::identify_col_type()
  {
    for (int c = 0; c <= used_col; ++c)
    {
      // INT_VARS have a goal, and usually a positive lower bound
      // DUMMY_VARS have no goal, and arbitrary bounds
      if (goals[c]==0)
      {
        col_type[c]=DUMMY_VAR;
        // EVEN_VARS have a coefficient of 2 in a single row
        if (col_rows[c].size()==1)
        {
          auto v = get_M(col_rows[c][0],c);
          if (abs(v) % 2 == 0)
          {
            col_type[c]=EVEN_VAR;
          }
        }
      }
      else
      {
        col_type[c]=INT_VAR;
      }
    }
  }

  void IncrementalIntervalAssignment::force_satisfied_with_slack(int r)
  {
    // already satisfied?
    if (row_sum(r)==0)
      return;

    auto &row = rows[r];
    
    // Is there already a slack variable we can use?
    //   E.g. the user may have already included one.
    for (size_t i = 0; i < row.cols.size(); ++i)
    {
      const auto c = row.cols[i];
      if (col_rows[c].size()==1) // require var is in single row so that we know we aren't making some other row unsatisfied
      {
        auto old_solution = col_solution[c];
        col_solution[c]=0;
        auto sum = row_sum(r);
        auto v = row.vals[i];
        if (sum % v==0)
        {
          col_solution[c] = -sum / v;
          assert( row_sum(r) == 0);
          return; // success
        }
        else
        {
          // restore solution, keep looking
          col_solution[c]=old_solution;
        }
      }
    }
    
    // Need to make a new slack variable
    int slack_var = next_dummy();
    row.cols.push_back(slack_var);
    row.vals.push_back(-1);
    col_rows[slack_var].push_back(r);
    col_type[slack_var] = DUMMY_VAR;
    goals[slack_var]=0.;
    
    // slack_var will start out of bounds, then we optimize to put it into bounds
    col_lower_bounds[slack_var]=0;
    col_upper_bounds[slack_var]=0;
    
    // assign slack variable a value
    int slack = row_sum(r);
    col_solution[slack_var]=slack;
    assert( row_sum(r) == 0);
  }

  bool IncrementalIntervalAssignment::convert_inequalities()
  {
    auto old_used_col = used_col;
    
    // convert inequalities into equalities with dummy variables
    int num_even=0, num_ineq=0, num_eq=0, num_new=0;
    const bool is_re_solving = solved_used_row>=0;
    for (int r = max(0,solved_used_row+1); r <= used_row; ++r)
    {
      const bool row_is_new = (is_re_solving && r>solved_used_row);
      if (row_is_new)
        ++num_new;
      
      const auto con=constraint[r];
      if (con == EVEN)
      {
        // convert constraint
        auto c = next_dummy();
        col_type[c]=EVEN_VAR;
        goals[c]=0.;
        rows[r].cols.push_back(c);
        rows[r].vals.push_back(-2);
        col_rows[c].push_back(r);
        col_lower_bounds[c]=2; // sum is at least 4
        // no upper bound
        
        constraint[r]=EQ;
        
        ++num_even;
        
        if (row_is_new)
        {
          int slack = row_sum(r);
          col_solution[c]=slack/2;
          if (row_sum(r)!=0)
          {
            // add *second* dummy to make it satisfied
            force_satisfied_with_slack(r);
          }
        }
      }
      else if (con==LE || con==GE)
      {
        // convert constraint
        auto c = next_dummy();
        col_type[c] = DUMMY_VAR;
        goals[c]=0.;
        rows[r].cols.push_back(c);
        const int c_coeff = (con==LE ? 1 : -1);
        rows[r].vals.push_back(c_coeff);
        col_rows[c].push_back(r);
        col_lower_bounds[c]=0;
        // no upper bound
        
        constraint[r]=EQ;
        
        ++num_ineq;
        
        if (row_is_new)
        {
          int slack = row_sum(r);
          col_solution[c] = -c_coeff * slack;
          assert( row_sum(r) == 0);
        }
        
      }
      else if (row_is_new && con==EQ)
      {
        ++num_eq;
        // for new eq constraints, introduce a slack variable so we start with an initial feasible solution.
        force_satisfied_with_slack(r);
      }
    }
    
    // debug
    if (result->log_debug)
    {
      if (num_new)
      {
        result->debug_message("%d new constraints.\n",num_new);
      }
      if (num_even)
      {
        result->debug_message("%d sum-even constraints converted to equalities with dummy variables.\n",num_even);
      }
      if (num_ineq)
      {
        result->debug_message("%d inequality constraints converted to equalities with dummy variables.\n",num_ineq);
      }
      if (num_eq)
      {
        result->debug_message("%d new equality constraints satisfied.\n",num_ineq);
      }
      const auto num_slack = used_col-old_used_col;
      if (num_slack)
      {
        result->debug_message("%d dummy (slack) variables introduced.\n",num_slack);
      }
      if (num_even || num_ineq || num_eq)
      {
        print_problem("full problem initial after inequality conversion");
      }
      else
      {
        result->debug_message("No inequalities or even constraints, no new equalities to satisfy, nothing to convert.\n");
      }
    }
    // debug
    
    return (num_even || num_ineq || num_eq);
  }
  
  bool IncrementalIntervalAssignment::solve(bool do_improve, bool first_time)
  {
    // debug
    result->info_message("Running incremental interval assignment.\n");
    // what options were turned on
    result->info_message("use HNF_always, HNF_goals, map_nullspace, best_improvement, better_rref_pivots, fixed_RatioR, "
                         "smart_rounding, Ratio_for_improvement, flip_blocker_signs"
                         " = %d%d%d%d%d%d%d%d%d\n",
                         use_HNF_always,
                         use_HNF_goals,
                         use_map_nullspace,
                         use_best_improvement,
                         use_better_rref_step1_pivots_in_sumeven_phase,
                         use_fixed_RatioR_straddle,
                         use_smart_adjust_for_fractions,
                         use_Ratio_for_improvement_queue,
                         use_Nullspace_sumeven_to_flip_sign_of_blocker);
    
    CpuTimer total_timer;
    double setup_time=0, map_time(0.);
    if (result->log_debug)
    {
      print_problem("full problem initial before setup");
      // print_problem_summary("full problem initial before setup");
    }
    // debug
    
    if (solved_used_row<0)
      first_time=true;
    else if (first_time)
    {
      solved_used_row = -1; // force it to ignore any prior solution
      solved_used_col = -1;
    }
    
    count_used_rows_cols();
    
    // sort rows, fill in columns
    const auto new_row_start=max(0,solved_used_row+1);
    sort_rows(new_row_start);
    fill_in_cols_from_rows(); // could be more efficient about this the second time we are solving

    identify_col_type();
    
    convert_inequalities();
    
    int new_row_min = -1, new_row_max = -1;
    if (first_time)
    {
      // derived-class-specific initialization
      assign_vars_goals(true);
      should_adjust_solution_tied_variables = true;
    }
    else // !first_time
    {
      new_row_min = solved_used_row+1;
      new_row_max = used_row;
      
      // don't init solution, don't adjust tied variables
      should_adjust_solution_tied_variables = false;
    }

    if (result->log_debug_time)
    {
      setup_time = total_timer.cpu_secs();
      result->info_message("IIA setup time = %f\n", setup_time );
    }
    
    if (infeasible_constraints())
    {
      result->info_message("problem is infeasible\n");
      return false;
    }
    
    // satisfied unless something failes later
    result->constraints_satisfied=true;
    result->bounds_satisfied=true;
    
    // solve mapping constraints, do these separately first for efficiency
    result->debug_message("Solving IIA mapping subproblems.\n");
    
    bool success = solve_phase(true, do_improve, new_row_min, new_row_max);
    if (0)
    {
      print_solution_summary("solution after mapping phase");
    }
    
    if (result->log_debug_time)
    {
      map_time = total_timer.cpu_secs();
      result->info_message("IIA mapping-phase solver time = %f\n", map_time );
    }
    
    // if the mapping problems were feasible
    if ( success )
    {
      result->debug_message("Solving IIA sum-even subproblems.\n");
      // set should_adjust_solution_tied_variables to false
      //   if it was true, then we already did the adjustments we should have when solving the mapping phase
      //   if we adjust again, it may "undo" the mapping solution and cause problems.
      should_adjust_solution_tied_variables = false;
      success = solve_phase(false, do_improve, new_row_min, new_row_max);
    }
    if (success)
    {
      hasBeenSolved=true;
      solved_used_row = used_row;
      solved_used_col = used_col;
    }
    else
    {
      hasBeenSolved=false;
    }
    result->solved = success;

    // debug & diagnositcs
    if (result->log_debug_time)
    {
      const double even_time = total_timer.cpu_secs();
      result->info_message("IIA sum-even-phase solver time = %f\n", even_time );
      const double total_time = setup_time+map_time+even_time;
      result->info_message("Total IIA solve time = %f\n", total_time);
    }
    if (result->log_debug)
    {
      print_problem("full problem solution");
      print_solution("full problem solution");
      // print_problem_summary("full problem solution");
      // print_solution_summary("full problem solution");
    }
    // debug
    
    return success;
  }
  
  int RowSparseInt::dot(const RowSparseInt &b)
  {
    int d(0), i(0),j(0);
    while (i<cols.size() && j<b.cols.size())
    {
      int r =   cols[i];
      int s = b.cols[j];
      // for cols in both rows
      if (r < s)
      {
        ++i;
      }
      else if (r > s)
        ++j;
      else
      {
        d += vals[i] * b.vals[j];
        ++i;
        ++j;
      }
    }
    return d;
  }

  // improvement since IMR paper
  void IncrementalIntervalAssignment::transfer_parent_nullspace_to_subs(vector<IncrementalIntervalAssignment*> &sub_problems)
  {
    if (result->log_debug)
    {
      print_problem("parent");
    }
    
    // make parent to sub column map, using sub to parent map
    const size_t mapsize = used_col+1;
    vector<IncrementalIntervalAssignment*> col_subproblem(mapsize,nullptr);
    parentCols.assign(mapsize,-1);
    for ( auto *sub_problem : sub_problems)
    {
      if (result->log_debug)
      {
        sub_problem->print_problem("sub_problem");
      }
      
      for (int sc=0; sc < sub_problem->parentCols.size(); ++sc )
      {
        int pc = sub_problem->parentCols[sc];
        if (pc>=0)
        {
          parentCols[ pc ] = sc;
          col_subproblem[ pc ] = sub_problem;
        }
      }
    }
    
    // gather relevant nullspace from this to sub problems
    for (auto &prow : Nullspace.rows)
    {
      if (result->log_debug)
      {
        result->info_message("Converting parent nullspace ");
        prow.print_row(result);
      }

      if (prow.cols.empty())
        continue;
      
      // determine which subproblem, and verify *all* the columns are in that subproblem
      IncrementalIntervalAssignment *sub_problem = col_subproblem[ prow.cols[0] ];
      for (int pc : prow.cols)
      {
        if ( col_subproblem[pc] != sub_problem )
        {
          sub_problem = nullptr;
          break;
        }
      }
      if (sub_problem==nullptr)
        continue;
      
      // convert from parent columns to sub_problem columns
      RowSparseInt srow; // sub nullspace row
      srow.vals = prow.vals;
      for (auto pc : prow.cols)
      {
        int sc = parentCols[pc];
        assert(sc>=0);
        srow.cols.push_back(sc);
      }
      
      // iterate over the sub_problem rows that have a variable in common with the nullspace row
      set<int> subrows;
      for (auto c : srow.cols)
      {
        subrows.insert( sub_problem->col_rows[c].begin(), sub_problem->col_rows[c].end() );
      }
      
      // verify row is in the nullspace of the subproblem. If not, figure out the dummy variables that need to be added to it to make it in the nullspace.
      
      // a spanning set of nullspace vectors includes one entry for each row that has a non-trivial contribution.
      // nullspace row = k_row * row + k_dummy * dummy_col
      struct nullrowinfo
      {
        int k_row, k_dummy, dummy_col, c_dummy;
        nullrowinfo() : k_row(0), k_dummy(0), dummy_col(-1), c_dummy(1) {}
        nullrowinfo(int k_row_in, int k_dummy_in, int dummy_col_in, int c_dummy_in) : k_row(k_row_in), k_dummy(k_dummy_in), dummy_col(dummy_col_in), c_dummy(c_dummy_in) {}
      };
      map < int, vector<nullrowinfo> > nullrows;
      bool ok(true);
      
      for (auto r : subrows) // sub constraint row
      {
        auto &row = sub_problem->rows[r];
        
        // calculate row-sum, plus gather any extra varibles in the subproblem row but not in the nullspace vector
        vector<int> extra_cols;
        int row_sum(0), i(0),j(0);
        while (i<row.cols.size() && j<srow.cols.size())
        {
          int rc =  row.cols[i];
          int sc = srow.cols[j];
          // for cols in both rows
          if (rc < sc)
          {
            extra_cols.push_back(rc);
            ++i;
          }
          else if (rc > sc)
            ++j;
          else
          {
            row_sum += row.vals[i] * srow.vals[j];
            ++i;
            ++j;
          }
        }
        while(i<row.cols.size())
        {
          extra_cols.push_back(row.cols[i]);
          ++i;
        }
        
        
        // row_sum==0 means its already in the nullspace of that row. Leave nullrowinfo empty.
        if (row_sum!=0)
        {
          auto &v = nullrows[r]; //creates
                                 // for any extra col, if it is in just the one row, then we're good, we can make a nullspace row for the sub
          for (auto dummy_col : extra_cols)
          {
            if (sub_problem->col_rows[dummy_col].size()==1)
            {
              // assert(col_rows[extra_cols[0]][0]) should be this row
              int coeff = row.get_val(dummy_col);
              int p = abs(lcm(coeff,row_sum));
              
              // ensure k_row > 0 for simpicity
              if (row_sum < 0) p = -p;
              int k_dummy =   -p / coeff;
              int k_row   =    p / row_sum;
              assert(k_row>0);
              // nullspace row = k_row * row + k_dummy * dummy_col
              v.emplace_back(nullrowinfo(k_row,k_dummy,dummy_col,coeff));
            }
          }
          if (v.empty())
          {
            // failed, we can't make a nullspace row this way
            ok = false;
            break;
          }
        }
      }
      
      // now build the nullspace vectors and verify
      if (ok)
      {
        // use the "smallest" dummy variable for every row to make a nullspace vector
        // there is probably a more sophisticated way to pick the one that has a small lcm, but try greedy for now.
        // also calculate k_row, the overall multiple of the nullspace row to use
        vector <int> ii;
        int k_row=1;
        for (auto &n : nullrows)
        {
          int min_i=0, min_kd=1, min_kk=k_row;
          bool min_sign_agrees=false;
          VarType min_var_type=UNKNOWN_VAR; // {INT_VAR, EVEN_VAR, DUMMY_VAR, UNKNOWN_VAR};
          for (int i=0; i<n.second.size();++i)
          {
            const int ki = n.second[i].k_row;
            const int kk = abs(lcm(k_row,ki));
            const int kd = n.second[i].k_dummy;
            const VarType ktype = col_type[ki]; // or should we be looking at the subproblem?

            // selection criteria:
            //   type sum-even var (this is the opposite of prefering small coefficient for multiplying row)
            //   type dummy var
            //   k_dummy is positve, dummy variable will increase with additions to the nullspace vector
            //   small coefficient for multiplying row
            //   small coefficient for multiplying k, i.e. kk
            
            // first time through, set min to the first entry
            bool pick_i  = (i==0);
            bool decided = (i==0);
            
            // prefer sum-even
            if (!decided)
            {
              if (ktype==EVEN_VAR && min_var_type!=EVEN_VAR)
              {
                pick_i=true;
                decided=true;
              }
              else if (ktype!=EVEN_VAR && min_var_type==EVEN_VAR)
              {
                pick_i=false;
                decided=true;
              }
            }
            // prefer dummy
            if (!decided)
            {
              if (ktype==DUMMY_VAR && min_var_type!=DUMMY_VAR)
              {
                pick_i=true;
                decided=true;
              }
              else if (ktype!=DUMMY_VAR && min_var_type==DUMMY_VAR)
              {
                pick_i=false;
                decided=true;
              }
            }
            // prefer coefficient sign
            if (!decided)
            {
              // k_row was always set to be positive
              // k_dummy, want it positive, too, as most dummy vars have a lower bound but no upper bound
              bool sign_agrees = (kd>0);
              if (sign_agrees && !min_sign_agrees)
              {
                pick_i = true;
                decided = true;
              }
              else if (!sign_agrees && min_sign_agrees)
              {
                pick_i = false;
                decided = true;
              }
            }
            // prefer small coefficient for row, small coeff for k
            if (!decided)
            {
              if (kk < min_kk || (kk==min_kk && abs(kd)<abs(min_kd)))
              {
                pick_i = true;
              }
              decided = true;
            }
            assert(decided);
            if (pick_i)
            {
              min_i =i;
              min_kk = kk;
              min_kd = kd;
              min_var_type = ktype;
            }
          }
          k_row=min_kk;
          ii.push_back(min_i);
        }
        
        // generate the subproblem nullspace row
        auto nrow = srow;
        nrow.multiply(k_row);
        
        // gather the dummy vars
        {
          int i=0;
          for (auto &n : nullrows)
          {
            auto &nd = n.second[ii[i]];
            nrow.cols.push_back( nd.dummy_col );
            const int k = k_row / nd.k_row;
            nrow.vals.push_back( k * nd.k_dummy );
            ++i;
          }
        }
        nrow.sort();
        row_simplify_vals(nrow.vals);
        
        if (result->log_debug)
        {
          result->info_message("Subproblem nullspace ");
          nrow.print_row(result);
        }
        
        // debug, double-check that the rowsum is now zero
        if (1)
        {
          for (auto r : subrows) // sub constraint row
          {
            auto &row = sub_problem->rows[r];
            auto row_sum = row.dot(nrow);
            
            if (result->log_debug)
            {
              result->debug_message("row %d row_sum=%d, ",r,row_sum);
              row.print_row(result);
            }
            assert(row_sum==0);
          }
        }
        
        sub_problem->Nullspace.push_row(nrow);
        
        // For each of the rows that had multiple isolated vars, generate a nullspace row for the pairs of one of those vars and the one we chose for above
        //   but only if both of those vars are dummy vars. We catch the INT_VARs compbinations already through the parent nullspace
        {
          int i=0;
          for (auto &n : nullrows)
          {
            if (n.second.size()>1)
            {
              // jj is one of the nullspace indices
              const int jj = ii[i];
              if (col_type[jj] == INT_VAR) continue;
              
              for (int j=0; j<n.second.size(); ++j)
              {
                if (j==jj) continue;
                if (col_type[j] == INT_VAR) continue;

                // make nullspace using n.second[j] and n.second[jj]
                auto cj  = n.second[ j].c_dummy;
                auto cjj = n.second[jj].c_dummy;

                auto k = lcm(cj,cjj);
                auto nj  =  k / cj;
                auto njj = -k / cjj;

                auto dj  = n.second[ j].dummy_col;
                auto djj = n.second[jj].dummy_col;

                RowSparseInt nn({dj,djj}, {nj,njj});
                nn.sort();
                row_simplify_vals(nn.vals);

                if (result->log_debug)
                {
                  result->info_message("Subproblem nullspace isolated pair ");
                  nn.print_row(result);
                }

                sub_problem->Nullspace.push_row(nn);
              }
            }
            ++i;
          }
        }
      } // ok
    }
    // gather relevant nullspace from this to sub problems
    
    // debug subproblem nullspaces
    if (result->log_debug)
    {
      for ( auto *sub_problem : sub_problems)
      {
        //sub_problem->print_problem("sub_problem");
        get_result()->info_message("sub_problem nullspace from parent");
        sub_problem->Nullspace.print_matrix();
      }
    }
  }


  
  bool IncrementalIntervalAssignment::solve_phase(bool map_only_phase, bool do_improve, int row_min, int row_max)
  {
    if (result->log_debug)
    {
      if (map_only_phase)
        result->debug_message("\nSolving interval matching 'map' phase.\n");
      else
        result->debug_message("\nSolving interval matching 'even' phase.\n");
      print_problem("\nEntire Interval Problem");
    }
    
    // hack for test_problem_augmented_nullspace_demo zzyk
    if (0 && !map_only_phase)
    {
      col_solution[0]=4;
      col_solution[1]=4;
      col_solution[2]=4;
      col_solution[3]=5;
    }

    // subdivide problems into independent components
    vector<IncrementalIntervalAssignment*> sub_problems;
    vector<int> orphan_cols;
    subdivide_problem( sub_problems, orphan_cols, !map_only_phase, row_min, row_max  );
    const auto num_subs = sub_problems.size();

    // copy and fix nullspaces from parent to sub
    if (!map_only_phase)
    {
      transfer_parent_nullspace_to_subs(sub_problems);
    }

    // timing
    double solve_time=0.;
    CpuTimer solve_timer;
    
    bool success = true;
    for ( size_t i=0; i < num_subs; ++i)
    {
      auto *sub_problem = sub_problems[i];
      
      if (result->log_debug)
      {
        result->debug_message("Subproblem %d of %d:\n", i+1, num_subs);
        sub_problem->print_problem_summary( "subproblem " + to_string(i+1) + " of " + to_string(num_subs) );
      }
      
      // In the second phase, for variables that have not been assigned anything, assign goals.
      //   This happens for a single paving face, where the variables were not in any mapping subproblem and so are uninitialized
      if (!map_only_phase)
      {
        sub_problem->assign_vars_goals(false);
      }
      
      if (!sub_problem->solve_sub(do_improve,map_only_phase))
      {
        success = false;
        result->warning_message("Subproblem %d of %d was infeasible.\n", i+1, num_subs );
      }
      
      // timing
      if (result->log_debug_time)
      {
        const auto solve_sub_time=solve_timer.cpu_secs();
        solve_time+=solve_sub_time;
        result->info_message("Subproblem %d of %d size %d X %d time = %f (%f total)\n\n",
                             i+1, num_subs, sub_problem->used_row+1, sub_problem->used_col+1, solve_sub_time, solve_time);
      }
      else if (result->log_debug)
      {
        result->info_message("Subproblem %d of %d size %d X %d done\n\n",
                             i+1, num_subs, sub_problem->used_row+1, sub_problem->used_col+1);
      }

    }

    // gather sub-problem soutions back into this problem, copy solution to permanent storage
    if ( success )
      gather_solutions( sub_problems, use_map_nullspace && map_only_phase, orphan_cols );
    delete_subproblems( sub_problems );
    
    return success;
  }
  
  IncrementalIntervalAssignment::~IncrementalIntervalAssignment()
  {
  }
  
  void IncrementalIntervalAssignment::print_row_iia(size_t r) const
  {
    if (parentProblem)
    {
      int pr = parentRows[r];
      result->info_message("  %d (p%d) row: ", (int) r, pr);
    }
    else
    {
      result->info_message("  %d row: ", (int) r);
    }
    auto &row=rows[r];
    auto &cols=row.cols;
    auto &vals=row.vals;
    for (size_t i=0;i<cols.size();++i)
    {
      auto v = vals[i];
      auto c = cols[i];
      result->info_message( "%+dx%d ", v, c );
    }
    
    switch (constraint[r])
    {
      case LE:
        result->info_message(" <= ");
        break;
      case GE:
        result->info_message(" >= ");
        break;
      case EQ:
        result->info_message(" = ");
        break;
      case EVEN:
        result->info_message(" = 2k (even) + ");
        break;
      default:
        result->info_message(" ?unknown_relation ");
    }
    result->info_message(" %d\n", rhs[r]);
  }
  
  void IncrementalIntervalAssignment::print_col_iia(size_t c) const
  {
    // col indices
    result->info_message("  x%d col", (int) c);
    result->info_message(" rows:");
    print_vec(result, col_rows[c], false);
    
    // values
    // expensive
    result->info_message(" vals:");
    print_vec(result, column_coeffs((int)c));
  }
  
  size_t MatrixSparseInt::num_nonzeros() const
  {
    size_t non_zeros=0;
    for (auto &row : rows)
      non_zeros+=row.cols.size();
    return non_zeros;
  }
  
  void IncrementalIntervalAssignment::print_problem_summary(string prefix) const
  {
    if (prefix.empty())
      result->info_message("\n");
    else
      result->info_message("%s ",prefix.c_str());
    result->info_message("IIA problem summary\n");
    result->info_message("  %d rows (0..%d used), %d cols (0..%d used), %lu non-zeros",
                         number_of_rows, used_row,
                         number_of_cols, used_col, (unsigned long) num_nonzeros());
    result->info_message(" %s\n", hasBeenSolved ? "SOLVED" : "unsolved");
  }
  
  void IncrementalIntervalAssignment::print_problem(string prefix) const
  {
    if (prefix.empty())
      result->info_message("\n");
    else
      result->info_message("%s ",prefix.c_str());
    result->info_message("IIA problem\n");
    result->info_message("  number rows = %d (0 to %d used)\n"
                         "  number cols = %d (0 to %d used)\n",
                         number_of_rows, used_row,
                         number_of_cols, used_col);
    result->info_message("  %lu non-zeros\n", (unsigned long)num_nonzeros());
    result->info_message("  %s\n", hasBeenSolved ? "SOLVED" : "unsolved");
    
    // print my matrix
    result->info_message("Rows\n");
    for (size_t r =0; (int) r <= used_row; ++r)
    {
      print_row_iia(r);
    }
    result->info_message("Column non-zeros\n");
    for (size_t c=0; (int) c <= used_col && c < col_rows.size(); ++c)
    {
      print_col_iia(c);
    }
    result->info_message("Column values\n");
    for (size_t c=0; (int) c <= used_col && c < col_rows.size(); ++c)
    {
      if (parentProblem)
      {
        int pc = parentCols[c];
        result->info_message("  x%d (p%d) col", (int) c, pc);
      }
      else
        result->info_message("  x%d col", (int) c);
      result->info_message(" bounds=[%s,%s], goal=%f, solution=%d",
                           (numeric_limits<int>::lowest() == col_lower_bounds[c] ? "-inf" : to_string(col_lower_bounds[c]).c_str()),
                           (numeric_limits<int>::max()    == col_upper_bounds[c] ? "+inf" : to_string(col_upper_bounds[c]).c_str()),
                           goals[c],
                           col_solution[c]);
      if (has_tied_variables)
      {
        if (tied_variables[c].size()==1)
        {
          const auto c2 = tied_variables[c].front();
          result->info_message(" <tied to x%d = %d>", c2, col_solution[c2]);
        }
        else if (tied_variables[c].size()>1)
        {
          result->info_message(" [ MasterTie of x");
          for (auto c2 : tied_variables[c])
          {
            result->info_message(" %d", c2);
          }
          result->info_message(" ]");
        }
      }
      result->info_message("\n");
    }
    result->info_message("end IA problem\n\n");
  }
  
  void IncrementalIntervalAssignment::print_problem(string prefix, const vector<int> &col_order) const
  {
    print_problem(prefix);
    result->info_message("^^^^Above columns (variables) are permuted into this order: ");
    print_vec(result, col_order);
    result->info_message("\n");
    
    // future, do something fancy with copying and reordering
    // auto copy = *this;
  }
  
  void IncrementalIntervalAssignment::print_solution_summary(string prefix) const
  {
    if (prefix.empty())
      result->info_message("\n");
    else
      result->info_message("%s ",prefix.c_str());
    result->info_message("IIA solution summary\n");
    result->info_message("  %s\n", hasBeenSolved ? "SOLVED" : "unsolved");
    result->info_message("Solution is %s.\n", (verify_full_solution(true,true) ? "feasible" : "BAD"));
    // might be "BAD" because the tied variables are not assigned values yet...
  }
  
  void IncrementalIntervalAssignment::print_solution(string prefix) const
  {
    if (prefix.empty())
      result->info_message("\n");
    else
      result->info_message("%s ",prefix.c_str());
    result->info_message("IIA solution\n");
    result->info_message("  %s\n", hasBeenSolved ? "SOLVED" : "unsolved");
    for (int c=0; c <= used_col && c < (int)col_rows.size(); ++c)
    {
      print_solution_col(c);
    }
    result->info_message("Solution is %s.\n", (verify_full_solution(true, true) ? "feasible" : "BAD"));
    result->info_message("end IA solution\n\n");
  }
  
  void IncrementalIntervalAssignment::print_solution_row(int r, string prefix) const
  {
    if (!prefix.empty()) result->info_message("%s ",prefix.c_str());
    auto &row = rows[r];
    for (size_t i=0; i<row.cols.size();++i)
    {
      auto c = row.cols[i];
      print_solution_col(c);
    }
  }
  
  void IncrementalIntervalAssignment::print_solution_col(int c) const
  {
    result->info_message("  x%d = %d", c, col_solution[c]);
    if (has_tied_variables && tied_variables[c].size()==1)
    {
      const auto c2 = tied_variables[c].front();
      result->info_message(" <tied to x%d = %d>", c2, col_solution[c2]);
    }
    result->info_message("\n");
  }
  
  void RowSparseInt::print_row(IAResultImplementation *result) const
  {
    result->info_message( "row: " );
    for (size_t i=0;i<cols.size();++i)
    {
      auto v = vals[i];
      auto c = cols[i];
      result->info_message( "%+dx%d ", v, c );
    }
    result->info_message( "\n" );
  }
  
  
  void MatrixSparseInt::print_matrix_summary(string prefix) const
  {
    if (prefix.empty())
      result->info_message("\n");
    else
      result->info_message("%s ",prefix.c_str());
    result->info_message("sparse integer matrix\n");
    result->info_message(" %lu rows, %lu cols, %lu non-zeros\n",
                         (unsigned long) rows.size(),
                         (unsigned long) col_rows.size(),
                         (unsigned long) num_nonzeros());
  }
  
  void MatrixSparseInt::print_matrix(string prefix) const
  {
    if (prefix.empty())
      result->info_message("\n");
    else
      result->info_message("%s ",prefix.c_str());
    result->info_message("sparse integer matrix\n");
    result->info_message("  number rows = %lu\n"
                         "  number cols = %lu\n",
                         (unsigned long) rows.size(),
                         (unsigned long) col_rows.size());
    //result->info_message("stubbed out\n"); return;
    result->info_message("Rows\n");
    for (size_t r =0; r < rows.size(); ++r)
    {
      result->info_message("  %lu ", (unsigned long) r);
      rows[r].print_row(result);
    }
    result->info_message("Column non-zeros\n");
    for (size_t c=0; c < col_rows.size(); ++c )
    {
      result->info_message("  %lu col: ", (unsigned long) c);
      print_vec(result, col_rows[c]);
    }
    result->info_message("end sparse integer matrix\n\n");
  }
  
  
  void IncrementalIntervalAssignment::print() const
  {
    print_problem("");
  }
  
  bool IncrementalIntervalAssignment::tiny_subspace(int c, MatrixSparseInt &M, MatrixSparseInt &MN)
  {
    if (col_rows[c].empty())
      return false;
    
    // choose full columns of this matrix until we have a submatrix with a nontrivial nullspace
    set <int> tiny_cs, tiny_rs;
    tiny_cs.insert(c);
    tiny_rs.insert(col_rows[c].begin(), col_rows[c].end());
    
    SetValuesTiny val_tiny_Q;
    val_tiny_Q.tiny_cs = &tiny_cs;
    const double threshold=0.;
    
    // priority queue with replacement
    QWithReplacement Q(this,val_tiny_Q,threshold);
    for (auto r : col_rows[c])
    {
      Q.update(rows[r].cols);
    }
    
    do
    {
      // pick next column to add
      QElement t;
      if (Q.tip_top(t))
        break;
      
      // add it
      tiny_cs.insert(t.c);
      auto old_tiny_row_size = tiny_rs.size();
      tiny_rs.insert(col_rows[t.c].begin(), col_rows[t.c].end());

      // good?
      if (old_tiny_row_size == tiny_rs.size() && // we didn't add a new row
          tiny_cs.size()     > tiny_rs.size())   // chance of non-degenerate nullspace
      {
        // if there is any row with only one variable, then continue
        
        // find nullspace
        vector<int> tiny_cols(tiny_cs.begin(),tiny_cs.end());
        vector<int> tiny_rows(tiny_rs.begin(),tiny_rs.end());
        IncrementalIntervalAssignment *tiny_problem = new_sub_problem(tiny_rows, tiny_cols);
        vector<int> t_rref_col_order;
        tiny_problem->rref_improve(t_rref_col_order);
        vector<int> t_cols_dep, t_cols_ind, t_rows_dep;
        MatrixSparseInt tM(result);
        int tMaxMrow;

        // nullspace
        tiny_problem->categorize_vars(t_rref_col_order, t_cols_dep, t_cols_ind, t_rows_dep);
        const auto worked = tiny_problem->create_nullspace(t_cols_dep, t_cols_ind, t_rows_dep, tM, tMaxMrow);
        if ( worked && tMaxMrow > 0)
        {
          tiny_problem->cull_nullspace_tied_variable_rows(tM, tMaxMrow);

          // zzyk todo: we need to check whether any of the rows of tM actually contain c.
          // If not, we need to keep going. If so, we should only save those rows, not all of them!
          
          // convert tM to space of original matrix
          for (int tr = 0; tr < tMaxMrow; ++tr)
          {
            auto &trow = tM.rows[tr];
            RowSparseInt prow;
            for (auto tc : trow.cols)
            {
              prow.cols.push_back( tiny_problem->parentCols[tc] );
            }
            prow.vals = trow.vals;
            M.push_row(prow);
            MN.push_row(prow);
          }
          delete tiny_problem;
          return true; // it worked, try these
        }
      }
    }
    while (!Q.empty() && tiny_cs.size() < min(used_col, used_col/2 + 2)); // i.e. while the tiny space is less than have the current matrix size
    return false; // couldn't find a subspace small enough
  }

  
  void IncrementalIntervalAssignment::recursively_add_edge(int int_var_column,
                                                           int do_sum_even,
                                                           vector<int> &sub_rows,
                                                           vector<int> &sub_cols,
                                                           vector<int> &sub_row_array,
                                                           vector<int> &sub_col_array )
  {
    // add int_var_column to list - make sure not already in the list
    assert( !sub_col_array[ int_var_column ] );
    sub_cols.push_back( int_var_column );
    sub_col_array[ int_var_column ] = true;
    
    auto &non_zeros = col_rows[int_var_column];
    for ( auto row : non_zeros )
    {
      // if row not already in sub-problem,
      if ( !sub_row_array[ row ] )
      {
        auto &row_non_zeros = rows[row].cols;
        
        // should we add this row?
        bool do_add = true;
        // skip sum-even rows if sum-even phase
        if ( !do_sum_even )
        {
          // row is sum-even?
          for ( auto cc : row_non_zeros )
          {
            if (col_type[cc]==EVEN_VAR)
            {
              // result->debug_message("IncrementalIntervalAssignment::recursively_add_edge: column %d is a sum-even so skipping row %d.\n", cc, row );
              do_add = false;
              break;
            }
          }
          
          // row is a small loop? add it anyway
          // in the mapping phase, add a sum-even constraint row if the curves in the loop are likely very few intervals,
          //   smaller than the minimum number of intervals implied by the sum-even variable lower bound.
          // Probably a better solution is to just place a lower-bound on the intervals for individual edges, e.g. when one curve is the whole loop.
          // Checking small loops slows down the mapping phase solution, but can make a difference in improving quality.
          if (!do_add) // && 0 to zzyk disable small loop checking
          {
            int min_edge_count=0;
            int expected_edge_count=0;
            int sum_even_min=0;
            for ( size_t k = 0; k < row_non_zeros.size(); ++k)
            {
              int cc = row_non_zeros[k];
              int vv = rows[row].vals[k];
              int lo = col_lower_bounds[cc];
              if ( col_type[cc]==INT_VAR )
              {
                // if the true_count is less than the variable bound, use the variable bound
                const int expected_count = max(lo, (int) lround(goals[cc]) );
                expected_edge_count += vv * expected_count;
                min_edge_count      += vv * lo;
              }
              else
              {
                // cc is the sum-even variable, should probably check that there is only one
                sum_even_min = vv*lo;
              }
            }
            // add in rhs
            int hard_count = rhs[ row ];
            min_edge_count -= hard_count;
            expected_edge_count -= hard_count;
            // if it is possible and likely that the sum might be too small, then add the constraint
            if ( abs(min_edge_count) < abs(sum_even_min) && //strict inequality
                // Above: We don't need this row for lower bounds if above is false.
                // Below: can skip the row on the bet that intervals won't decrease by more than a factor of 2.
                abs(expected_edge_count) <= 2*abs(sum_even_min) )
            {
              do_add = true;
            }
          }
        }
        
        // actually add row to sub problem
        if ( do_add )
        {
          sub_rows.push_back( row );
          sub_row_array[ row ] = true;
          
          // recursively add cross-columns
          for ( auto cross_column : row_non_zeros )
          {
            if ( cross_column != int_var_column && !sub_col_array[ cross_column ] )
            {
              recursively_add_edge( cross_column, do_sum_even, sub_rows, sub_cols, sub_row_array, sub_col_array );
            }
          }
        }
      }
    }
  }
  
  IncrementalIntervalAssignment* IncrementalIntervalAssignment::new_sub_problem(const vector <int> &sub_rows, const vector <int> &sub_cols)
  {
    auto sub_problem = new IncrementalIntervalAssignment(result);
    
    sub_problem->should_adjust_solution_tied_variables = should_adjust_solution_tied_variables;
    
    sub_problem->parentProblem = this;
    sub_problem->parentRows.resize( sub_rows.size() );
    sub_problem->parentCols.resize( sub_cols.size() );
    
    sub_problem->lastCopiedCol =  (int)sub_cols.size();
    sub_problem->number_of_cols = (int)sub_cols.size();
    sub_problem->number_of_rows = (int)sub_rows.size();
    
    sub_problem->freeze_problem_size();
    
    sub_problem->parentCols = sub_cols;
    sub_problem->parentRows = sub_rows;
    
    sub_problem->Nullspace.rows.clear();
    sub_problem->Nullspace.col_rows.clear();
    sub_problem->Nullspace.col_rows.resize(sub_problem->number_of_cols);

    // uses the parentCols
    copy_bounds_to_sub( sub_problem );
    copy_submatrix( sub_rows, sub_cols, sub_problem );
    
    return sub_problem;
  }
  
  void IncrementalIntervalAssignment::subdivide_problem(vector<IncrementalIntervalAssignment*> &sub_problems,
                                                        vector<int> &orphan_cols,
                                                        bool do_sum_even, int row_min, int row_max )
  {
    // debug
    result->debug_message("---- Subdividing into independent sub-problems ----\n");
    CpuTimer subdivide_timer;
    
    // seen_columns[i] = TRUE if column i is in *any* subproblem.
    vector<int> seen_columns( number_of_cols, 0 );
    
    // map from row (col) of subproblem to row (col) of parent
    vector <int> sub_rows, sub_cols;
    
    // flag if a row or column is in the current subproblem yet.
    vector<int> sub_row_array ( number_of_rows, 0 );
    vector<int> sub_col_array ( number_of_cols, 0 );
    
    // if row_min and row_max are specified, then we just care about subproblems involving the variables in those rows.
    const bool row_subset = (row_min!=-1 && row_max!=-1);
    
    for (int e = 0; e <= used_col; ++e)
    {
      
      // if column not already in a subproblem
      // mapping phase: gather everything connected to an int var that isn't in a sum-even row
      // sum_even_phase: gather everything connected to dummy var
      if ( seen_columns[e] == 0 &&
          ((!do_sum_even && col_type[e]==INT_VAR) ||
           ( do_sum_even && col_type[e]==EVEN_VAR) ) )
      {
        
        // gather data for new subproblem
        
        // gather rows and columns
        sub_rows.clear();
        sub_cols.clear();
        recursively_add_edge( e, do_sum_even, sub_rows, sub_cols, sub_row_array, sub_col_array );
        
        if (sub_rows.empty())
        {
          orphan_cols.push_back(e);
          continue;
        }
        
        sort(sub_rows.begin(),sub_rows.end());
        sort(sub_cols.begin(),sub_cols.end());
        
        // add in the seen columns, zero out the sub_row and sub_col arrays
        for ( auto col : sub_cols )
        {
          seen_columns [ col ] = 1;
          sub_col_array[ col ] = 0;
        }
        
        for ( auto row : sub_rows )
        {
          sub_row_array[row] = 0;
        }
        
        // does the subproblem contain any of the new rows? If not, discard and continue
        if (row_subset)
        {
          bool subproblem_matters = false;
          for (auto row : sub_rows)
          {
            if (row>=row_min && row<=row_max)
            {
              subproblem_matters=true;
              break;
            }
          }
          if (!subproblem_matters)
            continue;
        }
        
        // create new sub problem
        IncrementalIntervalAssignment *sub_problem = new_sub_problem(sub_rows, sub_cols);
        sub_problems.push_back( sub_problem );
        
        // debug
        if (result->log_debug)
        {
          string title = "subproblem ";
          title += to_string(sub_problems.size());
          sub_problem->print_problem_summary(title);
        }
        // debug
        
      }
    }
    if (result->log_debug_time)
    {
      const auto subdivide_time = subdivide_timer.cpu_secs();
      result->info_message("\nDivided into %d subproblems in time = %f\n\n", (int) sub_problems.size(), subdivide_time);
    }
    else if (result->log_debug)
    {
      result->debug_message("\nDivided into %d subproblems.\n\n", (int) sub_problems.size());
    }
  }
  
  void IncrementalIntervalAssignment::delete_subproblems( vector<IncrementalIntervalAssignment*> &sub_problems )
  {
    for ( auto *s : sub_problems )
    {
      delete s;
    }
    sub_problems.clear();
  }
  
  void IncrementalIntervalAssignment::gather_solutions( vector<IncrementalIntervalAssignment*> &sub_problems,
                                                       bool want_sub_nullspaces, vector<int> &orphan_cols )
  {
    if (want_sub_nullspaces)
    {
      Nullspace.rows.clear();
      Nullspace.col_rows.clear();
      Nullspace.col_rows.resize(number_of_cols);
    }
    for (auto *sub : sub_problems)
    {
      assert(sub->parentProblem == this);
      if ( !sub->get_is_solved() )
        continue;
      
      // for each column, copy the solution of the sub-problem into the parent
      for (int e=0; e < sub->lastCopiedCol; e++ )
      {
        const auto c = sub->parentCols[ e ];
        assert( c >= 0 && c < number_of_cols );
        col_solution[ c ] = sub->col_solution[e];
      }

      // copy the nullspace into the parent nullspace
      if (want_sub_nullspaces)
      {
        auto &sub_nullspace = sub->Nullspace;

        // debug
        if (result->log_debug)
        {
          result->info_message("sub nullspace, with tied variables");
          sub_nullspace.print_matrix();
        }
        // convert sub_problem nullspace to parent nullspace variables.
        for (auto &row : sub_nullspace.rows)
        {
          RowSparseInt prow;
          prow.vals = row.vals;
          for (auto c : row.cols)
          {
            prow.cols.push_back(sub->parentCols[c]);
          }
//          // tied variables were already put back into the sub-problem, so the above loop already took care of it
//          if (sub->has_tied_variables)
//          {
//            const auto old_size = prow.cols.size();
//            for (int i=0; i<row.cols.size();++i)
//            {
//              const int c = row.cols[i];
//              const int v = row.vals[i];
//              for (int j=1; j<sub->tied_variables[c].size();++j)
//              {
//                int cc = sub->tied_variables[c][j];
//                prow.cols.push_back(sub->parentCols[cc]);
//                prow.vals.push_back(v);
//              }
//            }
//            // no need to sort if we didn't append any tied variables
//            if (old_size<prow.cols.size())
//              prow.sort();
//          }
          Nullspace.push_row(prow);
        }
        
        // debug
        if (result->log_debug)
        {
          Nullspace.print_matrix("parent nullspace after gathering subproblem nullspaces");
        }
      }
    }
    auto num_gathered_nullrows = Nullspace.rows.size();
    
    // find elementary nullspace vectors of interval variables what weren't in any subproblem
    // e.g., two paving surfaces sharing a curve, or a curve just in one paving surface and nowhere else
    if (want_sub_nullspaces)
    {
      for (auto c : orphan_cols)
      {
        if (col_type[c]==INT_VAR) // skip sum-even and dummy variables, we find nullspace vectors for those later
        {
          Nullspace.push_row( RowSparseInt( {c}, {1} ) );
        }
      }
      // debug
      if (result->log_debug)
      {
        result->debug_message("parent nullspace after appending %d orphan variables", (int) orphan_cols.size());
        Nullspace.print_matrix();
      }
    }
    auto num_elementary_nullrows = Nullspace.rows.size() - num_gathered_nullrows;

    if (result->log_debug)
    {
      if (want_sub_nullspaces)
        result->info_message("%lu parent nullspace rows: %lu from %lu subproblem nullspaces and %lu from elementary rows.\n",
                             Nullspace.rows.size(), num_gathered_nullrows, sub_problems.size(), num_elementary_nullrows);
      else
        result->info_message("Did not want sub nullspaces\n");
    }
  }
  
  
  bool IncrementalIntervalAssignment::verify_full_solution(bool print_unsatisfied_constraints, bool pretend_solved) const
  {
    if (!pretend_solved && !get_is_solved())
    {
      if (print_unsatisfied_constraints)
      {
        result->info_message("Interval Problem has not even been solved yet, so we can't verify the feasibility of the solution.\n");
      }
      return false;
    }
    
    bool rc = true;
    
    int num_out_of_bounds=0;
    int worst_out_of_bound=0;
    int num_non_even=0;
    
    // check integer variables
    for (int c = 0; c <= used_col; ++c)
    {
      const int sol = col_solution[c];
      
      const int lower = col_lower_bounds[c];
      const int upper = col_upper_bounds[c];
      
      // < out of bounds?
      if (sol<lower)
      {
        if (print_unsatisfied_constraints)
        {
          result->info_message("Variable x%d = %d < %d, below lower bound.\n",
                               c,
                               sol,
                               lower);
          ++num_out_of_bounds;
          worst_out_of_bound=max(worst_out_of_bound,lower-sol);
          rc=false;
        }
        else
          return false;
      }
      
      // > out of bounds?
      if (sol>upper)
      {
        if (print_unsatisfied_constraints)
        {
          result->info_message("Variable x%d = %d > %d, above upper bound.\n",
                               c,
                               sol,
                               upper);
          worst_out_of_bound=max(worst_out_of_bound,sol-upper);
          ++num_out_of_bounds;
          rc=false;
        }
        else
          return false;
      }
    }
    
    // check constraints
    int num_constraints_unsatisfied=0;
    for (int row = 0; row < number_of_rows; ++row)
    {
      auto &cols = rows[row].cols;
      
      // ignore rows with no variables
      if (cols.empty())
        continue;
      
      const auto b = rhs[row];
      ConstraintType ctype = constraint[row];
      
      // if its a sum-even, check that the sum is even, rather than checking the value of the sum-even variable
      
      // add up lhs
      int lhs(0);
      for (size_t i = 0; i < cols.size(); ++i)
      {
        const auto c = cols[i];
        const auto v = rows[row].vals[i];
        lhs+=col_solution[c]*v;
      }
      
      // check
      bool is_satisfied(false);
      switch (ctype)
      {
          // currently all inequalities are converted to EQ so this is the only case
        case EQ:
          if (lhs==b)
            is_satisfied=true;
          break;
        case LE:
          if (lhs<b)
            is_satisfied=true;
          break;
        case GE:
          if (lhs>b)
            is_satisfied=true;
          break;
        case EVEN:
          if ( (lhs - b) % 2 == 0)
          {
            is_satisfied = true;
            ++num_non_even;
          }
          break;
        default:
          ;
      }
      if (!is_satisfied)
      {
        if (print_unsatisfied_constraints)
        {
          result->info_message("Constraint row %d is not satisfied.  It involves ", row );
          print_row_iia(row);
          result->info_message("(lhs %d)",lhs);
          switch (ctype)
          {
            case EQ:
              result->info_message(" != ");
              break;
            case LE:
              result->info_message(" !<= ");
              break;
            case GE:
              result->info_message(" !>= ");
              break;
            case EVEN:
              result->info_message(" !Even ");
              break;
            default:
              result->info_message(" ? ");
          }
          result->info_message("(%d rhs)\n", b);
          ++num_constraints_unsatisfied;
          rc=false;
        }
        else
        {
          return false;
        }
      }
    }
    if (print_unsatisfied_constraints)
    {
      if (num_constraints_unsatisfied) result->info_message("%d constraints unsatisfied.\n",num_constraints_unsatisfied);
      if (num_non_even) result->info_message("%d rows non-even.\n",num_non_even);
      if (num_out_of_bounds) result->info_message("%d variables out-of-bounds.\n",num_out_of_bounds);
      if (worst_out_of_bound) result->info_message("%d is how far the worst variable is from its bounds.\n",worst_out_of_bound);
    }
    return rc;
  }
  
  void MatrixSparseInt::print_map(const BlockingCols &amap) const
  {
    for (auto m : amap)
    {
      const auto c = m.first;
      const auto lo =  m.second.first;
      const auto hi = m.second.second;
      if (lo==block_min)
        result->info_message(" %d(-inf,%d)", c, hi);
      else if (hi==block_max)
        result->info_message(" %d(%d,inf)", c, lo);
      else
        result->info_message(" %d(%d,%d)", c, lo, hi);
    }
    result->info_message("\n");
  }

  void IncrementalIntervalAssignment::copy_me( IncrementalIntervalAssignment *target )
  {
    // Matrix
    target->result   = result;
    target->rows     = rows;
    target->col_rows = col_rows;
    
    // EqualsB
    target->rhs        = rhs;
    target->constraint = constraint;

    // IncrementalIntervalAssignment
    target->hasBeenSolved    = hasBeenSolved;
    target->number_of_rows   = number_of_rows;
    target->used_row         = used_row;
    target->solved_used_row  = solved_used_row;
    target->number_of_cols   = number_of_cols;
    target->used_col         = used_col;
    target->solved_used_col  = solved_used_col;
    target->col_lower_bounds = col_lower_bounds;
    target->col_upper_bounds = col_upper_bounds;
    target->col_solution     = col_solution;
    target->goals            = goals;
    target->col_type         = col_type;
    target->parentProblem    = parentProblem;
    target->parentRows       = parentRows;
    target->parentCols       = parentCols;
    target->lastCopiedCol    = lastCopiedCol;
    
    // don't do anything with result
  }
  
  char letter_diff(const QElement &qe1, const QElement &qe2, double &v1, double &v2)
  {
    if (qe1.valueA!=qe2.valueA) {v1=qe1.valueA; v2=qe2.valueA; return 'A';}
    if (qe1.valueB!=qe2.valueB) {v1=qe1.valueB; v2=qe2.valueB; return 'B';}
    if (qe1.valueC!=qe2.valueC) {v1=qe1.valueC; v2=qe2.valueC; return 'C';}
    /*  all three same, use A*/ {v1=qe1.valueA; v2=qe2.valueA; return '=';}
  }
  
  int IncrementalIntervalAssignment::is_better( vector<QElement> &qA, vector<QElement> &qB ) const
  {
    // 1 if A<B by lexicographic min max
    //   1 if A[i]<B[i] for some i and A[j]<=B[j] for all j>i
    
    assert(qA.size()==qB.size());
    if (qA.empty())
      return 0;
    
    // start at back, largest element
    for (int i = ((int) qA.size())-1; i>=0; --i)
    {
      if (qA[i] < qB[i])
      {
        if (result->log_debug)
        {
          double vA, vB;
          auto c = letter_diff(qA[i], qB[i], vA, vB);
          const int max_place = ((int) qA.size())-1;
          const int place = max_place-i;
          result->debug_message("X<Y value%c %7.5g < %7.5g at place %d of %d, for X:x%d and Y:x%d\n", c, vA, vB, place, max_place, qA[i].c, qB[i].c );
        }
        return 1;
      }
      else if (qB[i] < qA[i])
      {
        if (result->log_debug)
        {
          double vA, vB;
          auto c = letter_diff(qA[i], qB[i], vA, vB);
          const int max_place = ((int) qA.size())-1;
          const int place = max_place-i;
          result->debug_message("X>Y value%c %7.5g > %7.5g at place %d of %d, for X:x%d and Y:x%d\n", c, vA, vB, place, max_place, qA[i].c, qB[i].c );
        }
        return -1;
      }
    }
    // every element was equal
    return 0;
  }

  int truncate_trailing_zeros( vector<double> &R )
  {
    size_t i=0;
    for( ; i<R.size() && R[i]>0.; ++i ) {}
    auto R_trailing_zeros = R.size()-i;
    R.resize(i);
    return (int) R_trailing_zeros;
  }

  int IncrementalIntervalAssignment::solution_is_better_than_Y(vector<int> &Y,
                                                               bool print_summary, bool print_detail )
  {
    return solution_X_is_better_than_Y( col_solution, Y, print_summary, print_detail);
  }
  int IncrementalIntervalAssignment::solution_X_is_better_than_Y(vector<int> &X, vector<int> &Y,
                                                                 bool print_summary, bool print_detail )
  {
    // compute R for X and Y
    // sort by increasing R
    // see which is lexicographically better

    // number of columns, safe, min
    int nc = std::min( {num_cols(), (int)X.size(), (int)Y.size()} );
    std::vector<int> cols0n(nc);
    std::iota(cols0n.begin(), cols0n.end(), 0);

    vector<QElement> qX, qY;
    qX.reserve(nc);
    qY.reserve(nc);

    // X
    col_solution.swap(X);
    compute_quality_ratio(qX,cols0n);
    col_solution.swap(X);

    // Y
    col_solution.swap(Y);
    compute_quality_ratio(qY,cols0n);
    col_solution.swap(Y);

    // lex compare
    // set valuesB to 0, since that just has to do with the goals
    // no need to re-sort, as valuesB is just a tiebreaker
    for (auto &q : qX)
    {
      q.valueB = 0.;
      q.valueC = 0.;
    }
    for (auto &q : qY)
    {
      q.valueB = 0.;
      q.valueC = 0.;
    }
    auto old_debug = result->log_debug;
    if (print_summary)
    {
      result->log_debug=true;
    }
    auto rc = is_better(qX,qY);
    
    // print which was better
    if (print_summary)
    {
      if (rc==1)
      {
        result->info_message("solution R(X) < R(Y)\n");
      }
      else if (rc==-1)
      {
        result->info_message("solution R(X) > R(Y)\n");
      }
      else
      {
        assert(rc==0);
        result->info_message("solution R(X) == R(Y)\n");
      }
    }
    // print actual values of X,Y, qX, qY
    if (print_detail)
    {
      result->info_message("X: ");
      print_vec(result,X);
      result->info_message("Y: ");
      print_vec(result,Y);

      // largest to smallest
      vector<double> RX(nc,0), RY(nc,0);
      int i=nc;
      for (auto qe : qX)
      {
        RX[--i] = qe.valueA;
      }
      i=nc;
      for (auto qe : qY)
      {
        RY[--i] = qe.valueA;
      }
      
      // truncate tailing zeros to avoid printing them.
      int X_trailing_zeros = truncate_trailing_zeros( RX );
      int Y_trailing_zeros = truncate_trailing_zeros( RY );

      result->info_message("valuesA, largest to smallest.\n");
      result->info_message("R(X): ");
      print_vec(result,RX,false);
      if (X_trailing_zeros)
        result->info_message(" ... and %d trailing zeros", X_trailing_zeros);
      result->info_message("\n");
      result->info_message("R(Y): ");
      print_vec(result,RY,false);
      if (Y_trailing_zeros)
        result->info_message(" ... and %d trailing zeros", Y_trailing_zeros);
      result->info_message("\n");

    }
    result->log_debug=old_debug;
    return rc;
  }

  // zzyk, hack put this here for fast recompile
  IncrementalIntervalAssignment::IncrementalIntervalAssignment( IAResultImplementation *result_ptr ) : MatrixSparseInt(result_ptr), Nullspace(result_ptr)
  {
    // all on
    if (0)
    {
      use_HNF_always = false; // false then try rref solution first
      use_HNF_goals = true;  // false c=1
      use_map_nullspace = true; // false use M only during sum-even phase
      use_best_improvement = true; // false use first improvement vector
      
      use_better_rref_step1_pivots_in_sumeven_phase = true; // Default true. If false, use the same pivot rules as in the mapping phase
      use_fixed_RatioR_straddle = true; // Default true. If false, reintroduce a bug for calculating RatioR for solutions less than 1 away from goals.
      use_smart_adjust_for_fractions = true; // Default true. If false, then we always round upwards
      use_Ratio_for_improvement_queue = true; // Default true. If false, then use RatioR for selecting which variable
      use_Nullspace_sumeven_to_flip_sign_of_blocker = true; // Default true. If true, then when trying to increase x0, if "x0 - x1" is blocked by x1, then add a "2x1 + y" Nullspace row to flip the sign of x1
    }

    // all off
    if (1)
    {
      use_HNF_always = false; // false then try rref solution first
      use_HNF_goals = false;  // false c=1
      use_map_nullspace = false; // false use M only during sum-even phase
      use_best_improvement = false; // false use first improvement vector
      
      use_better_rref_step1_pivots_in_sumeven_phase = false; // Default true. If false, use the same pivot rules as in the mapping phase
      use_fixed_RatioR_straddle = false; // Default true. If false, reintroduce a bug for calculating RatioR for solutions less than 1 away from goals.
      use_smart_adjust_for_fractions = false; // Default true. If false, then we always round upwards
      use_Ratio_for_improvement_queue = false; // Default true. If false, then use RatioR for selecting which variable
      use_Nullspace_sumeven_to_flip_sign_of_blocker = false; // Default true. If true, then when trying to increase x0, if "x0 - x1" is blocked by x1, then add a "2x1 + y" Nullspace row to flip the sign of x1
    }
  }

} // namespace
