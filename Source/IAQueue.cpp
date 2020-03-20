// IAQueue.cpp
#include "IAQueue.h"

#include <numeric>
#include <cmath>

#include "IncrementalIntervalAssignment.h"
#include "IAResultImplementation.h"

namespace  IIA_Internal
{

  using std::swap;
  using std::min;
  
  bool SetValuesFn::update_values(const IncrementalIntervalAssignment &iia, QElement &qe)
  {
    // no solution previously?
    if (!qe.solution_set)
    {
      qe.solution_set=true;
      qe.solution = iia.get_solution(qe.c);
      set_values_implementation(iia, qe );
      return true;
    }
    // queue priority based on a different solution value?
    if (qe.solution!=iia.get_solution(qe.c))
    {
      QElement qe_old = qe;
      qe.solution = iia.get_solution(qe.c);
      set_values_implementation(iia, qe );
      return ( qe < qe_old ); // ||  qe_old < qe
    }
    return false;
  }
  
  void SetValuesFn::set_values(const IncrementalIntervalAssignment &iia, QElement &qe, int solution)
  {
    qe.solution=solution;
    qe.solution_set=true;
    set_values_implementation(iia,qe);
  }
  
  void SetValuesFn::set_values(const IncrementalIntervalAssignment &iia, QElement &qe)
  {
    qe.solution=iia.get_solution(qe.c);
    qe.solution_set=true;
    set_values_implementation(iia,qe);
  }
  
  
  void SetValuesBounds::set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe )
  {
    if (qe.c>-1)
    {
      assert(qe.solution_set);
      auto diff = iia.compute_out_of_bounds(qe.c,qe.solution);
      qe.valueA = abs(diff);
      qe.dx = (diff > 0 ? 1 :  (diff < 0 ? -1 : 0) );
      double glo,ghi;
      // qe.valueB = (iia.compute_tied_goals(qe.c,glo,ghi) ? abs(ghi): abs(glo));
      bool tied = iia.compute_tied_goals(qe.c,glo,ghi);
      if (tied)
      {
        // tied-to variables get the highest priority of all the tied variables
        qe.valueB = abs(ghi);
        qe.valueC = abs(glo);
      }
      else
      {
        qe.valueB = abs(glo);
        qe.valueC = 0;
      }
    }
    else
    {
      qe.valueA = 0;
      qe.valueB = 0;
      qe.valueC = 0;
    }
  }
  
  void SetValuesRatioR::set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe )
  {
    auto &c=qe.c;
    if (c>-1)
    {
      double glo,ghi;
      if (iia.compute_tied_goals(c,glo,ghi))
      {
        QElement qh = qe;
        set_values_goal(iia,glo,qe);
        set_values_goal(iia,ghi,qh);
        // return the *higher* priority
        if (qe<qh)
        {
          swap(qe,qh);
        }
      }
      else
      {
        set_values_goal(iia,glo,qe);
      }
    }
    else
    {
      qe.valueA = 0.;
      qe.valueB = 0.;
      qe.valueC = -2.;
    }
  }
  
  void SetValuesRatioR::set_values_goal(const IncrementalIntervalAssignment &iia, double g, QElement &qe )
  {
    assert(qe.solution_set);
    if (g==0.) // its a "don't care" dummy variable
    {
      qe.valueA=0.;
      qe.valueB=0.;
      qe.valueC=0.;
      return;
    }
    const int x_int = qe.solution;
    const double x = (double) x_int;
    if (x+1.<g)
    {
      qe.valueA = (x+1>1e-2 ? g/(x+1.) : 1000.*g*(abs(x+1)+1.));
      qe.dx = 1; // want increment
    }
    else if (x-1.>g)
    {
      qe.valueA = (x-1.)/g; // (g>1e-2 ? (x-1.)/g : 1000*(x-1)*(abs(g)+1.));
      qe.dx = -1; // desire decrement
    }
    else
    {
      double pre,post;
      if (x>g)
      {
        pre  = x/g; // (g>1e-2   ? x/g      : 1000*x*(abs(g)+1.));
        post = (x-1>1e-2 ? g/(x-1.) : 1000.*g-x);
        qe.dx = -1;
      }
      else
      {
        pre  = (x>1e-2 ? g/x      : 1000.*g*(abs(x)+1.));
        post = (x+1.)/g; // (g>1e-2 ? (x+1.)/g : 1000*(x+1)*(abs(g)+1.));
        qe.dx = 1;
      }
      qe.valueA = pre-post; // was pre-post+1, but then goal of 1 has solution 2 and 3 to have the same value!!!
    }
    // a value less than 1 means the value would get worse than it is currently
    qe.valueB = -g; // prefer to change curves with small goals first, if ratios are equal
    qe.valueC = 0.;
    
    // if we desire increment and we are already at the upper bound, set the priority to 0
    if (qe.dx>0)
    {
      if (x_int>=iia.get_upper(qe.c))
      {
        qe.valueA=0.;
        qe.valueB=0.;
        qe.valueC=0.;
      }
    }
    // if we desire decrement and we are already at the lower bound, set the priority to 0
    else
    {
      if (x_int<=iia.get_lower(qe.c))
      {
        qe.valueA=0.;
        qe.valueB=0.;
        qe.valueC=0.;
      }
    }
  }
  
  void SetValuesRatio::set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe )
  {
    auto &c=qe.c;
    if (c>-1)
    {
      double glo,ghi;
      if (iia.compute_tied_goals(c,glo,ghi))
      {
        QElement qh = qe;
        set_values_goal(iia,glo,qe);
        set_values_goal(iia,ghi,qh);
        // return the *higher* priority
        if (qe<qh)
        {
          swap(qe,qh);
        }
      }
      else
      {
        set_values_goal(iia,glo,qe);
      }
    }
    else
    {
      qe.valueA = 0.;
      qe.valueB = 0.;
      qe.valueC = -2.;
    }
  }
  
  void SetValuesRatio::set_values_goal(const IncrementalIntervalAssignment &iia, double g, QElement &qe )
  {
    assert(qe.solution_set);
    if (g==0.) // its a "don't care" dummy variable
    {
      qe.valueA=0.;
      qe.valueB=0.;
      qe.valueC=0.;
      return;
    }
    
    int x_int = qe.solution;
    double x = (double) x_int;
    if (x<g)
    {
      qe.valueA = (x>1e-2 ? g/x : 1000.*g*(abs(x)+1.));
      qe.dx = 1; // want increment
    }
    else
    {
      qe.valueA = x/g; // scaling away from zero is done within set_goal now. (g>1e-2 ? x/g : 1000*x*(abs(g)+1.));
      qe.dx = -1; // desire decrement
    }
    qe.valueB = -g;  // prefer to change curves with small goals first
    qe.valueC = 0.;
    
    // if we desire increment and we are already at the upper bound, set the priority to 0
    if (qe.dx>0)
    {
      if (x_int>=iia.get_upper(qe.c))
      {
        qe.valueA=0.;
        qe.valueB=0.;
        qe.valueC=0.;
      }
    }
    // if we desire decrement and we are already at the lower bound, set the priority to 0
    else
    {
      if (x_int<=iia.get_lower(qe.c))
      {
        qe.valueA=0.;
        qe.valueB=0.;
        qe.valueC=0.;
      }
    }
  }
  
  
  void SetValuesOneCoeff::set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe )
  {
    // valuesA = has only one coefficient, of value 1
    // valuesB = has a goal of 0 (is a slack variable)
    // valuesC = prefer variables with a larger goal (to be dependent variables, leaving the ones with small goals to be the independent ones)
    qe.valueA=0.;
    qe.valueB=0.;
    qe.valueC=0.;
    
    const int c = qe.c;
    if (iia.col_rows[c].size()==1)
    {
      const int cr = iia.col_rows[c][0];
      const int c_coeff = iia.rows[cr].get_val(c);
      if (abs(c_coeff)==1)
      {
        qe.valueA=1.;
        qe.valueB= (qe.valueC==0. ? 1 : 0);
        // pick based on lower goal, unsure what's best
        double glo, ghi;
        iia.compute_tied_goals(c,glo,ghi);
        qe.valueC=glo;
      }
    }
  }
  
  void SetValuesNumCoeff::set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe )
  {
    // in only one row? pick it
    // number of non-one coefficients in all its rows (set union), fewer is better
    //
    // larger goal
    
    // valuesA = -number of coefficients
    // valuesB = -smallest coefficient magnitude
    // valuesC = larger goal (because if these are in fewer rref rows then they end up in more nullspace vectors)
    
    // exceptions for different types of variables: edge; sum-even coeff 2, sum-even coeff 1; slack equality, slack inequality.
    
    const int c = qe.c;
    // valueA
    qe.valueA = - ((int) iia.col_rows[c].size()); // number of rows
    
    // valueAA
    if ( qe.valueA < -1) // if in more than one row...
    {
      // number of non-sum-even vars in its rows, or vars with a non-one coefficient
      std::set<int> sum_evens;
      for (auto r : iia.col_rows[c])
      {
        for (auto c : iia.rows[r].cols )
        {
          if (iia.col_type[c]==IncrementalIntervalAssignment::EVEN_VAR)
          {
            sum_evens.insert(c);
          }
        }
        for (auto v : iia.rows[r].vals )
        {
          if (abs(v)>1)
          {
            sum_evens.insert(c);
          }
        }
      }
      qe.valueA -= 2*sum_evens.size();
    }
    
    // valueB
    {
      int smallest_coeff = numeric_limits<int>::max();
      for (auto r : iia.col_rows[c])
      {
        const int c_coeff = abs(iia.rows[r].get_val(c));
        smallest_coeff = min( smallest_coeff, c_coeff );
      }
      qe.valueB = -smallest_coeff;
    }
    // valueC
    auto ctype = iia.col_type[c];
    if (ctype==IncrementalIntervalAssignment::INT_VAR)
    {
      double glo, ghi;
      iia.compute_tied_goals(c,glo,ghi);
      qe.valueC=glo; // prefer longer ones
    }
    else if (ctype==IncrementalIntervalAssignment::EVEN_VAR)  // sum-even
    {
      if (fabs(qe.valueB)!=1)
      {
        // sum-even with coefficient 2, don't want to pick it
        qe.valueA = numeric_limits<double>::lowest()/4; // will never get picked, even if it is in only one row
      }
      qe.valueC = numeric_limits<double>::max()/2; // highly desirable: weak bounds, don't care quality, OK if in many nullspace rows
    }
    else // constraint slack variable
    {
      // distinguish equality slack (bad to pick) from inequality slack (OK to pick)
      // equality slack
      if (iia.col_lower_bounds[c]==iia.col_upper_bounds[c])
      {
        // undesireable, highly constrained to a single value. want it to be independent so it's isolated in its own nullspace row, not dependent on any other slacks
        qe.valueA = numeric_limits<double>::lowest()/2; // this way it will never get picked, even if it is in only one row
      }
      qe.valueC = -1;
    }
  }

  void SetValuesCoeffRowsGoal::set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe )
  {
    // if smallest coeff != 1, then valuesA = -inf
    // valuesA = -number of rows it appears in, fewer is better
    // valuesB = -number of variables in its rows, smaller is better
    // valuesC = has a goal of 0 (is a slack variable), secondarily prefer variables with a larger goal (to be dependent variables, leaving the ones with small goals to be the independent ones)
    qe.valueA=0.;
    qe.valueB=0.;
    qe.valueC=0.;
    
    const int c = qe.c;
    
    int best_coeff=numeric_limits<int>::max();
    for (int r : iia.col_rows[c])
    {
      const int c_coeff = abs(iia.rows[r].get_val(c));
      if (c_coeff<best_coeff)
      {
        best_coeff=c_coeff;
      }
    }
    if (best_coeff != 1)
    {
      qe.valueA = numeric_limits<double>::lowest()/2;
      return;
    }
    
    // number of rows
    auto &rows = iia.col_rows[c];
    qe.valueA = - ((double)rows.size());
    
    // number vars of in its rows
    std::set<int> vars;
    for (auto r : rows )
    {
      vars.insert( iia.rows[r].cols.begin(), iia.rows[r].cols.end() );
    }
    qe.valueB = - ((double) vars.size());
    
    double glo, ghi;
    iia.compute_tied_goals(c,glo,ghi);
    qe.valueC = (glo==0. ? numeric_limits<double>::max()/2 : glo);
  }

  void SetValuesTiny::set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe )
  {
    //     for rows with no other cols, need at least 1
    //     among those, want cols with fewest non-zeros in other rows
    const int c = qe.c;
    // c already chosen?
    if (tiny_cs->find(c)!=tiny_cs->end())
    {
      qe.valueA = numeric_limits<int>::lowest();
      qe.valueB = 0.;
      qe.valueC = 0.;
      return;
    }
    
    int num_new = 0, num_needed = 0;
    for (auto r : iia.col_rows[c])
    {
      auto &cols = iia.rows[r].cols;
      workspace.resize(std::max(cols.size(),tiny_cs->size()));
      // count how many tiny_columns are in this row already
      auto it = set_intersection(cols.begin(), cols.end(),
                                 tiny_cs->begin(), tiny_cs->end(),
                                 workspace.begin() );
      workspace.resize(it - workspace.begin());
      auto num_tiny_cols_in_row_already = workspace.size();
      if (num_tiny_cols_in_row_already==1) // need something in this row, just have one variable
        ++num_needed; // need some column in this row
      else if (num_tiny_cols_in_row_already==0)
        ++num_new;
      // else 2 or more cols already
    }
    qe.valueA = num_needed > 0 ? numeric_limits<int>::max()/4 : 0;
    qe.valueB = -num_new;
    // pick sum-even and slack vars first
    //   otherwise, pick variables in few columns, since that's closer to rref
    const int num_rows = (int) iia.col_rows[c].size();
    qe.valueC = -num_rows;
  }
  
  bool QWithReplacement::tip_top( QElement &t)
  {
    iia->get_result()->debug_message("Q size %d ", (int) Q.size());
    if (Q.empty())
    {
      t = QElement();
      return true;
    }
    
    // end of set is the highest-priority item
    auto Qback = --Q.end();
    t = *Qback;
    Q.erase(Qback);
    elements.erase(t.c);
    
    if (iia->get_result()->log_debug)
    {
      iia->get_result()->info_message("tip top ");
      t.print(iia->get_result());
    }
    return false;
  }
  
  void QWithReplacement::update(const vector<int> &cols, bool changed_solution_only)
  {
    if (/* DISABLES CODE */ (0) && iia->get_result()->log_debug)
    {
      iia->get_result()->info_message("Q before update ");
      print();
    }
    
    for (auto c : cols)
    {
      auto it = elements.find(c);
      // column wasn't in the queue
      if (it==elements.end())
      {
        add(c);
      }
      // column was in the queue, and solution changed
      else if ( !changed_solution_only || it->second.solution != iia->get_solution(c) )
      {
        QElement qe;
        qe.c = c;
        val_fn->set_values(*iia,qe);
        // remove from Q
        Q.erase(it->second);
        // replace element
        if (qe.valueA>threshold)
        {
          it->second=qe;
          // update Q
          Q.insert(qe);
        }
        // remove element
        else
        {
          elements.erase(it);
        }
      }
    }
    
    if (/* DISABLES CODE */ (0) && iia->get_result()->log_debug)
    {
      iia->get_result()->info_message("Q after update ");
      print();
    }
    
  }
  
  void QElement::print(const IAResultImplementation *result) const
  {
    result->info_message("qe c:%d ",c);
    if (!solution_set)
      result->info_message("unset");
    else
      result->info_message("sol:%d",solution);
    result->info_message(" A:%f B:%f C:%f\n",valueA,valueB,valueC);
  }
  
  void QWithReplacement::print() const
  {
    iia->get_result()->info_message("QWithReplacement %lu elements, threshold %f\n", (unsigned long) elements.size(), threshold);
    //  for (auto it : elements)
    //  {
    //    result->info_message("%d ", it.first);
    //    it.second.print();
    //  }
    iia->get_result()->info_message("Q %lu set\n", (unsigned long) Q.size() );
    for (auto e : Q)
    {
      e.print(iia->get_result());
    }
  }
  
  void QWithReplacement::add(int c)
  {
    QElement qe;
    qe.c = c;
    val_fn->set_values(*iia,qe);
    if (qe.valueA>threshold)
    {
      elements[c]=qe;
      Q.insert(qe);
      // return true;
    }
    // return false;
  }
  
  void QWithReplacement::build(const vector<int> &qcol)
  {
    for (auto c : qcol )
    {
      add(c);
    }
  }
  
}
