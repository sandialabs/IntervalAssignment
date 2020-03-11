// IAQueue.h 
#ifndef IA_QUEUE_H
#define IA_QUEUE_H

#include <vector>
#include <map>
#include <set>

namespace IIA_Internal
{
  using std::vector;
  using std::map;
  using std::set;
  
  class IncrementalIntervalAssignment;
  class IAResultImplementation;
  
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
                       // c is last-resort arbitrary tiebreaker for less_for_sets
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
      return false;
      //      return
      //      ( valueA <  rhs.valueA ) ||
      //      ( valueA == rhs.valueA   && valueB <  rhs.valueB ) ||
      //      ( valueA == rhs.valueA   && valueB == rhs.valueB   && valueC <  rhs.valueC );
    }
    
    void print(const IAResultImplementation *result) const;
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
  
  // set_values function base class for assessing priority and quality of an interval setting for a curve
  // Updates value, valueB, valueC if intervals have changed.
  // return true if the new values are less than (lower priority than) the old values
  class SetValuesFn
  {
  public:
    // sets values based on current IIA.col_solution
    //   updates solution and returns true if solution has changed or was unset
    bool update_values(const IncrementalIntervalAssignment &iia, QElement &qe);
    // sets values based on current iia.col_solution
    void set_values(const IncrementalIntervalAssignment &iia, QElement &qe);
    // uses the passed in solution instead of iia's solution
    void set_values(const IncrementalIntervalAssignment &iia, QElement &qe, int solution);
  protected:
    // assumes qe.solution is already set, uses that instead of col_solution
    virtual void set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe ) = 0;
  };
  
  
  // specializations for differing quality/priority criteria
  
  
  // ratio  (current-increment)/goal                  for current>goal
  //                       goal/(current+increment)   for goal>current
  class SetValuesRatioR: public SetValuesFn
  {
  public:
    SetValuesRatioR() {}
  public:
    void set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe ) override;
  private:
    void set_values_goal(const IncrementalIntervalAssignment &iia, double g, QElement &qe );
  };
  
  // ratio  current/goal     for current>goal
  //           goal/current  for goal>current
  class SetValuesRatio: public SetValuesFn
  {
  public:
    SetValuesRatio() {}
  public:
    void set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe ) override;
    void set_values_goal(const IncrementalIntervalAssignment &iia, double g, QElement &qe );
  };
  
  // how far out-of-bounds is the variable value
  class SetValuesBounds: public SetValuesFn
  {
  public:
    SetValuesBounds() {}
  protected:
    void set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe ) override;
  };
  
  // used in RREF
  
  // for picking a column to eliminate in rref that has a single 1
  class SetValuesOneCoeff: public SetValuesFn
  {
  public:
    SetValuesOneCoeff() {}
  public:
    void set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe ) override;
  };
  
  // in order of increasing number of coefficients = number of rows the variable appears in
  class SetValuesNumCoeff: public SetValuesFn
  {
  public:
    SetValuesNumCoeff() {}
  public:
    void set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe ) override;
  };

  // in order of increasing number of coefficients = number of rows the variable appears in
  class SetValuesNumCoeffV2: public SetValuesFn
  {
  public:
    SetValuesNumCoeffV2() {}
  public:
    void set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe ) override;
  };

  // pick column to eliminate in rref, that has a small coeff, is in few rows, and is long
  class SetValuesCoeffRowsGoal: public SetValuesFn
  {
  public:
    SetValuesCoeffRowsGoal() {}
  public:
    void set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe ) override;
  };

// pick column to eliminate in rref, that has a small coeff, is in few rows, rows are short, and is long
class SetValuesCoeffRowsGoalVB: public SetValuesFn
{
public:
  SetValuesCoeffRowsGoalVB() {}
public:
  void set_values_implementation(const IncrementalIntervalAssignment &iia, QElement &qe ) override;
};

  // true if A<B by lexicographic min max
  bool is_better( vector<QElement> &qA, vector<QElement> &qB);
  
  // like a priority queue but you can update the priority of an element higher or lower. Based on set sorting.
  class QWithReplacement
  {
  public:
    
    QWithReplacement(const IncrementalIntervalAssignment *iia_in, SetValuesFn &f, double f_threshold)
    : val_fn(&f), threshold(f_threshold), iia(iia_in) {}
    
    void build(const vector<int> &cols);
    
    // returns true if done
    bool tip_top( QElement &t );
    
    // update the Q
    // cols are the indices of the columns that have changed
    // if changed_solution_only, then we only update elements whose col_solution has changed
    void update(const vector<int> &cols, bool changed_solution_only = true);
    
    bool empty() {return Q.empty();}
    
    // debug
    void print() const;
    
    size_t size() {return Q.size();}
    
  private:
    map<int,QElement> elements;
    set<QElement,QElement_less_for_sets> Q;
    SetValuesFn *val_fn = nullptr;
    double threshold = 0.0;
    const IncrementalIntervalAssignment *iia = nullptr;
    
    // conditionally add a q element for column c
    void add(int c);
  };
  
} // namespace

#endif
