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
    
    void print(IncrementalIntervalAssignment *iia);
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
  
  // set_values function. Updates value, valueB, valueC if intervals have changed.
  // return true if the new values are less than (lower priority than) the old values
  class SetValuesFn
  {
  public:
    // sets values based on current IIA.col_solution
    //   updates solution and returns true if solution has changed or was unset
    bool update_values(IncrementalIntervalAssignment &iia, QElement &qe);
    // sets values based on current iia.col_solution
    void set_values(IncrementalIntervalAssignment &iia, QElement &qe);
    // uses the passed in solution instead of iia's solution
    void set_values(IncrementalIntervalAssignment &iia, QElement &qe, int solution);
  protected:
    // assumes qe.solution is already set, uses that instead of col_solution
    virtual void set_values_implementation(IncrementalIntervalAssignment &iia, QElement &qe ) = 0;
    // virtual void print( const QElement &qe );
  };
  // specializations set the values based on differing criteria
  
  
  // ratio  (current-increment)/goal                  for current>goal
  //                       goal/(current+increment)   for goal>current
  class SetValuesRatioR: public SetValuesFn
  {
  public:
    SetValuesRatioR() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &iia, QElement &qe ) override;
  private:
    void set_values_goal(IncrementalIntervalAssignment &iia, double g, QElement &qe );
  };
  
  // ratio  current/goal     for current>goal
  //           goal/current  for goal>current
  class SetValuesRatio: public SetValuesFn
  {
  public:
    SetValuesRatio() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &iia, QElement &qe ) override;
    void set_values_goal(IncrementalIntervalAssignment &iia, double g, QElement &qe );
  };
  
  class SetValuesBounds: public SetValuesFn
  {
  public:
    SetValuesBounds() {}
  protected:
    void set_values_implementation(IncrementalIntervalAssignment &iia, QElement &qe ) override;
  };
  // also want a queue element for a constraint?
  
  // for picking a column to eliminate in rref that has a single 1
  class SetValuesOneCoeff: public SetValuesFn
  {
  public:
    SetValuesOneCoeff() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &iia, QElement &qe ) override;
  };
  
  // in order of increasing number of coefficients = number of rows the variable appears in
  class SetValuesNumCoeff: public SetValuesFn
  {
  public:
    SetValuesNumCoeff() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &iia, QElement &qe ) override;
  };
  
  
  // pick column to eliminate in rref, that has a small coeff, is in few rows, and is long
  class SetValuesCoeffRowsGoal: public SetValuesFn
  {
  public:
    SetValuesCoeffRowsGoal() {}
  public:
    void set_values_implementation(IncrementalIntervalAssignment &iia, QElement &qe ) override;
  };
  
  
  // true if A<B by lexicographic min max
  bool is_better( vector<QElement> &qA, vector<QElement> &qB);
  
  
  class QWithReplacement
  {
  public:
    
    QWithReplacement(IncrementalIntervalAssignment *iia_in, SetValuesFn &f, double f_threshold)
    : val_fn(&f), threshold(f_threshold), iia(iia_in) {}
    
    void build(const vector<int> &cols);
    
    // returns true if done
    bool tip_top( QElement &t );
    
    // update the Q
    // cols are the indices of the columns that have changed
    // col_solution is the new solution of the IncrementalIntervalAssignment
    void update(const vector<int> &cols);
    
    bool empty() {return Q.empty();}
    
    // debug
    void print();
    
    size_t size() {return Q.size();}
    
  private:
    map<int,QElement> elements;
    set<QElement,QElement_less_for_sets> Q;
    SetValuesFn *val_fn = nullptr;
    double threshold = 0.0;
    IncrementalIntervalAssignment *iia = nullptr;
    
    // conditionally add a q element for column c
    void add(int c);
  };
  
} // namespace

#endif
