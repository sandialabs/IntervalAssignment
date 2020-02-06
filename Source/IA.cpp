// IA.cpp

#include "IA.h"
#include "IncrementalIntervalAssignment.h"
#include <limits>

namespace IIA
{
  using IIA_Internal::IncrementalIntervalAssignment;
  
  IAResult::IAResult() :
  solved(false),
  constraints_satisfied(false),
  bounds_satisfied(false),
  optimized(false),
  error(false)
  {}
  
  IA::IA() : ia( new IncrementalIntervalAssignment()), free_row(0)
  {}
  
  IA::~IA()
  {
    delete ia;
    ia = nullptr;
  }
  
}

