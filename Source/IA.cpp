// IA.cpp

#include "IA.h"
#include "IncrementalIntervalAssignment.h"
#include <limits>

namespace IIA
{
  using IIA_Internal::IncrementalIntervalAssignment;

  // can't be inline because we need the non-stubbed version of IncrementalIntervalAssignment
  IA::IA() : ia( new IncrementalIntervalAssignment(&result)), free_row(0)
  {}
  
  IA::~IA()
  {
    delete ia;
    ia = nullptr;
  }
}

