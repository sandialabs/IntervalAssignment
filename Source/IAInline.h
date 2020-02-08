// IAInline.h

// namespace IIA

using IIA_Internal::IncrementalIntervalAssignment;

inline
const IAResult *IA::get_result() const
{
  return result;
}

inline
IAResult::IAResult() :
solved(false),
constraints_satisfied(false),
bounds_satisfied(false),
optimized(false),
error(false)
{}

