// IA.cpp

#include "IA.h"
#include "IAImplementation.h"
#include <limits>


IAResult::IAResult() : 
solved(false), 
constraints_satisfied(false), 
bounds_satisfied(false), 
optimized(false), 
error(false)
{}

IAResult::IAResult() : ia( new IAImplementation()), free_row(0)
{}

virtual ~IAResult()
{
	delete IAImplementation;
	IAImplementation = nullptr;
}



