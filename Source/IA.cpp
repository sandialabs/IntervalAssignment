// IA.cpp

#include "IA.h"
#include "IAImplementation.h"
#include <limits>

IAResult::IAResult()
{
	ia = new IAImplementation;
}

virtual ~IAResult()
{
	delete IAImplementation;
	IAImplementation = nullptr;
}



