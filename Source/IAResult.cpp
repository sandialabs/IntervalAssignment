// IAResult.cpp
#include "IAResult.h"

#include <iostream>

namespace IIA
{
  void print_log(IAResult &result)
  {
    for (auto &s : result.message_log)
    {
      std::cout << s;
    }
  }
}

