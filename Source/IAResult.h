// IAResult.h

#ifndef IARESULT_H
#define IARESULT_H

#include <vector>
#include <string>

namespace IIA
{
  
  class IAResult
  {
  public:
    IAResult() {}
    virtual ~IAResult() {}
    
    bool
    solved = false,
    constraints_satisfied = false,
    bounds_satisfied = false,
    optimized = false;
    
    // The remainder are mutable so caller can set them from IIA::IA::get_result()
    
    // Feedback that something went wrong. Caller can clear it.
    mutable bool error = false;
    
    // Caller can choose where text messages go, if anywhere.
    mutable std::ostream *message_log = nullptr;
    
    // Flags that control the verbocity of the message_log.
    mutable bool
    log_error = true,
    log_warning=true,
    log_info=false,
    log_debug=false;
  };
  
}

#endif
