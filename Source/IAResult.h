// IAResult.h

#include <vector>
#include <string>

namespace IIA
{
  
  class IAResult
  {
  public:
    IAResult();
    virtual ~IAResult();
    
    bool solved;
    bool constraints_satisfied;
    bool bounds_satisfied;
    bool optimized;
    bool error;
    int log_start; // first index of log entries for latest method call
    std::vector<std::string> message_log;
    
    // mutable so caller can set these flags directly from get_result
    mutable bool log_error=true, log_warning=true, log_info=false, log_debug=false;
  };

  void print_log(IAResult &result);
}
