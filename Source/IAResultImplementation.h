// IAResultImplementation.h
#ifndef IARESULT_IMPLEMENTATION_HPP 
#define IARESULT_IMPLEMENTATION_HPP

#include "IAResult.h"

#include <vector>
#include <set>
#include <map>

namespace IIA_Internal
{
  class IAResultImplementation : public IIA::IAResult
  {
    // io
  public:
    IAResultImplementation();
    virtual ~IAResultImplementation();
    
    // these call log_message
    void info_message(const char* format, ...) const;
    void warning_message(const char* format, ...) const;
    void error_message(const char* format, ...) const;
    void debug_message(const char* format, ...) const;
    
  protected:
    enum MessageType {INFO_MSG, WARNING_MSG, ERROR_MSG, DEBUG_MSG };
    void log_message(MessageType message_type, const std::string &prefix, const char* format, va_list args) const;
    // std::vector<char> message; // buffer for printing log messages
    mutable char message[2048]; // buffer for printing log messages
    mutable bool line_start = true;
    mutable MessageType last_message_type = INFO_MSG;
  };
  
  void print_vec(IAResultImplementation *result, const std::vector<int> &vec, bool lf = true);
  void print_vec(IAResultImplementation *result, const std::vector< std::pair<int,int> > &vec);
  void print_set(IAResultImplementation *result, const std::set< std::pair<int,int> > &aset);
  void print_map(IAResultImplementation *result, const std::map<int,int> &amap);
  
}

#endif

