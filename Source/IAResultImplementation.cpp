// IAResultImplementation.cpp
#include "IAResultImplementation.h"
#include <stdio.h>
#include <iostream>
#include <cstdarg>

namespace IIA_Internal
{
  
  IAResultImplementation::IAResultImplementation()
  {
    // message.reserve(2048);
  }
  
  IAResultImplementation::~IAResultImplementation()
  {
    if (message_log!=nullptr)
      message_log->flush();
  }
  
  
  void IAResultImplementation::info_message(const char* format, ...)
  {
    if (!message_log || !log_info)
      return;
    
    va_list args;
    va_start(args, format);
    log_message( INFO_MSG, "", format, args);
    va_end(args);
  }
  
  void IAResultImplementation::warning_message(const char* format, ...)
  {
    if (!message_log || !log_warning)
      return;
    
    va_list args;
    va_start(args, format);
    log_message( WARNING_MSG, "WARNING: ", format, args);
    va_end(args);
  }

  void IAResultImplementation::error_message(const char* format, ...)
  {
    if (!message_log || !log_error)
      return;
    
    va_list args;
    va_start(args, format);
    log_message( ERROR_MSG, "ERROR: ", format, args);
    va_end(args);
  }

  void IAResultImplementation::debug_message(const char* format, ...)
  {
    if (!message_log || !log_debug)
      return;
    
    va_list args;
    va_start(args, format);
    log_message( DEBUG_MSG, "DEBUG: ", format, args);
    va_end(args);
  }
  
  void IAResultImplementation::log_message(MessageType message_type, const std::string &prefix, const char* format, va_list args)
  {
    if (!message_log)
      return;
    
    // need a newline and prefix?
    if (last_message_type != message_type)
    {
      if (!line_start)
      {
        (*message_log) << "\n";
        line_start = true;
      }
      last_message_type = message_type;
    }
    
    // convert format, args to the message buffer
      // put formatted output into message buffer
//      message.clear();
//      auto len = vsnprintf(message.data(), message.capacity()-4, format, args);
//      print_ellipsis = (len > message.capacity()-6);
//      message.resize(len);
    auto len = vsnprintf(message, 2048-4, format, args);
    bool print_ellipsis = (len > 2048-6);
    
    // conditionally print prefix
    if (line_start)
    {
      (*message_log) << prefix;
    }
    
    if (len==0)
      return;
    
    // print message
    {
      for (int i=0; i < len; ++i)
        (*message_log) << message[i];
    }
    
    // print "..." if the message was too long to fit in the buffer
    if (print_ellipsis)
    {
      (*message_log) << "..." << "\n";
    }
    
    if (message[len-1] == '\n')
    {
      line_start = true;
      // flush?
      message_log->flush();
    }
  }
  
  
  
  void print_vec(IAResultImplementation *result, const std::vector<int> &vec, bool lf)
  {
    if (!result->message_log)
      return;
    for (auto v : vec)
      result->info_message(" %d",v);
    if (lf)
      result->info_message("\n");
  }

  void print_vec(IAResultImplementation *result, const std::vector< std::pair<int,int> > &vec)
  {
    if (!result->message_log)
      return;
    for (auto v : vec)
      result->info_message(" %d(%d)", v.first, v.second);
    result->info_message("\n");
  }

  void print_set(IAResultImplementation *result, const std::set< std::pair<int,int> > &aset)
  {
    if (!result->message_log)
      return;
    for (auto s : aset)
      result->info_message(" %d(%d)", s.first, s.second);
    result->info_message("\n");
  }
  
  void print_map(IAResultImplementation *result, const std::map<int,int> &amap)
  {
    if (!result->message_log)
      return;
    for (auto m : amap)
      result->info_message(" %d(%d)", m.first, m.second);
    result->info_message("\n");
  }

  
} // namespace
