#include "CpuTimer.h"

#include <time.h>

#if defined(_WIN32)
#  include <windows.h>
#endif

CpuTimer::CpuTimer()
{
  this->mCPUInitial = current_cpu_time();
  this->mCPU = this->mCPUInitial;
}

double CpuTimer::cpu_secs()
{
  const double now = this->current_cpu_time();
  const double delta = now - this->mCPU;
  this->mCPU = now;
  return delta;
}

double CpuTimer::cpu_elapsed()
{
  const double now = this->current_cpu_time();
  const double delta = now - this->mCPUInitial;
  return delta;
}

double CpuTimer::current_cpu_time()
{
#ifdef _WIN32
  FILETIME a,b,c,d;
  if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0)
  {
    return (double)(d.dwLowDateTime |
                    ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
  }
  else
  {
    return 0;
  }
#else
  return (double)clock() / CLOCKS_PER_SEC;
#endif
}
