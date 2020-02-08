#include "CpuTimer.h"

#include <chrono>
#include <time.h>

#if defined(_WIN32)
#  include <windows.h>
#endif

CpuTimer::CpuTimer()
{
  this->mCPUInitial = current_cpu_time();
  this->mCPU = this->mCPUInitial;

  this->mRealInitial = current_real_time();
  this->mReal = this->mRealInitial;
}

double CpuTimer::cpu_secs()
{
  double now = this->current_cpu_time();
  double delta = now - this->mCPU;
  this->mCPU = now;
  return delta;
}

double CpuTimer::clock_secs()
{
  double now = this->current_real_time();
  double delta = now - this->mReal;
  this->mReal = now;
  return delta;
}

double CpuTimer::cpu_elapsed()
{
  double delta;
  double now = this->current_cpu_time();
  delta = now - this->mCPUInitial;
  return delta;
}

double CpuTimer::clock_elapsed()
{
  double delta;
  double now = this->current_real_time();
  delta = now - this->mRealInitial;
  return delta;
}


double CpuTimer::current_real_time()
{
  std::chrono::high_resolution_clock::time_point t = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(t.time_since_epoch()).count() / 1.0e6;
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
