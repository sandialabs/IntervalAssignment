//- CpuTimer.h

#ifndef CPU_TIMER_HPP
#define CPU_TIMER_HPP

class CpuTimer
{
public:
  
  //- construction saves current system time
  CpuTimer();
  
  //- return CPU time in seconds since last call or construction.
  double cpu_secs();
  
  //- return CPU time in seconds since construction.
  double cpu_elapsed();
  
private:
  
  // the intial time at birth
  double mCPUInitial;
  // the last time computed
  double mCPU;
  
  // get the current cpu time
  static double current_cpu_time();
};

#endif
