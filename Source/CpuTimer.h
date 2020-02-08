//- Class: CpuTimer

#if !defined(CPU_TIMER_HPP)
#define CPU_TIMER_HPP

// #include "CubitUtilConfigure.h"

class CpuTimer
{
public:
  CpuTimer();			//- initialise to current system time


  double cpu_secs();		        //- return CPU time in seconds since last
                                //- call to cpu_secs();

  double cpu_elapsed();         //- return CPU time in seconds since 'birth'.

  double clock_secs();          //- return wall clock time in seconds since last
                                //- call to clock_secs();

  double clock_elapsed();       //- return wall clock time in seconds since 'birth'.

private:

  // the intial time at birth
  double mCPUInitial;
  // the last time computed
  double mCPU;

  // the real time at birth
  double mRealInitial;
  // the last real time computed
  double mReal;
  
  // get the current real time
  static double current_real_time();

  // get the current cpu time
  static double current_cpu_time();

};

#endif


