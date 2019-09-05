/* Clock.cpp -- (C) Mark Rodenkirch, December 2013

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>

#ifdef WIN32
   #include <windows.h>
#else
   #include <sys/resource.h>
   #include <sys/time.h>
   #include <assert.h>
   #ifdef __APPLE__
      #include <mach/mach.h>
   #endif
#endif

#include <inttypes.h>

#include "Clock.h"

Clock::Clock()
{
}

Clock::~Clock(void)
{
}

/* Time elapsed since some fixed base time. */
uint64_t  Clock::GetCurrentMicrosecond(void)
{
#ifdef _WIN32
   FILETIME ft;
   ULARGE_INTEGER ns100;
   GetSystemTimeAsFileTime(&ft);
   ns100.u.LowPart = ft.dwLowDateTime;
   ns100.u.HighPart = ft.dwHighDateTime;
   return (uint64_t)ns100.QuadPart/10;
#else
   struct timeval t;
   gettimeofday(&t,0);
   return (uint64_t) t.tv_sec*1000000 + t.tv_usec;
#endif
}

/* Time consumed by the current thread. */
uint64_t  Clock::GetThreadMicroseconds(void)
{
#ifdef _WIN32
   FILETIME ft_create, ft_exit, ft_kernel, ft_user;
   ULARGE_INTEGER ns100_user;
   GetThreadTimes(GetCurrentThread(), &ft_create, &ft_exit, &ft_kernel, &ft_user);
   ns100_user.u.LowPart = ft_user.dwLowDateTime;
   ns100_user.u.HighPart = ft_user.dwHighDateTime;
   return (uint64_t)(ns100_user.QuadPart)/10;
#elif defined(__APPLE__)
   // gleaned from https://github.com/lh3/misc/blob/master/sys/runit/runlib.c
   thread_t                thisThread;
   thread_info_data_t      thinfo;
   mach_msg_type_number_t  thread_info_count;
   thread_basic_info_t     tbi;
   
   thisThread = mach_thread_self();
   thread_info_count = TASK_BASIC_INFO_COUNT;
   
   kern_return_t error = thread_info(thisThread, THREAD_BASIC_INFO, (thread_info_t) thinfo, &thread_info_count);
                
   assert(error == KERN_SUCCESS);
   
   tbi = (thread_basic_info_t)thinfo;

   return (uint64_t) (tbi->user_time.seconds + tbi->system_time.seconds) * 1000000 
      + (tbi->user_time.microseconds + tbi->system_time.microseconds);
#else
   struct rusage r;
   getrusage(RUSAGE_THREAD,&r);
   return (uint64_t)(r.ru_utime.tv_sec + r.ru_stime.tv_sec)*1000000
      + (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
#endif
}

/* Total processor time consumed by all threads of the current process. */
uint64_t  Clock::GetProcessMicroseconds(void)
{
#ifdef _WIN32
   FILETIME ft_create, ft_exit, ft_kernel, ft_user;
   ULARGE_INTEGER ns100_user;
   GetProcessTimes(GetCurrentProcess(), &ft_create, &ft_exit, &ft_kernel, &ft_user);
   ns100_user.u.LowPart = ft_user.dwLowDateTime;
   ns100_user.u.HighPart = ft_user.dwHighDateTime;
   return (uint64_t)(ns100_user.QuadPart)/10;
#else
   struct rusage r;
   getrusage(RUSAGE_SELF,&r);
   return (uint64_t)(r.ru_utime.tv_sec + r.ru_stime.tv_sec)*1000000
      + (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
#endif
}

