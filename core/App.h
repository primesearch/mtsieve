/* App.h -- (c) Mark Rodenkirch, January 2018
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _APP_H
#define _APP_H

#define FRAMEWORK_VERSION     "2.2.3"

#define MAX_WORKERS           100
#define NO_WORKER             MAX_WORKERS + 99999
#define PMIN_MIN              1
#define PMAX_MAX_52BIT        (1ULL<<52)
#define PMAX_MAX_62BIT        (1ULL<<62)

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#include <stdio.h>
#include "main.h"
#include "Parser.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "GpuDevice.h"
#endif

#define MAX_PRIME_REPORT_COUNT   60

class App;

#include "Worker.h"
#include "SharedMemoryItem.h"

#include "../sieve/primesieve.hpp"

// Console output types
typedef enum { COT_OTHER = 1, COT_SIEVE } cotype_t;
typedef enum { AS_INITIALIZING, AS_RUNNING, AS_INTERRUPTED, AS_FINISHED } appstatus_t;
typedef enum { SS_NOT_STARTED, SS_SIEVING, SS_DONE } sievingstatus_t;

// Although declared here, this must be implemented by a child class of App
App *get_app(void);

typedef struct {
   uint64_t reportTimeUS;
   uint64_t primesTested;
} prime_report_t;

class App
{  
public:
   App(void);
   virtual ~App(void) = 0;

   void              Banner(void);
   
   virtual void      Help(void) = 0;
   virtual void      AddCommandLineOptions(std::string &shortOpts, struct option *longOpts) = 0;
   virtual parse_t   ParseOption(int opt, char *arg, const char *source) = 0;
   virtual void      ValidateOptions(void) = 0;
   
   uint32_t          GetCpuWorkSize(void) { return ii_CpuWorkSize; };
   bool              IsFixedCpuWorkSize(void) { return ib_FixedCpuWorkSize; };
   uint32_t          GetTotalWorkers(void) { return ii_TotalWorkerCount; };
   uint64_t          GetMaxPrimeForSingleWorker(void) { return il_MaxPrimeForSingleWorker; };
   
   void              SetRebuildNeeded(void) { ip_NeedToRebuild->SetValueNoLock(1); };
   
   uint32_t          GetCpuWorkerCount(void) { return ii_CpuWorkerCount; };
   uint32_t          GetGpuWorkerCount(void) { return ii_GpuWorkerCount; };
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          GetGpuPrimesPerWorker(void) { return ii_GpuWorkGroupSize * ii_GpuWorkGroups; };
   uint32_t          GetGpuWorkGroups(void) { return ii_GpuWorkGroups; };
   uint32_t          GetGpuWorkGroupSize(void) { return ii_GpuWorkGroupSize; };
   void              SetGpuWorkGroupSize(uint32_t gpuWorkGroupSize) { ii_GpuWorkGroupSize = gpuWorkGroupSize; };
      
   GpuDevice        *GetGpuDevice(void) { return ip_GpuDevice; }
   uint64_t          GetMinGpuPrime(void) { return il_MinGpuPrime; };
#endif
   
   uint64_t          GetMinPrime(void) { return il_MinPrime; };
   uint64_t          GetMaxPrime(void) { return il_MaxPrime; };
   
   void              ConvertNumberToShortString(uint64_t value, char *buffer);
   
   bool              IsSievingDone(void) { return (((sievingstatus_t) ip_SievingStatus->GetValueNoLock()) == SS_DONE); };
   bool              IsInterrupted(void) { return (((appstatus_t) ip_AppStatus->GetValueNoLock()) == AS_INTERRUPTED); };
   bool              IsRunning(void) { return (((appstatus_t) ip_AppStatus->GetValueNoLock()) == AS_RUNNING); };
   
   void              StopWorkers(void);
   void              Interrupt(const char *why);

   void              Run(void);

#ifdef __MINGW_PRINTF_FORMAT
   void              WriteToConsole(cotype_t consoleOutputType, const char *fmt, ...) __attribute__ ((format (__MINGW_PRINTF_FORMAT, 3, 4)));
   void              WriteToLog(const char *fmt, ...) __attribute__ ((format (__MINGW_PRINTF_FORMAT, 2, 3)));
#else
   void              WriteToConsole(cotype_t consoleOutputType, const char *fmt, ...) __attribute__ ((format (printf, 3, 4)));
   void              WriteToLog(const char *fmt, ...) __attribute__ ((format (printf, 2, 3)));
#endif
   
   void              TellAllWorkersToRebuild(void);

protected:
   virtual void      ResetFactorStats(void) = 0;
   
   virtual void      GetReportStats(char *reportStats, uint32_t maxStatsLength, double cpuUtilization) = 0;
   virtual void      LogStartSievingMessage(void) = 0;
   virtual void      Finish(const char *finishMethod, uint64_t elapsedTimeUS, uint64_t largestPrimeTested, uint64_t primesTested) = 0;
   virtual void      NotifyAppToRebuild(uint64_t largestPrimeTested) = 0;
   
   virtual void      PreSieveHook(void) = 0;
   virtual bool      PostSieveHook(void) = 0;

   void              SetBanner(std::string banner) { is_Banner = banner; };
   void              SetLogFileName(std::string logFileName);
   void              SetMaxPrimeForSingleWorker(uint64_t maxPrimeForSingleWorker)
   {
         if (il_MaxPrime < maxPrimeForSingleWorker)
            il_MaxPrimeForSingleWorker = il_MaxPrime;
         else
            il_MaxPrimeForSingleWorker = maxPrimeForSingleWorker;
   };
   
   void              SetAppMinPrime(uint64_t minPrime) { il_MinPrime = il_AppMinPrime = minPrime; };
   void              SetAppMaxPrime(uint64_t maxPrime) { il_MaxPrime = il_AppMaxPrime = maxPrime; };
   void              SetMinPrime(uint64_t minPrime);
   void              SetMaxPrime(uint64_t maxPrime, const char *why);
   void              SetMinGpuPrime(uint64_t minGpuPrime) { il_MinGpuPrime = minGpuPrime; };
   
   void              ParentHelp(void);
   void              ParentAddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParentParseOption(int opt, char *arg, const char *source);
   void              ParentValidateOptions(void);

   bool              IsRebuildNeeded(void) { return (ip_NeedToRebuild->GetValueNoLock() > 0); };
   
   void              GetWorkerStats(uint64_t &workerCpuUS, uint64_t &largestPrimeTestedNoGaps, uint64_t &largestPrimeTested, uint64_t &primesTested);
   uint64_t          GetLargestPrimeTested(bool finishedNormally);

   void              Sieve(void);
   
   virtual Worker   *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested) = 0;

   bool              ib_HaveCreatedWorkers;
   
   uint64_t          il_AppMinPrime;
   uint64_t          il_AppMaxPrime;
   uint64_t          il_MinPrime;
   uint64_t          il_MaxPrime;
   
   uint64_t          il_LargestPrimeSieved;
   
   // This represents the largest prime that must be tested by a single worker.  There is one
   // restriction, it must be a CPU worker.
   uint64_t          il_MaxPrimeForSingleWorker;

   uint64_t          il_SplitRangeSize;
   
   uint64_t          il_MinGpuPrime;
   uint64_t          il_StartSievingProcessUS;
   uint64_t          il_StartSievingUS;
   
   uint32_t          ii_CpuWorkSize;
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_GpuWorkGroupSize;
   uint32_t          ii_GpuWorkGroups;
#endif
   
   cotype_t          icot_LastConsoleOutputType;
                    
private:
   primesieve::iterator  ip_PrimeIterator;
   
   void              DeleteWorkers(void);
   void              CreateWorkers(uint64_t largestPrimeTested);
   
   uint64_t          PauseSievingAndRebuild(void);
   void              ReportStatus(void);
   uint32_t          GetNextAvailableWorker(bool useSingleThread, uint64_t &largestPrimeSieved);
   uint64_t          GetPrimesForWorker(uint32_t th);
   void              SetRebuildCompleted(void) { ip_NeedToRebuild->SetValueNoLock(0); };
   
   void              CheckReportStatus(void);
   
   void              Finish(void);
   void              GetPrimeStats(char *primeStats, uint64_t primesTested);

#ifdef USE_X86
   uint32_t          ii_SavedSseMode;
   uint16_t          ii_SavedFpuMode;
#endif

#if defined(USE_OPENCL) || defined(USE_METAL)
   GpuDevice        *ip_GpuDevice;
#endif

   SharedMemoryItem *ip_Console;
   SharedMemoryItem *ip_AppStatus;
   SharedMemoryItem *ip_SievingStatus;
   SharedMemoryItem *ip_NeedToRebuild;
   
   Worker          **ip_Workers;
   
   bool              ib_FixedCpuWorkSize;
   bool              ib_SetMinPrimeFromCommandLine;
   
   uint32_t          ii_CpuWorkerCount;
   uint32_t          ii_GpuWorkerCount;
   uint32_t          ii_TotalWorkerCount;
   
   std::string       is_LogFileName;
   std::string       is_Banner;
   
   // These represent a number of milli-seconds
   uint64_t          il_TotalClockTime;
   uint64_t          il_WorkerClockTime;

   // These are starting points from which other times are computed
   uint64_t          il_InitUS;
   time_t            it_StartTime;

   prime_report_t    ir_ReportStatus[MAX_PRIME_REPORT_COUNT];
   uint32_t          ii_LastStatusEntry;
   
   uint64_t          il_TotalSieveUS;

   time_t            it_ReportTime;
};                

#endif
