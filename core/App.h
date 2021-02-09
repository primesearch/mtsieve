/* App.h -- (c) Mark Rodenkirch, January 2018
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _APP_H
#define _APP_H

#define FRAMEWORK_VERSION     "2.2.0"

#define MAX_WORKERS           100
#define NO_WORKER             MAX_WORKERS + 99999
#define PMIN_MIN              1
#define PMAX_MAX_52BIT        (1ULL<<52)
#define PMAX_MAX_62BIT        (1ULL<<62)

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#include <stdio.h>
#include "main.h"
#include "Parser.h"

#ifdef HAVE_GPU_WORKERS
#include "../opencl/Device.h"
#endif

class App;

#include "Worker.h"
#include "SharedMemoryItem.h"

// Console output types
typedef enum { COT_OTHER = 1, COT_SIEVE } cotype_t;
typedef enum { AS_INITIALIZING, AS_RUNNING, AS_INTERRUPTED, AS_FINISHED } appstatus_t;
typedef enum { SS_NOT_STARTED, SS_SIEVING, SS_DONE } sievingstatus_t;

// Although declared here, this must be implemented by a child class of App
App *get_app(void);

class App
{  
public:
   App(void);
   virtual ~App(void) = 0;

   void              Banner(void);
   
   virtual void      Help(void) = 0;
   virtual void      AddCommandLineOptions(string &shortOpts, struct option *longOpts) = 0;
   virtual parse_t   ParseOption(int opt, char *arg, const char *source) = 0;
   virtual void      ValidateOptions(void) = 0;
   
   uint32_t          GetCpuWorkSize(void) { return ii_CpuWorkSize; };
   uint32_t          GetTotalWorkers(void) { return ii_TotalWorkerCount; };
   uint64_t          GetMaxPrimeForSingleWorker(void) { return il_MaxPrimeForSingleWorker; };
   
   void              SetRebuildNeeded(void) { ip_NeedToRebuild->SetValueNoLock(1); };
   
   uint32_t          GetCpuWorkerCount(void) { return ii_CpuWorkerCount; };
   uint32_t          GetGpuWorkerCount(void) { return ii_GpuWorkerCount; };
   
#ifdef HAVE_GPU_WORKERS
   uint32_t          GetGpuWorkSize(void) { return ii_GpuWorkGroupSize * ii_GpuWorkGroups; };
   uint32_t          GetGpuWorkGroups(void) { return ii_GpuWorkGroups; };
   void              SetGpuWorkGroupSize(uint32_t gpuWorkGroupSize) { ii_GpuWorkGroupSize = gpuWorkGroupSize; };
   
   Device           *GetDevice(void) { return ip_Device; }
   uint64_t          GetMinGpuPrime(void) { return il_MinGpuPrime; };
#endif
   
   uint64_t          GetMinPrime(void) { return il_MinPrime; };
   uint64_t          GetMaxPrime(void) { return il_MaxPrime; };
   
   void              ConvertNumberToShortString(uint64_t value, char *buffer);
   
   bool              IsSievingDone(void) { return (((sievingstatus_t) ip_SievingStatus->GetValueNoLock()) == SS_DONE); };
   bool              IsInterrupted(void) { return (((appstatus_t) ip_AppStatus->GetValueNoLock()) == AS_INTERRUPTED); };
   bool              IsRunning(void) { return (((appstatus_t) ip_AppStatus->GetValueNoLock()) == AS_RUNNING); };
   
   void              StopWorkers(void);
   void              Interrupt(void);

   void              Run(void);

   void              WriteToConsole(cotype_t consoleOutputType, const char *fmt, ...);
   void              WriteToLog(const char *fmt, ...);
   
   void              TellAllWorkersToRebuild(void);

protected:
   virtual void      ResetFactorStats(void) = 0;
   
   virtual void      GetReportStats(char *reportStats, double cpuUtilization) = 0;
   virtual void      LogStartSievingMessage(void) = 0;
   virtual void      Finish(const char *finishMethod, uint64_t elapsedTimeUS, uint64_t largestPrimeTested, uint64_t primesTested) = 0;
   virtual void      NotifyAppToRebuild(uint64_t largestPrimeTested) = 0;
   
   virtual void      PreSieveHook(void) = 0;
   virtual bool      PostSieveHook(void) = 0;

   void              SetBanner(string banner) { is_Banner = banner; };
   void              SetLogFileName(string logFileName);
   void              SetMaxPrimeForSingleWorker(uint64_t maxPrimeForSingleWorker) { il_MaxPrimeForSingleWorker = maxPrimeForSingleWorker; };
   
   void              SetAppMinPrime(uint64_t minPrime) { il_MinPrime = il_AppMinPrime = minPrime; };
   void              SetAppMaxPrime(uint64_t maxPrime) { il_MaxPrime = il_AppMaxPrime = maxPrime; };
   void              SetMinPrime(uint64_t minPrime);
   void              SetMaxPrime(uint64_t maxPrime, const char *why);
   void              SetMinGpuPrime(uint64_t minGpuPrime) { il_MinGpuPrime = minGpuPrime; };
   
   void              ParentHelp(void);
   void              ParentAddCommandLineOptions(string &shortOpts, struct option *longOpts);
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
   uint64_t          il_StartSievingUS;
   
   uint32_t          ii_CpuWorkSize;
   
#ifdef HAVE_GPU_WORKERS
   uint32_t          ii_GpuWorkGroupSize;
   uint32_t          ii_GpuWorkGroups;
#endif
   
   cotype_t          icot_LastConsoleOutputType;
                    
private:
   void              DeleteWorkers(void);
   void              CreateWorkers(uint64_t largestPrimeTested);
   
   uint64_t          PauseSievingAndRebuild(void);
   void              ReportStatus(void);
   uint32_t          GetNextAvailableWorker(bool useSingleThread, uint64_t &largestPrimeSieved);
   void              SetRebuildCompleted(void) { ip_NeedToRebuild->SetValueNoLock(0); };
   
   void              CheckReportStatus(void);
   
   void              Finish(void);
   
   uint32_t          ii_SavedSseMode;
   uint16_t          ii_SavedFpuMode;
   
#ifdef HAVE_GPU_WORKERS
   Device           *ip_Device;
#endif

   SharedMemoryItem *ip_Console;
   SharedMemoryItem *ip_AppStatus;
   SharedMemoryItem *ip_SievingStatus;
   SharedMemoryItem *ip_NeedToRebuild;
   
   Worker          **ip_Workers;
   
   bool              ib_SetMinPrimeFromCommandLine;
   
   uint32_t          ii_CpuWorkerCount;
   uint32_t          ii_GpuWorkerCount;
   uint32_t          ii_TotalWorkerCount;
   
   string            is_LogFileName;
   string            is_Banner;
   
   // These represent a number of milli-seconds
   uint64_t          il_TotalClockTime;
   uint64_t          il_WorkerClockTime;

   // These are starting points from which other times are computed
   uint64_t          il_InitUS;
   uint64_t          il_LastStatusReportUS;
   uint64_t          il_LastStatusPrimesTested;
   time_t            it_StartTime;

   uint64_t          il_TotalSieveUS;

   time_t            it_ReportTime;
};                

#endif
