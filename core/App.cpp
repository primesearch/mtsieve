/* App.cpp -- (C) Mark Rodenkirch, January 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>
#include <vector>
#include <stdarg.h>

#include "App.h"
#include "Clock.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/sse-asm-x86.h"

// Do not change this as some parts of the framework assume that this is set to 60 seconds
#define REPORT_SECONDS        60
#define REPORT_STRFTIME_FORMAT "ETC %Y-%m-%d %H:%M"
#define LOG_STRFTIME_FORMAT    "%Y-%m-%d %H:%M:%S"

App::App(void)
{   
   ii_SavedFpuMode = fpu_mod_init();
   ii_SavedSseMode = sse_mod_init();
   
   il_TotalClockTime = 0;
   il_WorkerClockTime = 0;
   il_InitUS = Clock::GetCurrentMicrosecond();
   ii_CpuWorkerCount = 0;
   ii_GpuWorkerCount = 0;

   ip_Console = new SharedMemoryItem("console");
   ip_AppStatus = new SharedMemoryItem("appstatus");
   ip_SievingStatus = new SharedMemoryItem("sievestatus");
   ip_NeedToRebuild = new SharedMemoryItem("rebuild");
   
   icot_LastConsoleOutputType = COT_OTHER;

   ip_AppStatus->SetValueNoLock(AS_INITIALIZING);
   ip_SievingStatus->SetValueNoLock(SS_NOT_STARTED);

   // The application cannot exceeds these max values
   il_AppMinPrime = PMIN_MIN;
   il_AppMaxPrime = PMAX_MAX_62BIT;
   
   il_MinPrime = il_AppMinPrime;
   il_MaxPrime = il_AppMaxPrime;
   il_MaxPrimeForSingleWorker = 0;
   il_LargestPrimeSieved = 0;
   
   ii_CpuWorkSize = 1000000;
   
   // We won't know this until we create a kernel in the GPU
   il_MinGpuPrime = 0;
   ib_HaveCreatedWorkers = false;
   ib_SetMinPrimeFromCommandLine = false;

   ip_Workers = (Worker **) xmalloc((MAX_WORKERS + 1) * sizeof(Worker *));
   
#ifdef HAVE_GPU_WORKERS
   ip_Device = new Device();
   ii_GpuWorkGroupSize = 0;
   ii_GpuWorkGroups = 10;
#endif
}

App::~App(void)
{
   DeleteWorkers();

   xfree(ip_Workers);
   
   delete ip_Console;

   sse_mod_fini(ii_SavedSseMode);
   fpu_mod_fini(ii_SavedFpuMode);
}

void App::Banner(void)
{
   printf("%s\n", is_Banner.c_str());
}

void App::ParentHelp(void)
{
   char  maxPrime[30];
   
   ConvertNumberToShortString(il_AppMaxPrime, (char *) maxPrime);
   
   printf("-p --pmin=P0          sieve start: P0 < p (default %" PRIu64")\n", il_AppMinPrime);
   printf("-P --pmax=P1          sieve end: p < P1 (default %s)\n", maxPrime);
   printf("-w --worksize=w       primes per chunk of work (default %u)\n", ii_CpuWorkSize);
   printf("-W --workers=W        start W workers (default %u)\n", ii_CpuWorkerCount);

#ifdef HAVE_GPU_WORKERS
   printf("-g --gpuworkgroups=g  work groups per call to GPU (default %u)\n", ii_GpuWorkGroups);
   printf("-G --gpuworkers=G     start G GPU workers (default %u)\n", ii_GpuWorkerCount);
   
   ip_Device->Help();
#endif
}

void  App::ParentAddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   shortOpts += "p:P:w:W:";

   AppendLongOpt(longOpts, "pmin",          required_argument, 0, 'p');
   AppendLongOpt(longOpts, "pmax",          required_argument, 0, 'P');
   AppendLongOpt(longOpts, "worksize",      required_argument, 0, 'w');
   AppendLongOpt(longOpts, "workers",       required_argument, 0, 'W');
   
#ifdef HAVE_GPU_WORKERS
   shortOpts += "g:G:";
   
   AppendLongOpt(longOpts, "gpuworkgroups", required_argument, 0, 'g');
   AppendLongOpt(longOpts, "gpuworkers",    required_argument, 0, 'G');
   
   ip_Device->AddCommandLineOptions(shortOpts, longOpts);
#endif
}

// Returns:
//    0 if the option is OK
//   -1 if the argument is invalid
//   -2 if the argument is out of range
//   99 if the argument is not supported by this module
parse_t App::ParentParseOption(int opt, char *arg, const char *source)
{
   parse_t      status = P_UNSUPPORTED;
   uint64_t     minPrime;

#ifdef HAVE_GPU_WORKERS
   status = ip_Device->ParseOption(opt, arg, source);
   
   if (status != P_UNSUPPORTED)
      return status;
#endif

   switch (opt)
   {     
      case 'p':
         status = Parser::Parse(arg, il_AppMinPrime, il_AppMaxPrime-1, minPrime);
         SetMinPrime(minPrime);
         ib_SetMinPrimeFromCommandLine = true;
         break;

      case 'P':
         status = Parser::Parse(arg, il_AppMinPrime+1, il_AppMaxPrime, il_MaxPrime);
         break;

      case 'w':
         status = Parser::Parse(arg, 10, 1000000000, ii_CpuWorkSize);
         break;

#ifdef HAVE_GPU_WORKERS
      case 'W':
         status = Parser::Parse(arg, 0, MAX_WORKERS, ii_CpuWorkerCount);
         break;

      case 'g':
         status = Parser::Parse(arg, 1, 1000000, ii_GpuWorkGroups);
         break;

      case 'G':
         status = Parser::Parse(arg, 0, MAX_WORKERS, ii_GpuWorkerCount);
         break;
#else
      case 'W':
         status = Parser::Parse(arg, 1, MAX_WORKERS, ii_CpuWorkerCount);
         break;
#endif
   }

   return status;
}

void App::ParentValidateOptions(void)
{
   if (ii_CpuWorkerCount + ii_GpuWorkerCount > MAX_WORKERS)
      FatalError("Too many workers are configured.  The limit is %u", MAX_WORKERS);
      
   if (ii_CpuWorkerCount + ii_GpuWorkerCount == 0)
   {
#ifdef HAVE_GPU_WORKERS
      ii_GpuWorkerCount = 1;
#else
      ii_CpuWorkerCount = 1;
#endif
   }

   // Ensure that the number of primes is divisble by 32
   // Although not required by all CPU sieves, many of them optimize
   // with groups of 4, 16, or 32 primes.
   while (ii_CpuWorkSize & 0x10)
      ii_CpuWorkSize++;
   
   if (il_MinPrime >= il_MaxPrime)
      FatalError("pmin must be less than pmax");

   ii_TotalWorkerCount = ii_CpuWorkerCount + ii_GpuWorkerCount;
   
#ifdef HAVE_GPU_WORKERS
   ip_Device->Validate();
#endif
}

void App::ConvertNumberToShortString(uint64_t value, char *buffer)
{
   if (value == PMAX_MAX_52BIT)
   {
      strcpy(buffer, "2^52");
      return;
   }
   
   if (value == PMAX_MAX_62BIT)
   {
      strcpy(buffer, "2^62");
      return;
   }
   
   if (value % 1000 != 0)
   {
      sprintf(buffer, "%" PRIu64"", value);
      return;
   }
   
   // Convert to scientific notation
   uint32_t e = 0;

   while (value % 10 == 0)
   {
      value /= 10;
      e++;
   }
   
   sprintf(buffer, "%" PRIu64"e%u", value, e);
}

void App::SetLogFileName(string logFileName)
{ 
   is_LogFileName = logFileName;
   
   FILE *fPtr = fopen(is_LogFileName.c_str(), "a");

   if (!fPtr)
      FatalError("Cannot open log file <%s>", is_LogFileName.c_str());
   
   fclose(fPtr);
}

void  App::WriteToLog(const char *fmt, ...)
{
   va_list  args;
   time_t   theTime = time(NULL);
   char     timeBuffer[32];
   
   if (!strftime(timeBuffer, sizeof(timeBuffer), LOG_STRFTIME_FORMAT, localtime(&theTime)))
      timeBuffer[0] = '\0';

   ip_Console->Lock();
   
   FILE *fPtr = fopen(is_LogFileName.c_str(), "a");

   if (!fPtr)
      FatalError("Cannot open log file <%s>", is_LogFileName.c_str());
   
   fprintf(fPtr, "%s: ", timeBuffer);
   
   va_start(args, fmt);
   vfprintf(fPtr, fmt, args);
   va_end(args);

   fprintf(fPtr, "\n");

   fclose(fPtr);

   ip_Console->Release();
}

void  App::WriteToConsole(cotype_t consoleOutputType, const char *fmt, ...)
{
   va_list args;

   ip_Console->Lock();

   // If the previous message was a status message and this message is not the
   // same type of status message, trigger a line feed.
   if (icot_LastConsoleOutputType != COT_OTHER && icot_LastConsoleOutputType != consoleOutputType)
      printf("\n");

   va_start(args, fmt);
   vprintf(fmt, args);
   va_end(args);

   // For some outputs to the console we want to trigger a line feed.  An example
   // of this would be an important result discovered by the application.  For
   // other outputs we just want to trigger a carriage return so that the next
   // line can overwrite the line we just output to the console.  An example of
   // this would be a status message that indicates progress.
   if (consoleOutputType == COT_OTHER)
      printf("\n");
   else
      printf("\r");

   fflush(stdout);

   icot_LastConsoleOutputType = consoleOutputType;

   ip_Console->Release();
}

void  App::StopWorkers(void)
{
   int64_t  count = 1;
   int64_t  iter = 0;

   // This tells the Workers to stop as soon as possible.
   ip_SievingStatus->SetValueNoLock(SS_DONE);

   count = 1;
   while (count)
   {
      CheckReportStatus();
         
      Sleep(100);

      count = 0;
      
      for (uint32_t ii=0; ii<=ii_TotalWorkerCount; ii++)
      {
         // ip_Worker[0] is the special CPU worker (if we need one)
         if (ii == 0 && ip_Workers[0] == NULL)
            continue;
         
         if (!ip_Workers[ii]->IsStopped())
            count++;
      }

      if (count && ++iter > 60000)
      {
         WriteToConsole(COT_OTHER, "%d workers didn't stop after 10 minutes", (int32_t) count);
         exit(0);
      }
   }
}

void  App::Interrupt(void)
{
   ip_AppStatus->SetValueNoLock(AS_INTERRUPTED);

   if (il_LargestPrimeSieved > il_MaxPrime)
      WriteToConsole(COT_OTHER, "CTRL-C accepted.  Threads will stop after sieving to %llu", il_MaxPrime);
   else
      WriteToConsole(COT_OTHER, "CTRL-C accepted.  Threads will stop after sieving to %llu", il_LargestPrimeSieved);
}

void  App::SetMinPrime(uint64_t minPrime)
{
   if (ib_SetMinPrimeFromCommandLine)
      return;
   
   il_MinPrime = minPrime;
}

// Overrice the max prime to be sieved so that we can guarantee
// that all remaining terms are prime.
void  App::SetMaxPrime(uint64_t maxPrime, const char *why)
{
   WriteToConsole(COT_OTHER, "Changing p_max to %" PRIu64".  %s.", maxPrime, why);
   
   il_MaxPrime = maxPrime;
}
 
void  App::Run(void)
{
   bool isDone;

   ValidateOptions();

   ii_SavedFpuMode = fpu_mod_init();
   ii_SavedSseMode = sse_mod_init();
   
   do
   {
      // This gives applications a chance to change any configurations prior to sieving.
      PreSieveHook();
    
      CreateWorkers(il_MinPrime);

#ifdef HAVE_GPU_WORKERS
      if (ii_GpuWorkerCount > 0)
      {
         if (ii_GpuWorkGroupSize == 0)
            FatalError("could not compute optimal size of GPU work group");
         
         WriteToConsole(COT_OTHER, "GPU primes per worker is %u", GetGpuWorkSize());
      }
#endif

      Sieve();
      
      isDone = PostSieveHook();

      if (!isDone)
         DeleteWorkers();
   } while (!isDone);
}

void  App::Sieve(void)
{
   uint32_t th;
   uint64_t sieveStartUS;
   bool     useSingleThread;
   
   LogStartSievingMessage();

   ip_AppStatus->SetValueNoLock(AS_RUNNING);
   ip_SievingStatus->SetValueNoLock(SS_SIEVING);
   
   il_LargestPrimeSieved = il_MinPrime - 1;
   il_StartSievingUS = Clock::GetCurrentMicrosecond();

   it_StartTime = time(NULL);
   it_ReportTime = it_StartTime + REPORT_SECONDS;
   
   useSingleThread = (il_LargestPrimeSieved < il_MaxPrimeForSingleWorker);
   
   while (il_LargestPrimeSieved < il_MaxPrime && IsRunning())
   {
      th = GetNextAvailableWorker(useSingleThread, il_LargestPrimeSieved);

      // Stop if we couldn't get a worker.  This should only happen if there is 
      // a problem with the application or if ip_AppStatus is not AS_RUNNING.
      if (th == NO_WORKER)
         break;

      sieveStartUS = Clock::GetThreadMicroseconds();

      il_LargestPrimeSieved = ip_Workers[th]->ProcessNextPrimeChunk(il_LargestPrimeSieved, il_MaxPrimeForSingleWorker);

      il_TotalSieveUS += (Clock::GetThreadMicroseconds() - sieveStartUS);

      // If we are using a single thread, then this will effectively block other Workers
      // from getting work until both this Worker is done and the largest prime tested for
      // this workr exceeds the max prime for a single CPU thread.
      if (th == 0 || useSingleThread)
      {
         while (!ip_Workers[th]->IsWaitingForWork(false) && !ip_Workers[th]->IsStopped())
         {
            CheckReportStatus();
            
            Sleep(10);
         }

         useSingleThread = (ip_Workers[th]->GetLargestPrimeTested() < il_MaxPrimeForSingleWorker);
      }
   }
   
   Finish();
}

uint32_t  App::GetNextAvailableWorker(bool useSingleThread, uint64_t &largestPrimeSieved)
{
   uint32_t  th;
   uint32_t  attempts = 0;

   while (true)
   {
      attempts++;
      
      // Return NO_WORKER to indicate that we are returning without selecting
      // a worker as we want to stop processing.
      if (!IsRunning())
         return NO_WORKER;

      // If rebuilding, then the largest prime tested might be smaller
      // than the largest prime sieved so we have to update it.
      if (IsRebuildNeeded())
         largestPrimeSieved = PauseSievingAndRebuild();
      
      CheckReportStatus();
      
      for (th=0; th<=ii_TotalWorkerCount; th++)
      {
         // ip_Worker[0] is the special CPU worker (if we need one).
         // It will only stop when we can switch to the other workers.
         if (th == 0)
         {
            if (ip_Workers[0] == NULL)
               continue;

            if (ip_Workers[0]->IsStopped())
               continue;
         
            if (ip_Workers[0]->IsWaitingForWork(true))
               return th;
               
            continue;
         }

         if (ip_Workers[th]->IsGpuWorker())
         {
            if (useSingleThread)
               continue;
            
            if (largestPrimeSieved < il_MinGpuPrime)
               continue;
         }
            
         if (ip_Workers[th]->IsWaitingForWork(true))
            return th;
      }
      
      // If we didn't find one, sleep
      Sleep(10);
   }
}

uint64_t  App::PauseSievingAndRebuild(void)
{
   uint64_t  largestPrimeTested;

   StopWorkers();
      
   // We can pass true because all workers are stopped which means that
   // they have completed sieving their respective range of primes.
   largestPrimeTested = GetLargestPrimeTested(true);
    
   NotifyAppToRebuild();
      
   DeleteWorkers();

   // Reset the sieving status since StopWorkers had changed to SS_DONE
   ip_SievingStatus->SetValueNoLock(SS_SIEVING);
   
   CreateWorkers(largestPrimeTested);
   
   ResetFactorStats();
         
   SetRebuildCompleted();
   
   return largestPrimeTested;
}

void  App::CheckReportStatus(void)
{
   time_t theTime = time(NULL);
   
   if (theTime > it_ReportTime)
   {
      ReportStatus();
      it_ReportTime = theTime + REPORT_SECONDS;
   }
}

void  App::DeleteWorkers(void)
{
   for (uint32_t ii=0; ii<=ii_TotalWorkerCount; ii++)
   {
      // ip_Worker[0] is the special CPU worker (if we need one)
      if (ii == 0 && ip_Workers[0] == NULL)
         continue;
            
      ip_Workers[ii]->CleanUp();

      delete ip_Workers[ii];
   }
}

void  App::CreateWorkers(uint64_t largestPrimeTested)
{
   uint32_t w, th;
   bool allWaiting = false;
   
   ip_Workers[0] = NULL;
   
   if (ii_CpuWorkerCount == 0 && largestPrimeTested < il_MinGpuPrime)
   {
      // Worker "0" is only created if all of the following conditions are met:
      //    no CPU workers specified
      //    the application supports GPU workers
      //    the next prime tested cannot be tested on a GPU
      WriteToConsole(COT_OTHER, "Creating CPU worker to use until p >= %" PRIu64"", il_MinGpuPrime);
      ip_Workers[0] = CreateWorker(0, false, largestPrimeTested);
   }

   th = 1;
   
   // This will create the workers and start executing them
   // Create the CPU workers first, then the GPU workers
   for (w=0; w<ii_CpuWorkerCount; w++, th++)
      ip_Workers[th] = CreateWorker(th, false, largestPrimeTested);
   
   for (w=0; w<ii_GpuWorkerCount; w++, th++)
      ip_Workers[th] = CreateWorker(th, true, largestPrimeTested);
      
   ib_HaveCreatedWorkers = true;
   
   // We can't start until all workers waiting for work
   while (!allWaiting)
   {
      Sleep(10);
      allWaiting = true;
      
      for (th=0; th<=ii_TotalWorkerCount; th++)
      {
         // ip_Worker[0] is the special CPU worker (if we need one)
         if (th == 0 && ip_Workers[0] == NULL)
            continue;
            
         if (!ip_Workers[th]->IsWaitingForWork(false))
            allWaiting = false;
      }
   }
   
   il_LastStatusReportUS = Clock::GetCurrentMicrosecond();
   il_LastStatusPrimesTested = 0;
}

void  App::Finish(void)
{
   uint64_t    largestPrimeTestedNoGaps, largestPrimeTested, primesTested;
   uint64_t    workerCpuUS;
   uint64_t    processCpuUS;
   uint64_t    elapsedTimeUS;
   double      cpuUtilization;
   const char *finishMethod = (IsInterrupted() ? "interrupted" : "completed");
   
   // This won't return until all workers have completed processing work assigned to them.
   StopWorkers();

   processCpuUS = Clock::GetProcessMicroseconds();

   elapsedTimeUS = Clock::GetCurrentMicrosecond() - il_StartSievingUS;

   GetWorkerStats(workerCpuUS, largestPrimeTestedNoGaps, largestPrimeTested, primesTested);

   // Since all threads finished normally, there are no gaps thus we use largestPrimeTested.
   WriteToConsole(COT_OTHER, "Sieve %s at p=%" PRIu64".", finishMethod, largestPrimeTested);
   
   cpuUtilization = ((double) processCpuUS) / ((double) elapsedTimeUS);
   
#ifdef HAVE_GPU_WORKERS
   uint64_t    processGpuUS = ip_Device->GetGpuMicroseconds();
   
   // This outputs statistics regarding CPU time, i.e. time the CPU (or cores) spent
   // on this program, thus not including time spent working on other processes
   WriteToConsole(COT_OTHER, "CPU time: %.2f sec. (%.2f sieving) (%.2f cores) GPU time: %.2f sec. ",
            processCpuUS/1000000.0,
            il_TotalSieveUS/1000000.0,
            cpuUtilization,
            processGpuUS/1000000.0);
#else
   // This outputs statistics regarding CPU time, i.e. time the CPU (or cores) spent
   // on this program, thus not including time spent working on other processes
   WriteToConsole(COT_OTHER, "CPU time: %.2f sec. (%.2f sieving) (%.2f cores)",
            processCpuUS/1000000.0,
            il_TotalSieveUS/1000000.0,
            cpuUtilization);
#endif
   


   Finish(finishMethod, elapsedTimeUS, largestPrimeTested, primesTested);

   ip_AppStatus->SetValueNoLock(AS_FINISHED);
}

void  App::ReportStatus(void)
{
   double   percentDone;
   bool     havePercentDone;
   double   primeRate;
   double   cpuUtilization;
   uint32_t primePrecision;
   const char  *primeRateUnit;
   struct tm   *finish_tm;
   char     childStats[200];
   char     finishTimeBuffer[32];
   uint64_t currentUS, workerCpuUS;
   uint64_t processCpuUS, elapsedTimeUS;
   uint64_t largestPrimeTestedNoGaps, largestPrimeTested, primesTested;
   time_t   finish_date;

   currentUS = Clock::GetCurrentMicrosecond();

   processCpuUS = Clock::GetProcessMicroseconds();
   
   elapsedTimeUS = Clock::GetCurrentMicrosecond() - il_StartSievingUS;
   
   GetWorkerStats(workerCpuUS, largestPrimeTestedNoGaps, largestPrimeTested, primesTested);

#ifdef HAVE_GPU_WORKERS
   // Treat time spent in GPU as if spent in CPU
   cpuUtilization = ((double) processCpuUS + ip_Device->GetGpuMicroseconds()) / ((double) elapsedTimeUS);
#else
   cpuUtilization = ((double) processCpuUS) / ((double) elapsedTimeUS);
#endif
   
   GetReportStats(childStats, cpuUtilization);
   
   // Compute the percentage of the range we have completed
   if (largestPrimeTestedNoGaps == 0)
      havePercentDone = false;
   else
   {
      havePercentDone = true;
      
      if (IsInterrupted())
         percentDone = (double)(largestPrimeTestedNoGaps-il_MinPrime)/(il_LargestPrimeSieved-il_MinPrime);
      else
      {
         percentDone = (double)(largestPrimeTestedNoGaps-il_MinPrime)/(il_MaxPrime-il_MinPrime);
         
         if (percentDone < 1.0 && il_MaxPrime == il_AppMaxPrime)
            havePercentDone = false;
      }
      
   }
   
   // Compute how many primes are tested per second since the last report
   primeRate = (double)(primesTested-il_LastStatusPrimesTested)/(currentUS-il_LastStatusReportUS);
   
   primeRateUnit = "M";
   if (primeRate < 1.0) primeRate *= 1000.0, primeRateUnit = "K";
   if (primeRate < 1.0) primeRate *= 1000.0, primeRateUnit = "";

   primePrecision = 0;
   if (primeRate < 1000.0) primePrecision = 1;
   if (primeRate < 100.0)  primePrecision = 2;
   if (primeRate < 10.0)   primePrecision = 3;
   
   // Calculate ETC.
   finishTimeBuffer[0] = '\0';
   if (havePercentDone && percentDone > 0.0)
   {
      finish_date = (time_t) (it_StartTime + (time(NULL)-it_StartTime)/percentDone);
      finish_tm = localtime(&finish_date);
      if (!finish_tm || !strftime(finishTimeBuffer, sizeof(finishTimeBuffer), REPORT_STRFTIME_FORMAT, finish_tm))
         finishTimeBuffer[0] = '\0';
   }

   if (strlen(childStats) > 0)
   {   
      if (!havePercentDone)
         WriteToConsole(COT_SIEVE, "  p=%" PRIu64", %.*f%s p/sec, %s                            ",
                        largestPrimeTestedNoGaps, primePrecision, primeRate, primeRateUnit, childStats);
      else
         WriteToConsole(COT_SIEVE, "  p=%" PRIu64", %.*f%s p/sec, %s, %.1f%% done. %s           ",
                        largestPrimeTestedNoGaps, primePrecision, primeRate, primeRateUnit,
                        childStats, 100.0*percentDone, finishTimeBuffer);
   }
   else
   {   
      if (!havePercentDone)
         WriteToConsole(COT_SIEVE, "  p=%" PRIu64", %.*f%s p/sec                                ",
                        largestPrimeTestedNoGaps, primePrecision, primeRate, primeRateUnit);
      else
         WriteToConsole(COT_SIEVE, "  p=%" PRIu64", %.*f%s p/sec, %%.1f%% done. %s              ",
                        largestPrimeTestedNoGaps, primePrecision, primeRate, primeRateUnit,
                        100.0*percentDone, finishTimeBuffer);
   }

   il_LastStatusReportUS = currentUS;
   il_LastStatusPrimesTested = primesTested;
}

// When the program finishes, then the largest prime tested across all workers
// tells us that all primes below that value have been tested.  This is because
// none of the workers can be terminated while working on a chunk of primes.
// If we are still sieving, the largest prime tested might not be correct
// as threads will finish out of sequence, but that is okay.
void  App::GetWorkerStats(uint64_t &workerCpuUS, uint64_t &largestPrimeTestedNoGaps, uint64_t &largestPrimeTested, uint64_t &primesTested)
{
   uint64_t workerLargestPrimeTested;
   
   workerCpuUS = 0;
   primesTested = 0;
   largestPrimeTested = 0;
   largestPrimeTestedNoGaps = PMAX_MAX_62BIT;
         
#ifdef HAVE_GPU_WORKERS
   // Treat time spent in GPU as if spent in CPU
   workerCpuUS += ip_Device->GetGpuMicroseconds();
#endif

   for (uint32_t ii=0; ii<=ii_TotalWorkerCount; ii++)
   {
      // ip_Worker[0] is the special CPU worker (if we need one)
      if (ii == 0)
         continue;

      ip_Workers[ii]->LockStats();

      workerLargestPrimeTested = ip_Workers[ii]->GetLargestPrimeTested();

      // Ignore worker if it hasn't done any work.
      if (workerLargestPrimeTested > 0)
      {
         // If this worker is waiting for work, then aasume that at least one worker
         // will be working on the next chunk, so use its stats instead.
         if (!ip_Workers[ii]->IsWaitingForWork(false))
         {
            // If there are multiple workers, this will be the largest prime tested
            // where we know that all primes less than this prime have been tested.
            if (workerLargestPrimeTested < largestPrimeTestedNoGaps)
               largestPrimeTestedNoGaps = workerLargestPrimeTested;
         }
         
         if (workerLargestPrimeTested > largestPrimeTested)
            largestPrimeTested = workerLargestPrimeTested;
      }
      
      primesTested += ip_Workers[ii]->GetPrimesTested();
      workerCpuUS += ip_Workers[ii]->GetWorkerCpuUS();

      ip_Workers[ii]->ReleaseStats();
   }
      
   if (ip_Workers[0] != NULL)
   {
      ip_Workers[0]->LockStats();

      if (largestPrimeTested == 0)
      {
         largestPrimeTestedNoGaps = ip_Workers[0]->GetLargestPrimeTested();
         largestPrimeTested = ip_Workers[0]->GetLargestPrimeTested();
      }
      
      primesTested += ip_Workers[0]->GetPrimesTested();
      workerCpuUS += ip_Workers[0]->GetWorkerCpuUS();

      ip_Workers[0]->ReleaseStats();
   }
   
   if (largestPrimeTestedNoGaps == PMAX_MAX_62BIT)
      largestPrimeTestedNoGaps = largestPrimeTested;
}

uint64_t  App::GetLargestPrimeTested(bool finishedNormally)
{
   uint64_t largestPrimeTested = (finishedNormally ? 0 : il_AppMaxPrime);
   uint64_t workerLargestPrimeTested;
   
   if (!ib_HaveCreatedWorkers)
      return il_MinPrime;
   
   for (uint32_t ii=0; ii<=ii_TotalWorkerCount; ii++)
   {
      // ip_Worker[0] is the special CPU worker (if we need one)
      if (ii == 0)
         continue;

      ip_Workers[ii]->LockStats();
      
      workerLargestPrimeTested = ip_Workers[ii]->GetLargestPrimeTested();
      
      // Ignore worker if it didn't do any work.
      if (workerLargestPrimeTested > 0)
      {
         // If the program finished normally, then the largest prime tested across all
         // workers tells us that all primes below that value have been tested.  If not,
         // then there might be gaps, so get the smallest prime from the workers that
         // has been tested.
         if (finishedNormally)
         {
            if (workerLargestPrimeTested > largestPrimeTested)
               largestPrimeTested = workerLargestPrimeTested;
         }
         else
         {
            if (workerLargestPrimeTested < largestPrimeTested)
               largestPrimeTested = workerLargestPrimeTested;
         }
      }
      
      ip_Workers[ii]->ReleaseStats();
   }
   
      
   // If only the special worker has done work, then we will get the larger prime from it.
   if (largestPrimeTested == 0 && ip_Workers[0] != NULL)
   {
      ip_Workers[0]->LockStats();

      largestPrimeTested = ip_Workers[0]->GetLargestPrimeTested();
      
      ip_Workers[0]->ReleaseStats();
   }
   
   return largestPrimeTested;
}

