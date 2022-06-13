/* Worker.cpp -- (C) Mark Rodenkirch, February 2012

   This is the parent class for the worker.  There will be an instance of
   this class for each executing thread.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <vector>
#include <assert.h>
#include "Worker.h"
#include "Clock.h"

#include "../sieve/primesieve.hpp"

#ifdef USE_X86
#include "../x86_asm/fpu-asm-x86.h"
#endif

#ifdef WIN32
   static DWORD WINAPI ThreadEntryPoint(LPVOID threadInfo);
#else
   static void *ThreadEntryPoint(void *threadInfo);
#endif

Worker::Worker(uint32_t myId, App *theApp)
{
   char        name1[30], name3[30];

   sprintf(name1, "thread_%d_stats", myId);
   sprintf(name3, "thread_%d_worker", myId);
   
   ip_StatsLocker = new SharedMemoryItem(name1);
   ip_WorkerStatus = new SharedMemoryItem(name3);

   ib_Initialized = false;

   ip_App = theApp;
   ii_MyId = myId;

   il_PrimesTested = 0;
   il_LargestPrimeTested = 0;
   il_WorkerCpuUS = 0;

   ib_GpuWorker = false;

   ii_WorkSize = ip_App->GetCpuWorkSize();
   
   il_PrimeList = NULL;

   ii_MiniChunkSize = 0;
   il_MinPrimeForMiniChunkMode = PMAX_MAX_62BIT;
   il_MaxPrimeForMiniChunkMode = PMAX_MAX_62BIT;
   
#ifdef WIN32
   // Ignore the thread handle return since the parent process won't suspend
   // or terminate the thread.
   CreateThread(0,                        // security attributes (0=use parents)
                0,                        // stack size (0=use parents)
                ThreadEntryPoint,
                this,
                0,                        // execute immediately
                0);                       // don't care about the thread ID
#else
   pthread_create(&ih_Thread, NULL, &ThreadEntryPoint, this);
   pthread_detach(ih_Thread);
#endif
}

Worker::~Worker()
{
   delete ip_StatsLocker;
   delete ip_WorkerStatus;
   
   if (il_PrimeList != NULL)
      xfree(il_PrimeList);
}

#ifdef WIN32
DWORD WINAPI ThreadEntryPoint(LPVOID threadInfo)
#else
static void *ThreadEntryPoint(void *threadInfo)
#endif
{
   Worker *worker;

   worker = (Worker *) threadInfo;
   
   while (!worker->IsInitialized())
      Sleep(10);

   worker->StartProcessing();

#ifdef WIN32
   return 0;
#else
   pthread_exit(0);
#endif
}

void  Worker::AllocatePrimeList(void)
{
   if (il_PrimeList != NULL)
      return;
   
   // Get a little extra space because we want to use 0 to end the list.
   il_PrimeList = (uint64_t *) xmalloc((ii_WorkSize + 10) * sizeof(uint64_t));
}

// This is executed in a thread that is not the main thread
void  Worker::StartProcessing(void)
{
   uint64_t startTime, endTime;

#ifdef USE_X86
   uint16_t savedFpuMode;
#endif
   
   AllocatePrimeList();
   
   SetStatusWaitingForWork();
   
   while (true)
   {
      if (IsStatusWaitingForWork())
      {
         if (ip_App->IsSievingDone())
            break;
      
         Sleep(1);
         continue;
      }
                  
      SetStatusWorking();
      
      startTime = Clock::GetThreadMicroseconds();

#ifdef USE_X86
      // This is so the worker classes don't need to do this.
      savedFpuMode = fpu_mod_init();
#endif

      if (ii_MiniChunkSize > 0 && 
         il_PrimeList[0] > il_MinPrimeForMiniChunkMode &&
         il_PrimeList[ii_PrimesInList-1] < il_MaxPrimeForMiniChunkMode)      
         TestWithMiniChunks();
      else
         TestMegaPrimeChunk();

#ifdef USE_X86
      fpu_mod_fini(savedFpuMode);
#endif
      
      // We need to lock while updating these variables as the main thread can read them.
      ip_StatsLocker->Lock();
      
      endTime = Clock::GetThreadMicroseconds();

      il_WorkerCpuUS += (endTime - startTime);
      
      ip_StatsLocker->Release();

#if defined(USE_OPENCL) || defined(USE_METAL)      
      // If this is the special CPU worker and we no longer need it,
      // then we can break out of this loop and stop this thread.
      if (ii_MyId == 0 && il_LargestPrimeTested > ip_App->GetMinGpuPrime())
         break;
#endif

      if (!ib_GpuWorker && il_LargestPrimeTested > 100000)
      {
         double seconds = (double) ((endTime - startTime)) / 1000000.0;
         
         if (seconds < 1.0)
         {            
            while (seconds < 1.0)
            {
               seconds *= 2;
               ii_WorkSize *= 2;
            }
            
            if (ii_MyId == 0)
               ip_App->WriteToConsole(COT_OTHER, "Increasing worksize to %u since each chunk is tested in less than a second", ii_WorkSize);
            
            xfree(il_PrimeList);
            il_PrimeList = NULL;
            
            AllocatePrimeList();
         }
         
         if (seconds > 5.0)
         {
            while (seconds > 5.0)
            {
               seconds /= 2;
               ii_WorkSize /= 2;
            }
            
            while (ii_WorkSize % 16)
               ii_WorkSize++;
                        
            if (ii_MyId == 0)
               ip_App->WriteToConsole(COT_OTHER, "Decreasing worksize to %u since each chunk needs more than 5 seconds to test", ii_WorkSize);
            
            xfree(il_PrimeList);
            il_PrimeList = NULL;
            
            AllocatePrimeList();
         }
         
      }

      SetStatusWaitingForWork();
   }

   SetStatusStopped();
}

void   Worker::SetMiniChunkRange(uint64_t minPrimeForMiniChunkMode, uint64_t maxPrimeForMiniChunkMode, uint32_t chunkSize)
{
   if (chunkSize < 2 || chunkSize > 128)
      FatalError("Invalid number for chunk size");
   
   ii_MiniChunkSize = chunkSize;
   il_MinPrimeForMiniChunkMode = minPrimeForMiniChunkMode;
   il_MaxPrimeForMiniChunkMode = maxPrimeForMiniChunkMode;   
}

void    Worker::TestWithMiniChunks(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t countInChunk = 0;
   uint64_t miniPrimeChunk[128];
   uint32_t pIdx = 0;
   
   while (pIdx < ii_PrimesInList && il_PrimeList[pIdx] < maxPrime)
   {
      countInChunk = 0;
      
      while (pIdx < ii_PrimesInList && countInChunk < ii_MiniChunkSize)
      {
         miniPrimeChunk[countInChunk] = il_PrimeList[pIdx];
         
         pIdx++;
         
         countInChunk++;
      }
   
      if (countInChunk == 0)
         break;

      // Fill the rest of the array with useful values
      while (countInChunk < ii_MiniChunkSize)
      {
         miniPrimeChunk[countInChunk] = miniPrimeChunk[countInChunk-1];
         countInChunk++;
      }
      
      TestMiniPrimeChunk(miniPrimeChunk);
      
      il_PrimesTested += countInChunk;
      il_LargestPrimeTested = miniPrimeChunk[countInChunk-1];
   }
}

// Determine if there is a value x such that x^2 = 2 (mod p).
// Note that this does not find x.  That is what findRoot() will do.
// Since findRoot() is more computationally expensive this can
// quickly eliminate p from further consideration.
bool  Worker::IsQuadraticResidue(uint64_t n, uint64_t p)
{
	int32_t	   j = 1;
	uint64_t	   swap;

	while (n > 1)
	{
		if (n & 1)
		{
			swap = n;
			n = p;
			p = swap;
			if (((n & 3) == 3) && ((p & 3) == 3))
				j = -j;
			n = n % p;
		}
		else
		{
			n >>= 1;
			if (((p & 7) == 3) || ((p & 7) == 5))
				j = -j;
		}
	}

   return (j == 1);
}

// Compute 1/a (mod p)  Assumes 1 < a < p
uint64_t Worker::InvMod64(uint64_t a, uint64_t p)
{
  /* Thanks to the folks at mersenneforum.org.
     See http://www.mersenneforum.org/showthread.php?p=58252. */

   uint64_t ps1, ps2, q, r, t, dividend, divisor;
   uint32_t parity;
 
   assert(a < p);

   if (a < 3)
    return (a < 2) ? a : (p+1)/2;

   q = p / a;
   r = p % a;

   assert(r > 0);

   dividend = a;
   divisor = r;
   ps1 = q;
   ps2 = 1;
   parity = 0;

   while (divisor > 1)
   {
      r = dividend - divisor;
      t = r - divisor;
      if (r >= divisor) {
      q += ps1; r = t; t -= divisor;
      if (r >= divisor) {
         q += ps1; r = t; t -= divisor;
         if (r >= divisor) {
            q += ps1; r = t; t -= divisor;
            if (r >= divisor) {
            q += ps1; r = t; t -= divisor;
            if (r >= divisor) {
               q += ps1; r = t; t -= divisor;
               if (r >= divisor) {
                  q += ps1; r = t; t -= divisor;
                  if (r >= divisor) {
                  q += ps1; r = t; t -= divisor;
                  if (r >= divisor) {
                     q += ps1; r = t;
                     if (r >= divisor) {
                        q = dividend / divisor;
                        r = dividend % divisor;
                        q *= ps1;
                     } } } } } } } } }
         q += ps2;
       parity = ~parity;
       dividend = divisor;
       divisor = r;
       ps2 = ps1;
       ps1 = q;
   }

   assert(0 < ps1);
   assert(ps1 < p);

   return (parity) ? ps1 : p - ps1;
}

// Compute 1/a (mod p)  Assumes 1 < a < 2^32, a < p
uint64_t Worker::InvMod32(uint32_t a, uint64_t p)
{
  /* Thanks to the folks at mersenneforum.org.
     See http://www.mersenneforum.org/showthread.php?p=58252. */

   uint64_t ps1, ps2, q;
   uint32_t r, t, dividend, divisor;
   int32_t  sign = -1;

   if (a <= 1)
      return 0;

   q = p / a;
   r = p % a;

   dividend = a;
   divisor = r;
   ps1 = q;
   ps2 = 1;

   while (divisor > 1)
   {
      r = dividend - divisor;
      t = r - divisor;
      if (r >= divisor) {
         q += ps1; r = t; t -= divisor;
         if (r >= divisor) {
            q += ps1; r = t; t -= divisor;
            if (r >= divisor) {
               q += ps1; r = t; t -= divisor;
               if (r >= divisor) {
                  q += ps1; r = t; t -= divisor;
                  if (r >= divisor) {
                     q += ps1; r = t; t -= divisor;
                     if (r >= divisor) {
                        q += ps1; r = t; t -= divisor;
                        if (r >= divisor) {
                           q += ps1; r = t; t -= divisor;
                           if (r >= divisor) {
                              q += ps1; r = t;
                              if (r >= divisor) {
                                 q = dividend / divisor;
                              r = dividend % divisor;
                              q *= ps1;
      } } } } } } } } }
      q += ps2;
      sign = -sign;
      dividend = divisor;
      divisor = r;
      ps2 = ps1;
      ps1 = q;
   }

   assert(0 < ps1);
   assert(ps1 < p);

   return (sign > 0) ? ps1 : p - ps1;
}
