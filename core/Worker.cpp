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
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/sse-asm-x86.h"

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

   il_LargestPrimeTested = 0;
   il_WorkerCpuUS = 0;
   ii_WorkSize = 0;

   il_PrimesTested = 0;

   ib_GpuWorker = false;

   ii_WorkSize = ip_App->GetCpuWorkSize();
   
#ifdef HAVE_GPU_WORKERS
   // This is only used by GPU workers.
   il_PrimeList = 0;
#endif

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

#ifdef HAVE_GPU_WORKERS
   if (il_PrimeList)
      xfree(il_PrimeList);
#endif
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

   worker->WaitForHandOff();

#ifdef WIN32
   return 0;
#else
   pthread_exit(0);
#endif
}

#ifdef HAVE_GPU_WORKERS
void  Worker::AllocatePrimeList(uint32_t workGroupSize)
{
   ip_App->SetGpuWorkGroupSize(workGroupSize);
   
   ii_WorkSize = ip_App->GetGpuWorkSize();
   
   il_PrimeList = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
}
#endif

// This is called by the main thread to determine if this worker
// is waiting for work.
bool  Worker::IsWaitingForWork(bool lockWorkerStatus)
{
   workerstatus_t status;
   
   ip_WorkerStatus->Lock();
   
   status = (workerstatus_t) ip_WorkerStatus->GetValueHaveLock();

   // If not waiting for work, then release the lock and return
   if (status != WS_WAITING_FOR_WORK)
   {
      ip_WorkerStatus->Release();
      return false;
   }

   // Only release the lock if we don't want to keep it.
   if (!lockWorkerStatus)
      ip_WorkerStatus->Release();

   return true;   
}

// This is called by the main thread when the thread is waiting for work.
uint64_t  Worker::ProcessNextPrimeChunk(uint64_t startFrom, uint64_t maxPrimeForChunk)
{
   uint64_t largestPrime;
   uint32_t primeCount;
   
   iv_Primes.clear();

   // Generate primes for this worker
   if (maxPrimeForChunk > 0 && maxPrimeForChunk > startFrom)
   {
      primeCount = primesieve::count_primes(startFrom+1, maxPrimeForChunk);
     
      // This assumes that worker threads need a count of primes that is a multiple of 4
      while (primeCount & 0x03)
         primeCount++;

      primesieve::generate_n_primes(primeCount, startFrom+1, &iv_Primes);
      
      largestPrime = maxPrimeForChunk;
   }
   else
   {
      primesieve::generate_n_primes(ii_WorkSize, startFrom+1, &iv_Primes);
   
      largestPrime = iv_Primes.back();
   }
   
   SetHasWorkToDo();
   
   ip_WorkerStatus->Release();
   
   // This could exceed maxPrime, but the TestPrimes() function will exit early if maxPrime is reached.
   return largestPrime;
}

// This is executed in a thread that is not the main thread
void  Worker::WaitForHandOff(void)
{
   uint64_t startTime;
   uint32_t savedSseMode;
   uint16_t savedFpuMode;
   
   SetWaitingForWork();
   
   while (true)
   {
      if (IsWaitingForWork(false))
      {
         if (ip_App->IsSievingDone())
            break;
      
         Sleep(10);
         continue;
      }
                  
      ip_WorkerStatus->SetValueNoLock(WS_WORKING);
      
      startTime = Clock::GetThreadMicroseconds();

      // This is so the worker classes don't need to do this.
      savedFpuMode = fpu_mod_init();
      savedSseMode = sse_mod_init();

#ifdef HAVE_GPU_WORKERS
      if (ib_GpuWorker)
      {
         vector<uint64_t>::iterator it = iv_Primes.begin();
         uint32_t idx = 0;
         
         while (it != iv_Primes.end())
         {
            il_PrimeList[idx] = *it;
            
            it++;
            idx++;
         }
      }
#endif

      if (ii_MiniChunkSize > 0 && 
         iv_Primes.front() > il_MinPrimeForMiniChunkMode &&
         iv_Primes.back() < il_MaxPrimeForMiniChunkMode)      
         TestWithMiniChunks();
      else
         TestMegaPrimeChunk();

      sse_mod_fini(savedSseMode);
      fpu_mod_fini(savedFpuMode);
      
      // We need to lock while updating these variables as the main thread can read them.
      ip_StatsLocker->Lock();

      il_WorkerCpuUS += (Clock::GetThreadMicroseconds() - startTime);
      
      ip_StatsLocker->Release();

#ifdef HAVE_GPU_WORKERS      
      // If this is the special CPU worker and we no longer need it,
      // then we can break out of this loop and stop this thread.
      if (ii_MyId == 0 && il_LargestPrimeTested > ip_App->GetMinGpuPrime())
         break;
#endif

      ip_WorkerStatus->SetValueNoLock(WS_WAITING_FOR_WORK);
   }

   ip_WorkerStatus->SetValueNoLock(WS_STOPPED);
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
   vector<uint64_t>::iterator it = iv_Primes.begin();
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t countInChunk = 0;
   uint64_t miniPrimeChunk[128];
   
   while (it != iv_Primes.end() && *it < maxPrime)
   {
      countInChunk = 0;
      
      while (it != iv_Primes.end() && countInChunk < ii_MiniChunkSize)
      {
         miniPrimeChunk[countInChunk] = *it;
         
         it++;
         
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
