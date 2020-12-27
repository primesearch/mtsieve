/* MultiFactorialWorker.cpp -- (C) Mark Rodenkirch, October 2012

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <math.h>
#include "MultiFactorialWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/avx-asm-x86.h"
#include "../core/inline.h"
#include "../core/MpArithVector.h"

extern "C" int mfsieve(uint32_t start, uint32_t mf, uint32_t minmax, uint64_t *P);
extern "C" int multifactorial(uint32_t start, uint32_t mf, uint32_t minmax, uint64_t *P);

MultiFactorialWorker::MultiFactorialWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_MultiFactorialApp = (MultiFactorialApp *) theApp;
   
   ii_MinN = ip_MultiFactorialApp->GetMinN();
   ii_MaxN = ip_MultiFactorialApp->GetMaxN();
   ii_MultiFactorial = ip_MultiFactorialApp->GetMultiFactorial();

   ib_Initialized = true;
}

void  MultiFactorialWorker::CleanUp(void)
{
}

void  MultiFactorialWorker::TestMegaPrimeChunk(void)
{
   if (ii_MultiFactorial == 1)
      TestFactorial();
   else
      TestMultiFactorial();
}

void  MultiFactorialWorker::TestFactorial(void)
{
   uint64_t  ps[4], maxPrime = ip_App->GetMaxPrime();
   vector<uint64_t>::iterator it = iv_Primes.begin();
   uint32_t  n;
   // if i <= n_pair then (i - 1) * i < p. Compute n! = (2 * 3) * (4 * 5) * ... * ((n - 1) * n)
   uint32_t n_pair = std::max(2u, std::min(ii_MinN, uint32_t(sqrt(double(*it))) & ~1u));
      
   while (it != iv_Primes.end())
   {
      ps[0] = *it;
      it++;
      
      ps[1] = *it;
      it++;
      
      ps[2] = *it;
      it++;
      
      ps[3] = *it;
      it++;

      MpArithVec mp(ps);

      const MpResVec one = mp.one();
      const MpResVec minus_one = mp.sub(mp.zero(), one);
      const MpResVec two = mp.add(one, one);
      const MpResVec four = mp.add(two, two);
      const MpResVec eight = mp.add(four, four);

      // ri = residue of i, rf = residue of i!
      MpResVec ri = one, rf = one;
      // residue of i * (i + 1), the step is (i + 2) * (i + 3) - i * (i + 1) = 4 * i + 6
      MpResVec r_ixip1 = mp.zero(), r_step = mp.add(four, two);

      // Factorial with pairs of numbers: i! = ((i - 1) * i) * (i - 2)!
      for (n = 2; n < n_pair; n += 2)
      {
         r_ixip1 = mp.add(r_ixip1, r_step);
         r_step = mp.add(r_step, eight);
         rf = mp.mul(rf, r_ixip1);
      }

      // Factorial: i! = i * (i - 1)!
      ri = mp.toMp(n_pair - 1);
      for (n = n_pair; n < ii_MinN; ++n)
      {
         ri = mp.add(ri, one);
         rf = mp.mul(rf, ri);
      }

      // Factorial and check if i! = +/-1
      for (; n <= ii_MaxN; ++n)
      {
         ri = mp.add(ri, one);
         rf = mp.mul(rf, ri);

         if (MpArithVec::at_least_one_is_equal(rf, one) | MpArithVec::at_least_one_is_equal(rf, minus_one))
         {
            for (size_t k = 0; k < VECTOR_SIZE; ++k)
            {
               if (rf[k] == one[k])
                  if (ip_MultiFactorialApp->ReportFactor(ps[k], n, -1))
                     VerifyFactor(true, ps[k], n, -1);
                  
               if (rf[k] == minus_one[k]) 
                  if (ip_MultiFactorialApp->ReportFactor(ps[k], n, +1))
                     VerifyFactor(true, ps[k], n, +1);
            }
         }
      }
      
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] >= maxPrime)
         break;
   }
}

void  MultiFactorialWorker::TestMultiFactorial(void)
{
   uint64_t  ps[4], qs[4], mfrs[4];
   uint64_t  mOne[4], pOne[4];
   uint64_t  ri[4], rf[4];
   uint64_t  maxPrime = ip_App->GetMaxPrime();
   uint32_t  maxNFirstLoop = ii_MinN - ii_MultiFactorial;
   uint32_t  n, startN;
   
   vector<uint64_t>::iterator it = iv_Primes.begin();
      
   while (it != iv_Primes.end())
   {
      ps[0] = *it;
      it++;
      
      ps[1] = *it;
      it++;
      
      ps[2] = *it;
      it++;
      
      ps[3] = *it;
      it++;

      for (uint32_t idx=0; idx<4; idx++)
      {
         pOne[idx] = mmmOne(ps[idx]);
         mOne[idx] = mmmSub(0, pOne[idx], ps[idx]);
         qs[idx] = mmmInvert(ps[idx]);
         mfrs[idx] = mmmN(ii_MultiFactorial, ps[idx]);
      }
      
      for (startN=1; startN<=ii_MultiFactorial; startN++)
      {
         // If startN is odd and mf is even, then i!mf is always odd, thus
         // startN!mf+1 and startN!mf-1 are always even.
         if (!(ii_MultiFactorial & 1) && (startN & 1))
            continue;
            
         for (uint32_t idx=0; idx<4; idx++)
         {
            if (startN > ps[idx])
               ri[idx] = mmmN(startN%ps[idx], ps[idx]);
            else
               ri[idx] = mmmN(startN, ps[idx]);
            rf[idx] = ri[idx];
         }
 
         // At this time we have:
         //    ri = residual of startN (mod p)
         //    rf = residual of startN!mf (mod p)
         
         n = startN + ii_MultiFactorial;
         for (; n<maxNFirstLoop; n+=ii_MultiFactorial)
         {
            ri[0] = mmmAdd(ri[0], mfrs[0], ps[0]);
            rf[0] = mmmMulmod(rf[0], ri[0], ps[0], qs[0]);
            
            ri[1] = mmmAdd(ri[1], mfrs[1], ps[1]);
            rf[1] = mmmMulmod(rf[1], ri[1], ps[1], qs[1]);
            
            ri[2] = mmmAdd(ri[2], mfrs[2], ps[2]);
            rf[2] = mmmMulmod(rf[2], ri[2], ps[2], qs[2]);
            
            ri[3] = mmmAdd(ri[3], mfrs[3], ps[3]);
            rf[3] = mmmMulmod(rf[3], ri[3], ps[3], qs[3]);

            //if (ps[0] == 110125033) printf("n=%u %llu %llu\n", n, ri[0], rf[0]);
            //if (ps[1] == 110125033) printf("n=%u %llu %llu\n", n, ri[1], rf[1]);
            //if (ps[2] == 110125033) printf("n=%u %llu %llu\n", n, ri[2], rf[2]);
            //if (ps[3] == 110125033) printf("n=%u %llu %llu\n", n, ri[3], rf[3]);
         }
         
         // At this time we have:
         //    ri = residual of mn (mod p)
         //    rf = residual of mn!mf (mod p)
         // where mn is the max n < ii_MinN for this starting n
         for (; n <=ii_MaxN; n+=ii_MultiFactorial)
         {
            ri[0] = mmmAdd(ri[0], mfrs[0], ps[0]);
            rf[0] = mmmMulmod(rf[0], ri[0], ps[0], qs[0]);
            
            ri[1] = mmmAdd(ri[1], mfrs[1], ps[1]);
            rf[1] = mmmMulmod(rf[1], ri[1], ps[1], qs[1]);
            
            ri[2] = mmmAdd(ri[2], mfrs[2], ps[2]);
            rf[2] = mmmMulmod(rf[2], ri[2], ps[2], qs[2]);
            
            ri[3] = mmmAdd(ri[3], mfrs[3], ps[3]);
            rf[3] = mmmMulmod(rf[3], ri[3], ps[3], qs[3]);
            
            //if (ps[0] == 110125033) printf("n=%u %llu %llu\n", n, ri[0], rf[0]);
            //if (ps[1] == 110125033) printf("n=%u %llu %llu\n", n, ri[1], rf[1]);
            //if (ps[2] == 110125033) printf("n=%u %llu %llu\n", n, ri[2], rf[2]);
            //if (ps[3] == 110125033) printf("n=%u %llu %llu\n", n, ri[3], rf[3]);
            
            for (uint32_t idx=0; idx<4; idx++)
            {
               if (rf[idx] == pOne[idx])
                  if (ip_MultiFactorialApp->ReportFactor(ps[idx], n, -1))
                     VerifyFactor(true, ps[idx], n, -1);
                  
               if (rf[idx] == mOne[idx])
                  if (ip_MultiFactorialApp->ReportFactor(ps[idx], n, +1))
                     VerifyFactor(true, ps[idx], n, +1);
            }
         }
      }
            
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] >= maxPrime)
         break;
   }
}
void  MultiFactorialWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("MultiFactorialWorker::TestMiniPrimeChunk not implemented");
}

bool  MultiFactorialWorker::VerifyFactor(bool badFactorIsFatal, uint64_t p, uint32_t theN, int32_t theC)
{
   uint64_t rem = 1;
   int32_t  n = (int32_t) theN;
   bool     termIsPrime = true;
   
   fpu_push_1divp(p);
   
   while (n > 1)
   {      
      if (rem * n > p + 1)
         termIsPrime = false;
      
      rem = fpu_mulmod(rem, n, p);
      
      n -= ii_MultiFactorial;
   }
   
   fpu_pop();
      
   if (rem == +1)
   {
      if (termIsPrime)
         ip_MultiFactorialApp->ReportPrime(p, theN, -1);
   
      return true;
   }
   
   if (p == rem + 1)
   {
      if (termIsPrime)
         ip_MultiFactorialApp->ReportPrime(p, theN, +1);

      return true;
   }   

   if (badFactorIsFatal)
   {
      if (ii_MultiFactorial == 1)
         FatalError("%" PRIu64" is not a factor of %u!%+d", p, theN, theC);
      else
         FatalError("%" PRIu64" is not a factor of %u!%u%+d", p, theN, ii_MultiFactorial, theC);
   }
   
   return false;
}
