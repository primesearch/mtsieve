/* MultiFactorialWorker.cpp -- (C) Mark Rodenkirch, October 2012

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <math.h>
#include "MultiFactorialWorker.h"
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
   uint32_t  n;
   
   // if i <= n_pair then (i - 1) * i < p. Compute n! = (2 * 3) * (4 * 5) * ... * ((n - 1) * n)
   uint32_t  n_pair = std::max(2u, std::min(ii_MinN, uint32_t(sqrt(double(il_PrimeList[0])))) & ~1u);

   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
      MpArithVec mp(ps);

      const MpResVec pOne = mp.one();
      const MpResVec mOne = mp.sub(mp.zero(), pOne);
      const MpResVec two = mp.add(pOne, pOne);
      const MpResVec four = mp.add(two, two);
      const MpResVec eight = mp.add(four, four);

      // ri = residue of i, rf = residue of i!
      MpResVec ri = pOne, rf = pOne;
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
      ri = mp.nToRes(n_pair - 1);
      for (n = n_pair; n < ii_MinN; ++n)
      {
         ri = mp.add(ri, pOne);
         rf = mp.mul(rf, ri);
      }

      // Factorial and check if i! = +/-1
      for (; n <= ii_MaxN; ++n)
      {
         ri = mp.add(ri, pOne);
         rf = mp.mul(rf, ri);

         if (MpArithVec::at_least_one_is_equal(rf, pOne) || MpArithVec::at_least_one_is_equal(rf, mOne))
         {
            for (size_t k = 0; k < VECTOR_SIZE; ++k)
            {
               if (rf[k] == pOne[k])
                  ip_MultiFactorialApp->ReportFactor(ps[k], n, -1);
                  
               if (rf[k] == mOne[k]) 
                  ip_MultiFactorialApp->ReportFactor(ps[k], n, +1);
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
   uint64_t  ps[4], maxPrime = ip_App->GetMaxPrime();
   uint32_t  maxNFirstLoop = ii_MinN - ii_MultiFactorial;
   uint32_t  n, startN;
   
   uint32_t  pIdx = 0;
   
   while (pIdx < ii_PrimesInList)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];

      pIdx += 4;

      MpArithVec mp(ps);

      const MpResVec pOne = mp.one();
      const MpResVec mOne = mp.sub(mp.zero(), pOne);
      const MpResVec resMf = mp.nToRes(ii_MultiFactorial);

      for (startN=1; startN<=ii_MultiFactorial; startN++)
      {
         // If startN is odd and mf is even, then i!mf is always odd, thus
         // startN!mf+1 and startN!mf-1 are always even.
         if (!(ii_MultiFactorial & 1) && (startN & 1))
            continue;

         MpResVec ri = mp.nToRes(startN);
         MpResVec rf = ri;
         
         // At this time we have:
         //    ri = residual of startN (mod p)
         //    rf = residual of startN!mf (mod p)
         
         n = startN + ii_MultiFactorial;
         for (; n<maxNFirstLoop; n+=ii_MultiFactorial)
         {
            ri = mp.add(ri, resMf);
            rf = mp.mul(rf, ri);
         }
         
         // At this time we have:
         //    ri = residual of mn (mod p)
         //    rf = residual of mn!mf (mod p)
         // where mn is the max n < ii_MinN for this starting n
         for (; n <=ii_MaxN; n+=ii_MultiFactorial)
         {
            ri = mp.add(ri, resMf);
            rf = mp.mul(rf, ri);

            if (MpArithVec::at_least_one_is_equal(rf, pOne) || MpArithVec::at_least_one_is_equal(rf, mOne))
            {
               for (size_t k = 0; k < VECTOR_SIZE; ++k)
               {
                  if (rf[k] == pOne[k])
                     ip_MultiFactorialApp->ReportFactor(ps[k], n, -1);
                     
                  if (rf[k] == mOne[k]) 
                     ip_MultiFactorialApp->ReportFactor(ps[k], n, +1);
               }
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
