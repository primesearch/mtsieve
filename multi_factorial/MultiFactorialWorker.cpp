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
   ip_Terms = ip_MultiFactorialApp->GetTerms();

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
   
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
      MpArithVec mp(ps);

      const MpResVec pOne = mp.one();
      const MpResVec mOne = mp.sub(mp.zero(), pOne);

      MpResVec resRem = pOne;
      MpResVec resBase = pOne;
      MpResVec resTemp = pOne;
      uint32_t power = 0;
      uint32_t tIdx = 0;
      
      while (ip_Terms[0].power[tIdx] > 0)
      {
         resBase = mp.nToRes(ip_Terms[0].base[tIdx]);
            
         // If this base has the same power as the previous base, just muliply
         // We will do exponentiation before we multiply by resRem
         if (ip_Terms[0].power[tIdx] == power)
         {
            resTemp = mp.mul(resTemp, resBase);
            tIdx++;
            continue;
         }

         if (power != 0)
         {
            // resRem = resTemp^power * resRem
            resTemp = mp.pow(resTemp, power);
            resRem = mp.mul(resRem, resTemp);
         }
         
         power = ip_Terms[0].power[tIdx];
         resTemp = resBase;
         tIdx++;
      }

      if (power != 0)
      {
         if (power > 1)
            resTemp = mp.pow(resTemp, power);

         resRem = mp.mul(resRem, resTemp);
      }
      
      n = ii_MinN - 1;
      MpResVec resN = mp.nToRes(n);
      
      // At this point resRem = (n-1)! and resN = (n-1)
      while (n < ii_MaxN)
      {
         n++;
         
         resN = mp.add(resN, pOne);
         resRem = mp.mul(resRem, resN);

         if (MpArithVec::at_least_one_is_equal(resRem, pOne, mOne))
         {
            for (size_t k = 0; k < VECTOR_SIZE; ++k)
            {
               if (resRem[k] == pOne[k])
                  ip_MultiFactorialApp->ReportFactor(ps[k], n, -1);
                  
               if (resRem[k] == mOne[k]) 
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
   uint32_t  n;
   
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
      const MpResVec resAdd = mp.nToRes(ii_MultiFactorial);

      // If ii_Multifactorial == 2 then mf=0 = 2*4*6*... and mf=1 = 1*3*5*...
      // If ii_Multifactorial == 3 then mf=0 = 3*6*9*... and mf=1 = 4*7*10*... and mf=2 = 5*8*11*...
      for (uint32_t mf=0; mf<ii_MultiFactorial; mf++)
      {
         // If ii_Multifactorial is even and mf is odd then 
         // n!ii_Multifactorial+1 and n!ii_Multifactorial-1 are always even
         // when mf is odd, so we do not need to go any further.
         if (!(ii_MultiFactorial & 1) && (mf & 1))
            continue;

         MpResVec resRem = pOne;
         MpResVec resBase = pOne;
         MpResVec resTemp = pOne;
         uint32_t power = 0;
         uint32_t tIdx = 0;
         
         while (ip_Terms[mf].power[tIdx] > 0)
         {
            resBase = mp.nToRes(ip_Terms[mf].base[tIdx]);
               
            // If this base has the same power as the previous base, just muliply
            // We will do exponentiation before we multiply by resRem
            if (ip_Terms[mf].power[tIdx] == power)
            {
               resTemp = mp.mul(resTemp, resBase);
               tIdx++;
               continue;
            }

            if (power != 0)
            {
               // resRem = resTemp^power * resRem
               resTemp = mp.pow(resTemp, power);
               resRem = mp.mul(resRem, resTemp);
            }
            
            power = ip_Terms[mf].power[tIdx];
            resTemp = resBase;
            tIdx++;
         }

         if (power != 0)
         {
            if (power > 1)
               resTemp = mp.pow(resTemp, power);

            resRem = mp.mul(resRem, resTemp);
         }
         
         n = ii_MinN - 1;
         while (n % ii_MultiFactorial != mf)
            n--;

         MpResVec resN = mp.nToRes(n);
         
         // At this point resRem = (n-1)! and resN = (n-1)
         // where n is the largest n less than ii_MinN for this mf.
         while (n < ii_MaxN)
         {
            n += ii_MultiFactorial;
            resN = mp.add(resN, resAdd);
            resRem = mp.mul(resRem, resN);

            if (MpArithVec::at_least_one_is_equal(resRem, pOne, mOne))
            {
               for (size_t k = 0; k < VECTOR_SIZE; ++k)
               {
                  if (resRem[k] == pOne[k])
                     ip_MultiFactorialApp->ReportFactor(ps[k], n, -1);
                     
                  if (resRem[k] == mOne[k]) 
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
