/* KBBWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "KBBWorker.h"
#include "../core/MpArithVector.h"

KBBWorker::KBBWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_KBBApp = (KBBApp *) theApp;
   
   il_K = ip_KBBApp->GetK();

   // Give a little extra space
   ii_BaseCount = ip_KBBApp->GetMaxB() - ip_KBBApp->GetMinB() + 10;
   
   ip_Bases = (uint32_t *) xmalloc(ii_BaseCount * sizeof(uint32_t));
   
   il_NextBaseBuild = 0;
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  KBBWorker::CleanUp(void)
{
}

void  KBBWorker::TestMegaPrimeChunk(void)
{
   uint64_t ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t idx;

   for (uint32_t pIdx=0; pIdx<ii_WorkSize; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
      // Every once in a while rebuild the list of base as it will have fewer entries
      // which will speed up testing for the next range of p.
      if (ps[0] > il_NextBaseBuild)
      {
         memset(ip_Bases, 0, ii_BaseCount * sizeof(uint32_t));
         
         ip_KBBApp->GetBases(ip_Bases);
         
         il_NextBaseBuild = (ps[3] << 1);
      }

      MpArithVec mp(ps);

      const MpResVec resK = mp.nToRes(il_K);
      
      const MpResVec pOne = mp.one();
      const MpResVec mOne = mp.sub(mp.zero(), pOne);
   
      for (idx=0; idx<ii_BaseCount; idx++)
      {
         if (ip_Bases[idx] == 0)
            break;
         
         MpResVec res = mp.nToRes(ip_Bases[idx]);
         
         res = mp.pow(res, ip_Bases[idx]);
         
         res = mp.mul(res, resK);
                  
         if (MpArithVec::at_least_one_is_equal(res, pOne))
         {
            for (size_t k = 0; k < VECTOR_SIZE; ++k)
            {
               if (res[k] == pOne[k])
                  ip_KBBApp->ReportFactor(ps[k], ip_Bases[idx], -1);
            }
         }
         
         if (MpArithVec::at_least_one_is_equal(res, mOne))
         {
            for (size_t k = 0; k < VECTOR_SIZE; ++k)
            {
               if (res[k] == mOne[k]) 
                  ip_KBBApp->ReportFactor(ps[k], ip_Bases[idx], +1);
            }
         }
      }
      
      if (ps[3] > maxPrime)
         break;
   }
   
   SetLargestPrimeTested(ps[3], 4);
}

void  KBBWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("KBBWorker::TestMiniPrimeChunk not implemented");
}
