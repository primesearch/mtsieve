/* XYYXSparseWorker.cpp -- (C) Mark Rodenkirch, September 2012

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "XYYXSparseWorker.h"
#include "../core/MpArithVector.h"

XYYXSparseWorker::XYYXSparseWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{  
   ip_XYYXApp = (XYYXApp *) theApp;
   
   ib_IsPlus = ip_XYYXApp->IsPlus();
   ib_IsMinus = ip_XYYXApp->IsMinus();
   
   ip_Terms = 0;
   il_NextTermsBuild = 0;
   
   ib_Initialized = true;
}

void  XYYXSparseWorker::CleanUp(void)
{
   if (ip_Terms)
      xfree(ip_Terms);
}

void  XYYXSparseWorker::TestMegaPrimeChunk(void)
{
   uint64_t ps[4];
   
   uint64_t maxPrime = ip_App->GetMaxPrime();
   
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
            // Every once in a while rebuild the term lists as it will have fewer entries
      // which will speed up testing for the next range of p.
      // Unlike the GPU, we will put all terms into one group.
      if (ps[0] > il_NextTermsBuild)
      {
         if (ip_Terms)
            xfree(ip_Terms);
         
         ip_Terms = ip_XYYXApp->GetSparseTerms();
         
         il_NextTermsBuild = (ps[3] << 1);
      }
      
      MpArithVec mp(ps);
      MpResVec resBaseX;
      MpResVec resBaseY;
      MpResVec resXexpY;
      MpResVec resYexpX;

      const MpResVec zero = mp.zero();

      for (uint32_t idx=0; ; idx++)
      {
         if (ip_Terms[idx].x == 0)
            break;

         if (ip_Terms[idx].haveFactor)
            continue;
         
         resBaseX = mp.nToRes(ip_Terms[idx].x);
         resBaseY = mp.nToRes(ip_Terms[idx].y);
            
         resXexpY = mp.pow(resBaseX, ip_Terms[idx].y);
         resYexpX = mp.pow(resBaseY, ip_Terms[idx].x);
         
         if (ib_IsMinus)
         {
            for (size_t k = 0; k < VECTOR_SIZE; ++k)
            {
               if (resXexpY[k] == resYexpX[k])
                  ip_XYYXApp->ReportFactor(ps[k], ip_Terms[idx].x, ip_Terms[idx].y);
            }
         }
         else
         {
            resYexpX = mp.add(resXexpY, resYexpX);
               
            for (size_t k = 0; k < VECTOR_SIZE; ++k)
            {
               if (resYexpX[k] == zero[k])
                  ip_XYYXApp->ReportFactor(ps[k], ip_Terms[idx].x, ip_Terms[idx].y);
            }
         }
      }

      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] > maxPrime)
         break;
   }
}
void  XYYXSparseWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("XYYXSparseWorker::TestMiniPrimeChunk not implemented");
}
