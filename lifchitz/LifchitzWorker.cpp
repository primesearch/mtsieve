/* LifchitzWorker.cpp -- (C) Mark Rodenkirch, April 2024

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "LifchitzWorker.h"

#define BIT(base)       (base - ii_MinBase)

LifchitzWorker::LifchitzWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{  
   ip_LifchitzApp = (LifchitzApp *) theApp;

   uint32_t minX = ip_LifchitzApp->GetMinX();
   uint32_t maxX = ip_LifchitzApp->GetMaxX();
   uint32_t minY = ip_LifchitzApp->GetMinY();
   uint32_t maxY = ip_LifchitzApp->GetMaxY();
   
   ii_MinBase = (minX > minY ? minY : minX);
   ii_MaxBase = (maxX > maxY ? maxX : maxY);

   iv_Bases.resize(1 + ii_MaxBase - ii_MinBase);
   std::fill(iv_Bases.begin(), iv_Bases.end(), false);

   ip_Remainders = (MpResVec *) xmalloc((1 + ii_MaxBase - ii_MinBase) * sizeof(MpResVec));
   
   il_NextTermsBuild = 0;
   
   ib_Initialized = true;
}

void  LifchitzWorker::CleanUp(void)
{
   xfree(ip_Remainders);
}

void  LifchitzWorker::TestMegaPrimeChunk(void)
{
   uint32_t base;
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
      if (ps[3] > il_NextTermsBuild)
      {
         BuildTerms();
      
         il_NextTermsBuild = (ps[3] << 2);
      }

      MpArithVec  mp(ps);
   
      for (base=ii_MinBase; base<=ii_MaxBase; base++)
      {
         if (iv_Bases[BIT(base)])
         {
            MpResVec resBase = mp.nToRes(base);
            MpResVec rem = mp.pow(resBase, base);
            ip_Remainders[BIT(base)] = mp.resToN(rem);
         }
      }

      uint64_t idx = 0;
      
      while (ip_Terms[idx].x > 0)
      {
         for (uint32_t pIdx=0; pIdx<4; pIdx++)
         {
            uint64_t xPowX = ip_Remainders[BIT(ip_Terms[idx].x)][pIdx];
            uint64_t yPowY = ip_Remainders[BIT(ip_Terms[idx].y)][pIdx];
            
            if (ip_Terms[idx].sign == -1)
            {
               if (xPowX == yPowY)
                  ip_LifchitzApp->ReportFactor(ps[pIdx], ip_Terms[idx].x, ip_Terms[idx].y, -1);
            }
            else {
               if (xPowX + yPowY == ps[pIdx])
                  ip_LifchitzApp->ReportFactor(ps[pIdx], ip_Terms[idx].x, ip_Terms[idx].y, +1);
            }
         }
         
         idx++;
      }
      
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] > maxPrime)
         break;
   }
}

void  LifchitzWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("LifchitzWorker::TestMiniPrimeChunk not implemented");
}

void   LifchitzWorker::BuildTerms(void)
{
   ip_Terms = ip_LifchitzApp->GetTerms();
   uint64_t idx = 0;
   
   std::fill(iv_Bases.begin(), iv_Bases.end(), false);
   
   while (ip_Terms[idx].x > 0)
   {
      iv_Bases[BIT(ip_Terms[idx].x)] = true;
      iv_Bases[BIT(ip_Terms[idx].y)] = true;
      
      idx++;
   }
}
