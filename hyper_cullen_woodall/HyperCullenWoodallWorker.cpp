/* HyperCullenWoodallWorker.cpp -- (C) Mark Rodenkirch

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "HyperCullenWoodallWorker.h"

#define X_INDEX(x)  ((x) - ii_MinX)
#define Y_INDEX(y)  ((y) - ii_MinY)

#define MAX_POWERS   10

HyperCullenWoodallWorker::HyperCullenWoodallWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{ 
   ip_HyperCullenWoodallApp = (HyperCullenWoodallApp *) theApp;
   
   ib_IsPlus = ip_HyperCullenWoodallApp->IsPlus();
   ib_IsMinus = ip_HyperCullenWoodallApp->IsMinus();
   
   ii_MinX = ip_HyperCullenWoodallApp->GetMinX();
   ii_MaxX = ip_HyperCullenWoodallApp->GetMaxX();
   ii_MinY = ip_HyperCullenWoodallApp->GetMinY();
   ii_MaxY = ip_HyperCullenWoodallApp->GetMaxY();

   ii_XCount = ip_HyperCullenWoodallApp->GetXCount();
   ii_YCount = ip_HyperCullenWoodallApp->GetYCount();
   
   il_NextTermsBuild = 0;
  
   ib_Initialized = true;
}

void  HyperCullenWoodallWorker::CleanUp(void)
{
}

void  HyperCullenWoodallWorker::TestMegaPrimeChunk(void)
{
   uint64_t ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t idx, x, y;
   MpResVec xyRes, yxRes;
   
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
      // Every once in a while rebuild the term lists as it will have fewer entries
      // which will speed up testing for the next range of p.
      if (ps[0] > il_NextTermsBuild)
      {
         ip_HyperCullenWoodallApp->BuildTerms();
         
         ip_xTerms = ip_HyperCullenWoodallApp->GetXTerms();
         ip_yTerms = ip_HyperCullenWoodallApp->GetYTerms();
               
         il_NextTermsBuild = (ps[3] << 2);
      }
      
      MpArithVec mp(ps);
      MpResVec   mpOne = mp.one();
      MpResVec   mpPm1 = mp.sub(mp.zero(), mp.one());
      
      // Compute x^y for all x and y
      ComputeResidues(mp, ip_xTerms, ii_MinY);
      
      // Compute y^x for all x and y
      ComputeResidues(mp, ip_yTerms, ii_MinX);
      
      base_t *xPtr = ip_xTerms;
      base_t *yPtr;

      while (xPtr->base > 0)
      {
         x = xPtr->base;
         
         for (idx=0; idx<xPtr->powerCount; idx++)
         {
            y = xPtr->powers[idx];
            xyRes = xPtr->residues[idx];
            
            yPtr = &ip_yTerms[y-ii_MinY];

            yxRes = yPtr->residues[x-ii_MinX];
            
            MpResVec res = mp.mul(xyRes, yxRes);
            
            for (size_t k = 0; k < VECTOR_SIZE; ++k)
            {
               if (res[k] == mpOne[k])
                  ip_HyperCullenWoodallApp->ReportFactor(ps[k], x, y, -1);
                  
               if (res[k] == mpPm1[k]) 
                  ip_HyperCullenWoodallApp->ReportFactor(ps[k], x, y, +1);
            }
         }
         
         xPtr++;
      }

      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] > maxPrime)
         break;
   }
}

void  HyperCullenWoodallWorker::ComputeResidues(MpArithVec mp, base_t *terms, uint32_t minPower)
{  
   base_t   *bPtr = terms;
   MpResVec  mpPowers[MAX_POWERS+1];
   MpResVec  mpBase;
   MpResVec  mpRes;
   uint32_t  idx, prevPower, powerDiff;
   
   while (bPtr->base > 0)
   {
      mpBase = mp.nToRes(bPtr->base);
      
      mpPowers[0] = mp.one();
      mpPowers[1] = mpBase;
      
      // Multiply successive terms by the base (mod p)
      for (idx=2; idx<=MAX_POWERS; idx++)
          mpPowers[idx] = mp.mul(mpPowers[idx-1], mpBase);
      
      mpRes = mp.pow(mpBase, bPtr->powers[0]);
      
      if (bPtr->indexedByPower)
         bPtr->residues[bPtr->powers[0] - minPower] = mpRes;
      else
         bPtr->residues[0] = mpRes;

      prevPower = bPtr->powers[0];

      for (idx=1; idx<bPtr->powerCount; idx++)
      {
         powerDiff = bPtr->powers[idx] - prevPower;
         
         while (powerDiff > MAX_POWERS) 
         {
            mpRes = mp.mul(mpRes, mpPowers[MAX_POWERS]);
            powerDiff -= MAX_POWERS;
         };
         
         if (powerDiff > 0)
            mpRes = mp.mul(mpRes, mpPowers[powerDiff]);
            
         if (bPtr->indexedByPower)
            bPtr->residues[bPtr->powers[idx] - minPower] = mpRes;
         else
            bPtr->residues[idx] = mpRes;
      
         prevPower = bPtr->powers[idx];
      }

      bPtr++;
   }
}

void  HyperCullenWoodallWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("TestMiniPrimeChunk has not been implemented");
}

