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

   ii_MinX = ip_LifchitzApp->GetMinX();
   ii_MaxX = ip_LifchitzApp->GetMaxX();
   ii_MinY = ip_LifchitzApp->GetMinY();
   ii_MaxY = ip_LifchitzApp->GetMaxY();
   
   ii_MinBase = (ii_MinX > ii_MinY ? ii_MinY : ii_MinX);
   ii_MaxBase = (ii_MaxX > ii_MaxY ? ii_MaxX : ii_MaxY);

   ip_Terms = NULL;
   il_NextTermsBuild = 0;
   iv_Bases.resize(1 + ii_MaxBase - ii_MinBase);
   std::fill(iv_Bases.begin(), iv_Bases.end(), false);
   
   ii_Elements = 1 + ii_MaxX - ii_MinX;
   
   if (ii_Elements < 8)
      ii_Elements = 8;
   
   ii_HashSize = (1 << 11);
   
   while (ii_HashSize < ii_Elements/0.65)
      ii_HashSize <<= 1;
   
   ii_HashM1 = ii_HashSize - 1;
   
   ip_Remainders = (MpResVec *) xmalloc((1 + ii_MaxBase - ii_MinBase), sizeof(MpResVec), "remainders");
   ip_Hashes = (hash_t *) xmalloc(2 * ii_HashSize, sizeof(hash_t), "hashTable");
   ip_Residues = (MpRes *) xmalloc(ii_Elements, sizeof(MpRes), "hashElements");
   
   ib_Initialized = true;
}

void  LifchitzWorker::CleanUp(void)
{
   xfree(ip_Hashes);
   xfree(ip_Residues);
   
   if (ip_Terms != NULL)
      xfree(ip_Terms);
}

void  LifchitzWorker::TestMegaPrimeChunk(void)
{
   if (il_PrimeList[0] <= ii_MaxX)
      TestInitialChunk();
   else
      TestLaterChunk();
}

void  LifchitzWorker::TestInitialChunk(void)
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

      for (uint64_t idx=0; ip_Terms[idx].x>0; idx++)
      {
         for (uint32_t pIdx=0; pIdx<4; pIdx++)
         {
            uint64_t xPowX = ip_Remainders[BIT(ip_Terms[idx].x)][pIdx];
            uint64_t yPowY = ip_Remainders[BIT(ip_Terms[idx].y)][pIdx];
            
            if (ip_Terms[idx].signs & M_ONE)
            {
               if (xPowX == yPowY)
                  ip_LifchitzApp->ReportFactor(ps[pIdx], ip_Terms[idx].x, ip_Terms[idx].y, -1, idx);
            }
            
            if (ip_Terms[idx].signs & P_ONE)
            {
               if (xPowX + yPowY == ps[pIdx])
                  ip_LifchitzApp->ReportFactor(ps[pIdx], ip_Terms[idx].x, ip_Terms[idx].y, +1, idx);
            }
         }
      }
      
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] > maxPrime)
         break;

      // Stop ASAP if the user hit ^C since each execution of this method can take a long time to complete.
      if (ip_LifchitzApp->IsInterrupted())
         break;
   }
}

void  LifchitzWorker::TestLaterChunk(void)
{
   uint64_t thePrime;
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t x, y;
   uint32_t slot, emptySlot;
   
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx++)
   {
      memset(ip_Hashes, 0, ii_HashSize * sizeof(hash_t));
      memset(ip_Residues, 0, ii_Elements * sizeof(MpRes));
      
      thePrime = il_PrimeList[pIdx+0];
     
      MpArith  mp(thePrime);
      MpRes    resBase, resPow;

      emptySlot = ii_HashSize;
      
      for (uint32_t x=ii_MinX; x<=ii_MaxX; x++)
      {
         resBase = mp.nToRes(x);
         resPow = mp.pow(resBase, x);

         ip_Residues[x-ii_MinX] = resPow;
         
         slot = resPow & ii_HashM1;

         if (ip_Hashes[slot].x == 0)
         {
            ip_Hashes[slot].x = x;
            ip_Hashes[slot].next = 0;
         }
         else
         {
            // Create a change of x where we insert this x into the chain for this residue
            ip_Hashes[emptySlot].x = ip_Hashes[slot].x;
            ip_Hashes[emptySlot].next = ip_Hashes[slot].next;
            
            ip_Hashes[slot].x = x;
            ip_Hashes[slot].next = emptySlot;
            
            emptySlot++;
         } 
      }
      
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         resBase = mp.nToRes(y);
         resPow = mp.pow(resBase, y);
         
         slot = resPow & ii_HashM1;

         // We need to meet the following criteria to know we have a match:
         //    ip_Hashes[slot].x is not 0
         //    x > y for that entry
         //    ip_Residues[x] = resPow
         while (ip_Hashes[slot].x > y)
         {
            x = ip_Hashes[slot].x;
            
            if (ip_Residues[x-ii_MinX] == resPow)
               ip_LifchitzApp->ReportFactor(thePrime, x, y, -1);
            
            slot = ip_Hashes[slot].next;
            
            if (slot == 0)
               break;
         }

         resPow = mp.sub(0, resPow);

         slot = resPow & ii_HashM1;
         
         while (ip_Hashes[slot].x > y)
         {
            x = ip_Hashes[slot].x;
            
            if (ip_Residues[x-ii_MinX] == resPow)
               ip_LifchitzApp->ReportFactor(thePrime, x, y, +1);
            
            slot = ip_Hashes[slot].next;
            
            if (slot == 0)
               break;
         }
      }

      SetLargestPrimeTested(thePrime, 1);
      
      if (thePrime > maxPrime)
         break;

      // Stop ASAP if the user hit ^C since each execution of this method can take a long time to complete.
      if (ip_LifchitzApp->IsInterrupted())
         break;
   }
}

void  LifchitzWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("LifchitzWorker::TestMiniPrimeChunk not implemented");
}

void   LifchitzWorker::BuildTerms(void)
{ 
   if (ip_Terms != NULL)
      xfree(ip_Terms);

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
