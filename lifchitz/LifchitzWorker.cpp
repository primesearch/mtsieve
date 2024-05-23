/* LifchitzWorker.cpp -- (C) Mark Rodenkirch, April 2024

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "LifchitzWorker.h"
#include "../core/MpArith.h"

LifchitzWorker::LifchitzWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{  
   ip_LifchitzApp = (LifchitzApp *) theApp;

   ii_MinX = ip_LifchitzApp->GetMinX();
   ii_MaxX = ip_LifchitzApp->GetMaxX();
   ii_MinY = ip_LifchitzApp->GetMinY();
   ii_MaxY = ip_LifchitzApp->GetMaxY();
   
   ii_Elements = 1 + ii_MaxX - ii_MinX;
   
   if (ii_Elements < 8)
      ii_Elements = 8;
   
   ii_HashSize = (1 << 11);
   
   while (ii_HashSize < ii_Elements/0.65)
      ii_HashSize <<= 1;
   
   ii_HashM1 = ii_HashSize - 1;
      
   ip_Hashes = (hash_t *) xmalloc(2 * ii_HashSize, sizeof(hash_t), "hashTable");
   ip_Residues = (MpRes *) xmalloc(ii_Elements, sizeof(MpRes), "hashElements");
   
   ib_Initialized = true;
}

void  LifchitzWorker::CleanUp(void)
{
   xfree(ip_Hashes);
   xfree(ip_Residues);
}

void  LifchitzWorker::TestMegaPrimeChunk(void)
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
