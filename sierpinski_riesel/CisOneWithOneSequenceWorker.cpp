/* CisOneWithOneSequenceWorker.cpp -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "CisOneWithOneSequenceWorker.h"
#include "CisOneSequenceHelper.h"
#include "../core/inline.h"
#include "../core/MpArith.h"

#define NO_PARITY 999

//#define IS_DEBUG
#define DEBUG_PRIME 1

#define N_TERM(q, i, j)      ((ii_SieveLow + (j) + (i)*babySteps)*ii_BestQ + q)

CisOneWithOneSequenceWorker::CisOneWithOneSequenceWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);
   
   ip_CisOneHelper = (CisOneSequenceHelper *) appHelper;
   
   
   // Everything we need is done in the constuctor of the parent class
   ib_Initialized = true;
}

void  CisOneWithOneSequenceWorker::CleanUp(void)
{
   delete ip_HashTable;
   
   xfree(resBD);
   xfree(resBJ);
   xfree(resX);
}

void  CisOneWithOneSequenceWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
{
   ii_BestQ = bestQ;   
   ii_SieveLow = ii_MinN / ii_BestQ;

   resX = (MpRes *) xmalloc((POWER_RESIDUE_LCM+5) * sizeof(MpRes));
   resBD = (MpRes *) xmalloc((ii_SubsequenceCount+4)*sizeof(MpRes));
   resBJ = (MpRes *) xmalloc((ii_SubsequenceCount+4)*sizeof(MpRes));
   
   ip_Legendre = ip_FirstSequence->legendrePtr;
   
   ib_UseLegendreTables = ip_CisOneHelper->UseLegendreTables();
   ip_DivisorShifts = ip_CisOneHelper->GetDivisorShifts();
   ii_PrlCount = ip_CisOneHelper->GetPrlCount();
   ip_PrlIndices = ip_CisOneHelper->GetPrlIndices();
   
   ip_HashTable = new HashTable(ip_CisOneHelper->GetMaxBabySteps());
}

void  CisOneWithOneSequenceWorker::TestMegaPrimeChunk(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t p;
   uint64_t invBase;
   uint32_t k, orderOfB, ssCount;
   uint32_t babySteps, giantSteps;
   uint32_t i, j;
   sp_t     parity;

   vector<uint64_t>::iterator it = iv_Primes.begin();
   
   while (it != iv_Primes.end())
   {
      p = *it;
      it++;

      parity = GetParity(p);

      if (parity != SP_NO_PARITY)
      {
         MpArith mp(p);

         // compute 1/base (mod p)
         invBase = invmod64(ii_Base, p);

         MpRes resB = mp.nToRes(ii_Base);
         MpRes resInvBase = mp.nToRes(invBase);
      
         ssCount = SetupDiscreteLog(mp, resB, resInvBase, parity);
   
         if (ssCount > 0)
         {
            ip_HashTable->Clear();
   
            babySteps = ip_Subsequences[ssCount-1].babySteps;
            giantSteps = ip_Subsequences[ssCount-1].giantSteps;
 
            orderOfB = BabySteps(mp, resB, resInvBase, babySteps);
            
            if (orderOfB > 0)
            {
               // If orderOfB > 0, then this is all the information we need to
               // determine every solution for this p, so no giant steps are neede
               for (k=0; k<ssCount; k++)
               {
                   j = ip_HashTable->Lookup(resBD[k]);

                   while (j < babySteps * giantSteps)
                   {
                      ip_SierpinskiRieselApp->ReportFactor(p, ip_FirstSequence, N_TERM(ip_Qs[k], 0, j), true);
                      
                      j += orderOfB;
                   }
               }
            }
            else
            {
               // First giant step
               for (k=0; k<ssCount; k++)
               {
                  j = ip_HashTable->Lookup(resBD[k]);

                  if (j != HASH_NOT_FOUND)
                     ip_SierpinskiRieselApp->ReportFactor(p, ip_FirstSequence, N_TERM(ip_Qs[k], 0, j), true);
               }
  
               // Remaining giant steps
               if (giantSteps > 1)
               {
                  MpRes resBQM = mp.pow(resBexpQ, babySteps);
                  
                  for (i=1; i<giantSteps; i++)
                  {
                     for (k=0; k<ssCount; k++)
                     {
                        resBD[k] = mp.mul(resBD[k], resBQM);
                        
                        j = ip_HashTable->Lookup(resBD[k]);

                        if (j != HASH_NOT_FOUND)
                           ip_SierpinskiRieselApp->ReportFactor(p, ip_FirstSequence, N_TERM(ip_Qs[k], i, j), true);
                     }
                  }
               }
            }
         }
      }
      
      SetLargestPrimeTested(p, 1);
      
      if (p >= maxPrime)
         return;
   }
}
void  CisOneWithOneSequenceWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CisOneWithOneSequenceWorker::TestMiniPrimeChunk not implemented");
}

sp_t   CisOneWithOneSequenceWorker::GetParity(uint64_t p)
{
   // Mixed parity sequences
   if (ip_FirstSequence->parity == SP_MIXED)
   {
      bool qr_m1, qr_p1;
      
      if (ib_UseLegendreTables)
      {
         uint32_t qr_mod = (p/2) % ip_Legendre->mod;
                  
         qr_m1 = (ip_Legendre->dualParityMapM1[L_BYTE(qr_mod)] & L_BIT(qr_mod));
         qr_p1 = (ip_Legendre->dualParityMapP1[L_BYTE(qr_mod)] & L_BIT(qr_mod));
      }
      else
      {
         int32_t sym = legendre(ip_FirstSequence->kcCore, p);
         
         qr_m1 = (sym == 1);
         qr_p1 = (sym == legendre(ii_Base, p));
      }
      
      if (qr_m1)
         return (qr_p1 ? SP_MIXED : SP_EVEN);

      return (qr_p1 ? SP_ODD : SP_NO_PARITY);
   }
   
   // Single parity sequences
   bool qr;
   
   if (ib_UseLegendreTables)
   {
      uint32_t qr_mod = (p/2) % ip_Legendre->mod;
         
      qr = (ip_Legendre->oneParityMap[L_BYTE(qr_mod)] & L_BIT(qr_mod));
   }
   else
   {
      int32_t sym = legendre(ip_FirstSequence->kcCore, p);
      
      if (ip_FirstSequence->parity == SP_EVEN)
         qr = (sym == 1);
      else
         qr = (sym == legendre(ii_Base, p));
   }
   
   if (qr)
      return ip_FirstSequence->parity;
      
   return SP_NO_PARITY;
}

// This function builds the list ii_CSSList[] of subsequences (k*b^d)*(b^Q)^m+c for
// which p may be a factor (-ckb^d is a quadratic/cubic/quartic/quintic
// residue with respect to p) and initialises the table D64[] with the
// values -c/(k*b^d) (mod p). As a side effect, bQ is set to the value b^Q
// (mod p) for use later in bsgs64(). Returns the number of subsequences listed in ii_CSS[].
uint32_t  CisOneWithOneSequenceWorker::SetupDiscreteLog(MpArith mp, MpRes resBase, uint64_t resInvBase, sp_t parity)
{
   uint64_t   negCK, pShift, p = mp.p();
   uint32_t   idx;
   uint32_t   h, r;
   int16_t    shift;
   MpRes      resNegCK;

   /* neg_ck <-- -k/c (mod p) == -ck (mod p) */
   if (p < ip_FirstSequence->k)
      negCK = ip_FirstSequence->k % p;
   else 
      negCK = ip_FirstSequence->k;
  
   if (ip_FirstSequence->c > 0)
      negCK = p - negCK;
   
   idx = (p/2) % (POWER_RESIDUE_LCM/2);
   shift = ip_DivisorShifts[idx];
   
   resNegCK = mp.nToRes(negCK);
   
   if (shift == 0)
   {
      r = ip_PrlIndices[1];
      h = 0;
   
      // p = 1 (mod 2) is all we know, check for quadratic residues only.
      idx = ip_FirstSequence->congruentQIndices[CSS_INDEX(parity, r, h)];
   
      if (idx == 0)
         return 0;
         
      ip_Qs = &ip_FirstSequence->congruentQs[idx];

      idx = ip_FirstSequence->congruentLadderIndices[CSS_INDEX(parity, r, h)];
      ip_Ladders = &ip_FirstSequence->congruentLadders[idx];
            
      // For each subsequence (k*b^d)*(b^Q)^(n/Q)+c, compute -c/(k*b^d) (mod p)
      return BuildHashTableAndClimbLadder(mp, resBase, resNegCK);
   }
   
   if (shift > 0)
   {
      // p = 1 (mod s), where s is not a power of 2. Check for r-th power
      // residues for each prime power divisor r of s.
      pShift = p / shift;
   }
   else
   {
      // p = 1 (mod 2^s), where s > 1. Check for r-th power residues for each divisor
      // r of s. We handle this case seperately to avoid computing p/s using plain division.
      pShift = p >> (-shift);
      shift = 1 << (-shift);
   }
     
   resX[0] = mp.one();
   resX[1] = mp.pow(resInvBase, pShift);
  
   for (r=1; resX[r] != resX[0]; r++)
      resX[r+1] = mp.mul(resX[r], resX[1]);

   if (shift % r != 0)
      FatalError("SetupDiscreteLog issue, shift %% xIdx != 0 (%u %% %u = %u)", shift, r, shift % r);
    
   // resX[r] <- (-ck)^((p-1)/shift)
   resX[r] = mp.pow(resNegCK, pShift);

   // Find h such that resX[h] = resX[r], i.e. (-ckb^h)^((p-1)/r)=1 (mod p), or h=r if not found
   for (h=0; resX[r] != resX[h]; h++)
      ;

   // If no h was found, then there is nothing further to do.
   if (h == r)
      return 0;

   r = ip_PrlIndices[r];
      
   idx = ip_FirstSequence->congruentQIndices[CSS_INDEX(parity, r, h)];
   
   if (idx == 0)
      return 0;
         
   ip_Qs = &ip_FirstSequence->congruentQs[idx];
   
   idx = ip_FirstSequence->congruentLadderIndices[CSS_INDEX(parity, r, h)];
   ip_Ladders = &ip_FirstSequence->congruentLadders[idx];
      
   // -ckb^d is an r-th power residue for at least one term (k*b^d)*(b^Q)^(n/Q)+c of this subsequence
   return BuildHashTableAndClimbLadder(mp, resBase, resNegCK);
}

// Assign BJ64[i] = b^i (mod p) for each i in the ladder.
// Return b^Q (mod p).
uint32_t  CisOneWithOneSequenceWorker::BuildHashTableAndClimbLadder(MpArith mp, MpRes resBase, MpRes resNegCK)
{
   uint32_t  i, j, idx, lLen, qLen;

   lLen = *ip_Ladders;
   
   // Precompute b^d (mod p) for 0 <= d <= Q, as necessary
   resX[0] = mp.one();
   resX[1] = resBase;
   resX[2] = mp.mul(resX[1], resX[1]);
   
   i = 2;
   for (j=0; j<lLen; j++)
   {
      idx = ip_Ladders[j+1];
      
      resX[i+idx] = mp.mul(resX[i], resX[idx]);
      
      i += idx;
   }

   qLen = *ip_Qs;
   ip_Qs++;
   
   resBexpQ = resX[ii_BestQ];

   for (j=0; j<qLen; j++)
      resBD[j] = mp.mul(resX[ip_Qs[j]], resNegCK);

   return qLen;
}

uint32_t  CisOneWithOneSequenceWorker::BabySteps(MpArith mp, MpRes resBase, MpRes resInvBase, uint32_t babySteps)
{
   MpRes    resInvBaseExpQ;
   MpRes    resBJ;
   MpRes    firstResBJ;
   uint32_t j;
      
   // b <- inv_b^Q (mod p)
   resInvBaseExpQ = mp.pow(resInvBase, ii_BestQ);
  
   firstResBJ = resBJ = mp.pow(resInvBaseExpQ, ii_SieveLow);

   for (j=0; j<babySteps; j++)
   { 
      ip_HashTable->Insert(resBJ, j);
      
      resBJ = mp.mul(resBJ, resInvBaseExpQ);
            
      if (resBJ == firstResBJ)
         return j + 1;
   }
   
   return 0;
}
