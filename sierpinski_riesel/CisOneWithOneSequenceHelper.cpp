/* CisOneWithOneSequenceHelper.cpp -- (C) Mark Rodenkirch, May 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <assert.h>
#include "../core/inline.h"
#include "SierpinskiRieselApp.h"
#include "CisOneWithOneSequenceHelper.h"
#include "CisOneWithOneSequenceWorker.h"
#ifdef HAVE_GPU_WORKERS
#include "CisOneWithOneSequenceGpuWorker.h"
#endif

CisOneWithOneSequenceHelper::CisOneWithOneSequenceHelper(App *theApp, uint64_t largestPrimeTested) : CisOneSequenceHelper(theApp, largestPrimeTested)
{
   // The parent's constuctor does much of the work
  
   theApp->WriteToConsole(COT_OTHER, "Sieving one sequence where abs(c) = 1 for p >= %llu", largestPrimeTested);
}

Worker  *CisOneWithOneSequenceHelper::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   AbstractWorker *theWorker;

   // Note that GenericWorker inherits from Worker.  This will not
   // only create the worker, but also start it.
   
#ifdef HAVE_GPU_WORKERS
   if (gpuWorker)
      theWorker = new CisOneWithOneSequenceGpuWorker(id, ip_App, this);
   else
#endif
      theWorker = new CisOneWithOneSequenceWorker(id, ip_App, this);
      
   theWorker->Prepare(largestPrimeTested, ii_BestQ);

   return theWorker;
}

uint32_t    CisOneWithOneSequenceHelper::RateQ(uint32_t Q, uint32_t s)
{
   uint32_t work, W[POWER_RESIDUE_LCM+1];
   uint32_t i;
   //uint32_t r, babySteps, giantSteps;

   assert(Q % 2 == 0);
   assert(Q % BASE_MULTIPLE == 0);
   assert(LIMIT_BASE % Q == 0);

   for (i=0; i<=POWER_RESIDUE_LCM; i++)
      W[i] = 0;

   work = 0;
   for (i = 2; i <= POWER_RESIDUE_LCM; i += 2)
   {
      if (POWER_RESIDUE_LCM % i == 0)
      {
         W[i] = (uint32_t) EstimateWork(Q, (s+i-1)/i);

         // giantSteps are very expensive compared to other loops, so this might
         // always be the best choice even if it means more memory is needed
         //uint32_t r = 1 + ii_MaxN/Q - ii_MinN/Q;
         
         //ChooseSteps(r, s, babySteps, giantSteps);
         
         //W[i] = giantSteps;
      }

      if (gcd32(i+1,POWER_RESIDUE_LCM) == 1)
         work += W[gcd32(i, POWER_RESIDUE_LCM)];
   }
   
   return work;
}

// giantSteps are expensive compared to other loops, so it should
// have a height weight than the others
#define BABY_WORK    1.0    // 1 mulmod, 1 insert
#define GIANT_WORK   3.0    // 1 mulmod, 1 lookup
#define EXP_WORK     0.3    // 1 mulmod
#define SUBSEQ_WORK  1.4    // 1 mulmod, 1 lookup (giant step 0)
                               
// Q = q, s = number of subsequences.
double    CisOneWithOneSequenceHelper::EstimateWork(uint32_t Q, uint32_t s)
{
   uint32_t babySteps, giantSteps;
   double   work;

   ChooseSteps(Q, s, babySteps, giantSteps);

   work = babySteps*BABY_WORK + s*(giantSteps-1)*GIANT_WORK + Q*EXP_WORK + s*SUBSEQ_WORK;
   
   return work;
}

// For sequences with single parity terms (parity=+/-1):
// Set seq_mod and seq_map[0] for sequence k*b^n+c so that bit (p/2)%mod
// of seq_map[0] is set if and only if
// (-ck/p)=1 for sequences with all n even,
// (-bck/p)=1 for sequences with all n odd.
// 
// For sequences with mixed parity terms (parity=0):
// Set seq_mod, seq_map[0] and seq_map[1] for sequence k*b^n+c so that
// bit (p/2)%mod of seq_map[0] is set if and only if (-ck/p)=1,
// bit (p/2)%mod of seq_map[1] is set if and only if (-bck/p)=1.
// 
// In the worst case each table for k*b^n+c could be 4*b*k bits long.
bool   CisOneWithOneSequenceHelper::BuildLegendreTableForSequence(seq_t *seq)
{
   int64_t     r;
   uint64_t    mod, ssqfb;
   legendre_t *legendre;
   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) ip_App;

   mod = seq->squareFreeK;
   ssqfb = srApp->GetSquareFreeBase();
   
   r = seq->kcCore;
   
   switch (seq->parity)
   {
      // odd n, test for (-bck/p)==1
      case SP_ODD: 
         mod *= ssqfb;
         r *= ssqfb;
         // Fall through

      // even n, test for (-ck/p)==1
      case SP_EVEN: 
         if ((r < 0 && (-r) % 4 != 3) || (r > 0 && r % 4 != 1))
            mod *= 2;
             
         if ((2*mod+1) >= INT32_MAX)
         {
            ip_App->WriteToConsole(COT_OTHER, "Cannot use Legendre tables because square-free part of k is too large");
            return false;
         }
         
        //if (cache_file_name != NULL)
        //{
        //  if (read_legendre_cache(r,m,ii_Parity,&map0,NULL) == 0)
        //    if ((map0 = make_bitmap(m,NULL)) != NULL)
        //    {
        //      build_legendre_table(r,m,ii_Parity,map0);
        //      write_legendre_cache(r,m,ii_Parity,map0,NULL);
        //    }
        //  seq_map[0] = map0;
        //}
        //else
         legendre = (legendre_t *) xmalloc(sizeof(legendre_t));

         legendre->mapSize = L_BYTES(mod);
         legendre->oneParityMap = (uint8_t *) xmalloc(legendre->mapSize * sizeof(uint8_t));
         
         BuildLegendreMap(mod, r, "single", legendre->oneParityMap);
         
         seq->legendrePtr = legendre;
         break;

      // odd and even n, test for (-ck/p)==1 and (-bck/p)==1
      default:
         mod = mod*2*ssqfb;
      
         if ((2*mod+1) >= INT32_MAX)
         {
            ip_App->WriteToConsole(COT_OTHER, "Cannot use Legendre tables because square-free part of k is too large");
            return false;
         }
         
         if (abs(r) * ssqfb >= INT32_MAX)
         {
            ip_App->WriteToConsole(COT_OTHER, "Cannot use Legendre tables because square-free part of b is too large");
            return false;
         }
         
     //if (cache_file_name != NULL)
     //{
     //  if (read_legendre_cache(r,m,0,&map0,&map1) == 0)
     //    if ((map0 = make_bitmap(m,NULL)) != NULL)
     //    {
     //      if ((map1 = make_bitmap(m,NULL)) != NULL)
     //      {
     //        build_legendre_table(r,m,1,map0);
     //        build_legendre_table(r*b,m,-1,map1);
     //        write_legendre_cache(r,m,0,map0,map1);
     //      }
     //      else
     //      {
     //        free(map0);
     //        map0 = NULL;
     //      }
     //    }
     //  seq_map[0] = map0;
     //  seq_map[1] = map1;
     //}
     //else
     
         legendre = (legendre_t *) xmalloc(sizeof(legendre_t));
         
         legendre->mapSize = L_BYTES(mod);
         legendre->dualParityMapP1 = (uint8_t *) xmalloc(legendre->mapSize * sizeof(uint8_t));
         legendre->dualParityMapM1 = (uint8_t *) xmalloc(legendre->mapSize * sizeof(uint8_t));

         BuildLegendreMap(mod, r,       "even", legendre->dualParityMapP1);
         BuildLegendreMap(mod, r*ssqfb, "odd",  legendre->dualParityMapM1);
         
         seq->legendrePtr = legendre;
         break;
   }
   
   legendre->mod = mod;
   return true;
}

void    CisOneWithOneSequenceHelper::BuildLegendreMap(uint32_t size, int64_t r, const char *which, uint8_t *legendreMap)
{
   uint32_t i;
   uint32_t reset = 0;

   for (i=0; i<size; i++)
   {
      // Provide progress when building a large table
      if (reset == 10000000)
      {
         ip_App->WriteToConsole(COT_OTHER, "Building %s parity Legendre symbol table: %.1f%% done", which, 100.0*i/size);
         reset = 0;
      }

      if (jacobi(r, 2*i+1) == 1)
         legendreMap[L_BYTE(i)] |= L_BIT(i);

      reset++;
   }
}

// Build tables sc_lists[i][r] of pointers to lists of subsequences
// whose terms k*b^m+c satisfy m = j (mod r)
void  CisOneWithOneSequenceHelper::MakeSubseqCongruenceTables(seq_t *seq)
{
   uint32_t   ssIdx, h, r, len[SP_COUNT];
   uint16_t   tempQs[SP_COUNT][POWER_RESIDUE_LCM];
   sp_t       parity = ip_FirstSequence->parity;
   
   seq->congruentIndexCount = SP_COUNT * (ii_PrlCount + 1) * (POWER_RESIDUE_LCM + 1);
   seq->congruentQIndices = (uint32_t *) xmalloc(seq->congruentIndexCount * sizeof(uint32_t));
   seq->congruentLadderIndices = (uint32_t *) xmalloc(seq->congruentIndexCount * sizeof(uint32_t));
   
   // We will not use the first element so that an index of 0 means "does not exist"
   ii_TempQIndex = ii_TempLadderIndex = 1;
   for (r=1; r<=POWER_RESIDUE_LCM; r++)
   {
      if (POWER_RESIDUE_LCM % r != 0)
         continue;

      for (h=0; h<r; h++)
      {
         len[0] = len[1] = len[2] = 0;
         
         for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
         {
            if (CongruentTerms(ssIdx, r, h))
            {
               uint16_t q = ip_Subsequences[ssIdx].q;
               
               // odd and even
               if (parity == SP_MIXED)
               {
                  tempQs[SP_MIXED][len[SP_MIXED]++] = q;
                     
                  if (ip_Subsequences[ssIdx].q%2 == 1)
                     tempQs[SP_ODD][len[SP_ODD]++] = q;

                  if (ip_Subsequences[ssIdx].q%2 == 0)
                     tempQs[SP_EVEN][len[SP_EVEN]++] = q;
               }
               
               if (parity == SP_ODD)
                  tempQs[SP_ODD][len[SP_ODD]++] = q;

               if (parity == SP_EVEN)
                  tempQs[SP_EVEN][len[SP_EVEN]++] = q;
            }
         }

         CopyQsAndMakeLadder(seq, SP_EVEN,  r, h, tempQs[SP_EVEN],  len[SP_EVEN]);
         CopyQsAndMakeLadder(seq, SP_ODD,   r, h, tempQs[SP_ODD],   len[SP_ODD]);
         CopyQsAndMakeLadder(seq, SP_MIXED, r, h, tempQs[SP_MIXED], len[SP_MIXED]);
      }
   }

   seq->congruentQCount = ii_TempQIndex;
   seq->congruentLadderCount = ii_TempLadderIndex;
   
   seq->congruentQs = (uint16_t *) xmalloc(seq->congruentQCount * sizeof(uint16_t));
   seq->congruentLadders = (uint16_t *) xmalloc(seq->congruentLadderCount * sizeof(uint16_t));
   
   mempcpy(seq->congruentQs, ip_TempQs, seq->congruentQCount * sizeof(uint16_t));
   mempcpy(seq->congruentLadders, ip_TempLadders, seq->congruentLadderCount * sizeof(uint16_t));
}

// Return true iff subsequence h of k*b^n+c has any terms with n%a==b.
bool  CisOneWithOneSequenceHelper::CongruentTerms(uint32_t ssIdx, uint32_t a, uint32_t b)
{
   uint32_t g, m;

   g = gcd32(a, ii_BestQ);
   
   if (b % g == ip_Subsequences[ssIdx].q % g)
   {
      for (m=ii_MinM; m<=ii_MaxM; m++)
      {
         if (ip_Subsequences[ssIdx].mTerms[m-ii_MinM])
            if ((m*ii_BestQ + ip_Subsequences[ssIdx].q) % a == b)
               return true;
      }
   }

   return false;
}

bool  CisOneWithOneSequenceHelper::CopyQsAndMakeLadder(seq_t *seq, sp_t parity, uint32_t r, uint32_t h, uint16_t *qList, uint32_t qListLen)
{
   if (qListLen == 0)
      return false;

   r = ip_PrlIndices[r];

   uint32_t cssIndex = CSS_INDEX(parity, r, h);
   
   // These will contain the index to the qListLen/qList and ladderLen/ladder for this parity, r, and h
   seq->congruentQIndices[cssIndex] = ii_TempQIndex;
   seq->congruentLadderIndices[cssIndex] = ii_TempLadderIndex;


   ip_TempQs[ii_TempQIndex] = qListLen;
   ii_TempQIndex++;

   if (qListLen > ii_MaxQs)
      ii_MaxQs = qListLen;
   
   memcpy(&ip_TempQs[ii_TempQIndex], qList, qListLen * sizeof(uint16_t));
   
   ii_TempQIndex += qListLen;
   
   MakeLadder(seq, qList, qListLen);

   return true;
}
         
void   CisOneWithOneSequenceHelper::MakeLadder(seq_t *seq, uint16_t *qList, uint32_t qListLen)
{
   uint32_t  i, j, k, a;
   uint8_t   tempQs[LIMIT_BASE+1];

   ii_LadderCount++;
   
   for (i=0; i<ii_BestQ; i++)
      tempQs[i] = 0;

   tempQs[ii_BestQ] = 1;
   
   for (i=0, a=1; i<qListLen; i++, a++)
      tempQs[qList[i]] = 1;

   for (i = 0; i < 3; i++)
   {
      if (tempQs[i] == 1)
        a--;
      tempQs[i] = 2;
   }

   while (a > 0)
   {
      for (i=3, j=2; i<=ii_BestQ; i++)
      {
         if (tempQs[i] == 2)
            j = i;
         else
            if (tempQs[i] == 1)
               break;
      }
      
      assert(i <= ii_BestQ);

      if (tempQs[i-j] == 2)
      {
         /* We can use an existing rung */
         tempQs[i] = 2;
         a--; 
      }
      else
      {
         /* Need to create a new rung */
         k = MIN(i-j,(i+1)/2); 
         assert(tempQs[k]==0);
         tempQs[k] = 1;
         a++;
         
         /* Need to re-check rungs above the new one */
         for (k++; k<=j; k++) 
            if (tempQs[k] == 2)
            {
               tempQs[k] = 1;
               a++;
            }
      }
   }

   a = 1;
   for (i=3; i <= ii_BestQ; i++)
      if (tempQs[i] == 2)
         a++;

   ip_TempLadders[ii_TempLadderIndex] = a;
   ii_TempLadderIndex++;

   if (a > ii_MaxLadderRungs)
      ii_MaxLadderRungs = a;

   j = 2;
   for (i=3; i<=ii_BestQ; i++)
      if (tempQs[i] == 2)
      {
         assert(tempQs[i-j]==2);
         ip_TempLadders[ii_TempLadderIndex] = i - j;
         ii_TempLadderIndex++;
         j = i;
      }
}
