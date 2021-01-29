/* CisOneSequenceHelper.cpp -- (C) Mark Rodenkirch, January 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <assert.h>
#include "CisOneSequenceHelper.h"
#include "CisOneWithOneSequenceWorker.h"
#include "../core/inline.h"

#define NBIT(n)         ((n) - ii_MinN)
#define MBIT(m)         ((m) - ii_MinM)

CisOneSequenceHelper::CisOneSequenceHelper(App *theApp, uint64_t largestPrimeTested) : AbstractSequenceHelper(theApp, largestPrimeTested)
{
   // The parent's constuctor does much of the work

   uint32_t shifts = POWER_RESIDUE_LCM / 2;
   uint32_t i, r;
   int32_t  shift, divide;
   
   ii_DivisorShifts = (int16_t *) xmalloc(shifts * sizeof(int16_t));

   for (i=0; i<shifts; i++)
   {
      r = gcd32(2*i, POWER_RESIDUE_LCM);
      
      divide = r;
      for (shift = 0; divide % 2 == 0; shift++)
         divide /= 2;
         
      if (divide > 1)
      {
         ii_DivisorShifts[i] = r;
         continue;
      }

      if (shift > 1)
         ii_DivisorShifts[i] = -shift;
      else
         ii_DivisorShifts[i] = 0;
   }
   
   ii_MaxLadderRungs = 0;
   ii_MaxQs = 0;
}

void   CisOneSequenceHelper::CleanUp(void)
{
   uint32_t i, j;
   for (i=0; i<SP_SIZE; i++)
   {
      for (j=0; j<POWER_RESIDUE_LCM; j++)
      {
         if (ii_CssQs[i][j] != NULL)
         {
            xfree(ii_CssQs[i][j]);
            xfree(ii_CSSLadders[i][j]);
         }
      }
      
      xfree(ii_CssQs[i]);
      xfree(ii_CSSLadders[i]);
   }
}

void        CisOneSequenceHelper::LastChanceLogicBeforeSieving(void)
{
   uint32_t ssIdx, babySteps, giantSteps;

   ii_MaxBabySteps = 0;
   
   for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
   {
      ChooseSteps(ii_BestQ, ssIdx+1, babySteps, giantSteps);
      
      if (ii_SieveHigh >= ii_SieveLow + babySteps*giantSteps)
         FatalError("LastChanceLogicBeforeSieving miscomputed the steps");
         
      ip_Subsequences[ssIdx].babySteps = babySteps;
      ip_Subsequences[ssIdx].giantSteps = giantSteps;

      if (ii_MaxBabySteps < babySteps)
         ii_MaxBabySteps = babySteps;
   }
   
   MakeSubseqCongruenceTables();

   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) ip_App;
   
   ib_UseLegendreTables = srApp->UseLegendreTables();

   if (ib_UseLegendreTables)
      ib_UseLegendreTables = BuildLegendreTables(srApp->GetLegendreFileName());
}

uint32_t    CisOneSequenceHelper::FindBestQ(uint32_t &expectedSubsequences)
{
   uint32_t   i = 0, j, n;
   uint32_t   bit;
   uint32_t   S[NDIVISORS];
   uint32_t   W[NDIVISORS];
   bool       R[LIMIT_BASE];
   seq_t     *seq;

   seq = ip_FirstSequence;
   do
   {
      for (j=0; j<NDIVISORS; j++)
         S[j] = 0;

      for (j=0; j<LIMIT_BASE; j++)
         R[j] = false;
      
      bit = NBIT(ii_MinN);

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seq->nTerms[bit])
            R[n%LIMIT_BASE] = true;
         
         bit++;
      }
 
      i = 0;
      for (j=0; j<NDIVISORS; j++)
      {
         if (NDIVISORS % (j+1) == 0)
         {
            S[j] += CountResidueClasses((j+1)*BASE_MULTIPLE, LIMIT_BASE, R);
            
            W[j] = RateQ((j+1)*BASE_MULTIPLE, S[j]);
            
            if (W[j] < W[i])
               i = j;
         }
      }

      seq = (seq_t *) seq->next;
   } while (seq != NULL);

   expectedSubsequences = S[i];
   
   return (i+1)*BASE_MULTIPLE;
}

// Build tables sc_lists[i][j] of pointers to lists of subsequences
// whose terms k*b^m+c satisfy m = j (mod r)
void  CisOneSequenceHelper::MakeSubseqCongruenceTables(void)
{
   uint32_t   ssIdx, h, i, j, k, len[SP_SIZE];
   uint16_t   tempQs[SP_SIZE][POWER_RESIDUE_LCM];
   uint16_t  *qList, qLen;
   sp_t       parity = ip_FirstSequence->parity;
   uint64_t   preBuildCpuBytes;
   uint64_t   postBuildCpuBytes;

   preBuildCpuBytes = GetCpuMemoryUsage();
   
   for (uint32_t i=0; i<3; i++)
   {
      ii_CssQs[i] = (uint16_t ***) xmalloc(POWER_RESIDUE_LCM * sizeof(uint16_t **));
      ii_CSSLadders[i] = (uint16_t ***) xmalloc(POWER_RESIDUE_LCM * sizeof(uint16_t **));
   }
   
   for (j=1; j<=POWER_RESIDUE_LCM; j++)
   {
      if (POWER_RESIDUE_LCM % j != 0)
         continue;
            
      if (parity == SP_MIXED)
      {
         ii_CssQs[SP_MIXED][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
         ii_CSSLadders[SP_MIXED][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
         
         ii_CssQs[SP_EVEN][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
         ii_CSSLadders[SP_EVEN][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
         
         ii_CssQs[SP_ODD][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
         ii_CSSLadders[SP_ODD][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
      }

      if (parity == SP_EVEN)
      {
         ii_CssQs[SP_EVEN][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
         ii_CSSLadders[SP_EVEN][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
      }
      
      if (parity == SP_ODD)
      {
         ii_CssQs[SP_ODD][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
         ii_CSSLadders[SP_ODD][j] = (uint16_t **) xmalloc(j*sizeof(uint16_t *));
      }
      
      for (k=0; k<j; k++)
      {
         len[0] = len[1] = len[2] = 0;
         
         for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
         {
            if (CongruentTerms(ssIdx, j, k))
            {
               int16_t q = ip_Subsequences[ssIdx].q;
               
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
         
         for (i=0; i<SP_SIZE; i++)
         {
            if (len[i] == 0)
               continue;

            qLen = len[i];

            if (qLen > ii_MaxQs)
               ii_MaxQs = qLen;
            
            qList = (uint16_t *) xmalloc((qLen+5)*sizeof(uint16_t));
            qList[0] = qLen;
            
            for (h=0; h<len[i]; h++)
               qList[h+1] = tempQs[i][h];

            ii_CssQs[i][j][k] = qList;
            ii_CSSLadders[i][j][k] = MakeLadder(tempQs[i], len[i]);
         }
      }
   }
   
   postBuildCpuBytes = GetCpuMemoryUsage();
   
   ip_App->WriteToConsole(COT_OTHER, "%llu bytes used for congruence tables.  Max Qs = %u.  Max ladder size = %u", postBuildCpuBytes - preBuildCpuBytes, ii_MaxQs, ii_MaxLadderRungs);
}

// Return true iff subsequence h of k*b^n+c has any terms with n%a==b.
bool  CisOneSequenceHelper::CongruentTerms(uint32_t ssIdx, uint32_t a, uint32_t b)
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

uint16_t *CisOneSequenceHelper::MakeLadder(uint16_t *qs, uint32_t qLen)
{
   uint32_t  i, j, k, a;
   uint16_t *ladder;
   uint8_t   qList[LIMIT_BASE+1];

   for (i=0; i<ii_BestQ; i++)
      qList[i] = 0;

   qList[ii_BestQ] = 1;
   
   for (i=0, a=1; i<qLen; i++, a++)
      qList[qs[i]] = 1;

   for (i = 0; i < 3; i++)
   {
      if (qList[i] == 1)
        a--;
      qList[i] = 2;
   }

   while (a > 0)
   {
      for (i=3, j=2; i<=ii_BestQ; i++)
      {
         if (qList[i] == 2)
            j = i;
         else
            if (qList[i] == 1)
               break;
      }
      
      assert(i <= ii_BestQ);

      if (qList[i-j] == 2)
      {
         /* We can use an existing rung */
         qList[i] = 2;
         a--; 
      }
      else
      {
         /* Need to create a new rung */
         k = MIN(i-j,(i+1)/2); 
         assert(qList[k]==0);
         qList[k] = 1;
         a++;
         
         /* Need to re-check rungs above the new one */
         for (k++; k<=j; k++) 
            if (qList[k] == 2)
            {
               qList[k] = 1;
               a++;
            }
      }
   }

   a = 1;
   for (i=3; i <= ii_BestQ; i++)
      if (qList[i] == 2)
         a++;

   ladder = (uint16_t *) xmalloc((a+1)*sizeof(uint16_t));
   ladder[0] = a;
   
   if (a > ii_MaxLadderRungs)
      ii_MaxLadderRungs = a;

   j = 2;
   k = 1;
   for (i=3; i<=ii_BestQ; i++)
      if (qList[i] == 2)
      {
         assert(qList[i-j]==2);
         ladder[k] = i-j;
         j = i;
         k++;
      }
   
   assert(k == a);
         
   return ladder;
}

bool   CisOneSequenceHelper::BuildLegendreTables(string legendreFileName)
{
   uint64_t     preBuildCpuBytes;
   uint64_t     postBuildCpuBytes;
   legendre_t  *legendre;
   bool         builtLegendreTables;
   seq_t       *seq;
      
   preBuildCpuBytes = GetCpuMemoryUsage();
   
   seq = ip_FirstSequence;
   do
   {
      legendre = (legendre_t *) xmalloc(sizeof(legendre_t));
      
      seq->legendrePtr = legendre;

      builtLegendreTables = BuildLegendreTableForSequence(seq, legendre);
      
      if (!builtLegendreTables)
         break;

      seq = (seq_t *) seq->next;
   } while (seq != NULL);

   postBuildCpuBytes = GetCpuMemoryUsage();
   
   if (builtLegendreTables)
      ip_App->WriteToConsole(COT_OTHER, "%llu bytes used for Legendre tables", postBuildCpuBytes - preBuildCpuBytes);

   return builtLegendreTables;
}
