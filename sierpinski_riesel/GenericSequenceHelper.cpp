/* GenericSequenceHelper.cpp -- (C) Mark Rodenkirch, May 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "GenericSequenceHelper.h"
#include "GenericWorker.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "GenericGpuWorker.h"
#endif

#define NBIT(n)         ((n) - ii_MinN)
#define MBIT(m)         ((m) - ii_MinM)

GenericSequenceHelper::GenericSequenceHelper(App *theApp, uint64_t largestPrimeTested) : AbstractSequenceHelper(theApp, largestPrimeTested)
{
   // The parent's constuctor does the work
   theApp->WriteToConsole(COT_OTHER, "Sieving with generic logic for p >= %" PRIu64"", largestPrimeTested);
}

Worker   *GenericSequenceHelper::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   AbstractWorker *theWorker;

   // Note that GenericWorker inherits from Worker.  This will not
   // only create the worker, but also start it.
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   if (gpuWorker)
      theWorker = new GenericGpuWorker(id, ip_App, this);
   else
#endif
      theWorker = new GenericWorker(id, ip_App, this);
      
   theWorker->Prepare(largestPrimeTested, ii_BestQ);

   return theWorker;
}

uint32_t    GenericSequenceHelper::FindBestQ(uint32_t &expectedSubsequences)
{
   uint32_t      bit, i, j, k, n;
   uint32_t      userBestQ;
   uint32_t      bestQ, Q = 2880;  // DEFAULT_LIMIT_BASE from the old code
   std::vector<bool>  R;
   choice_bc_t  *S = 0;
   seq_t        *seqPtr;

   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) ip_App;
   userBestQ = srApp->GetUserBestQ();

   k = ComputePossibleQ(Q, R, S);
   
   R.resize(Q, false);
   S = (choice_bc_t *) xmalloc(k, sizeof(choice_bc_t), "S");

   for (i = 0; i < k; i++)
      S[i].subseqs = 0;

   seqPtr = ip_FirstSequence;
   do
   {
      for (j=0; j<Q; j++)
         R[j] = 0;
      
      bit = NBIT(ii_MinN);

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seqPtr->nTerms[bit])
            R[n%Q] = 1;
         
         bit++;
      }
      
      ForEachDivisor(seqPtr, Q, R, S);
      
      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
      
   j = 0;
   for (i = 0; i < k; i++)
   {      
      S[i].work = EstimateWork(S[i].div, S[i].subseqs, true);
            
      if (S[i].work < S[j].work)
         j = i;
   }
   
   for (i = 0; i < k; i++)
   {
      if (userBestQ == S[i].div)
         j = i;
   }

   bestQ = S[j].div;
   expectedSubsequences = S[j].subseqs;
   
   if (userBestQ > 0 && bestQ != userBestQ)
      srApp->WriteToConsole(COT_OTHER, "Q of %u is not possible.  Using Q = %u", userBestQ, bestQ);

   xfree(S);
 
   return bestQ;
}

uint32_t GenericSequenceHelper::ComputePossibleQ(uint32_t Q, std::vector<bool> R, choice_bc_t *S)
{
   uint32_t A[9], P[9], M[9], i, k, t;

   k = FindMultiplicities(Q, P, M);

   for (i = 0, t = 1; i < k; i++)
   {
      A[i] = 0;
      t *= (M[i]+1);
   }

   return t;
}

uint32_t GenericSequenceHelper::ForEachDivisor(seq_t *seq, uint32_t Q, std::vector<bool> R, choice_bc_t *S)
{
   uint32_t A[9], P[9], M[9], d, i, j, k, t, subseqs;

   k = FindMultiplicities(Q, P, M);

   for (i = 0, t = 1; i < k; i++)
   {
      A[i] = 0;
      t *= (M[i]+1);
   }

   uint32_t q = 0;
   double qWork = 999999999.0;
   
   for (i = 0; i < t; i++)
   {
      for (j = 0; j < k; A[j++] = 0)
         if (++A[j] <= M[j])
            break;
         
      for (j = 0, d = 1; j < k; j++)
         d *= pow32(P[j], A[j]);

      subseqs = CountResidueClasses(d, Q, R);

      double work = EstimateWork(d, subseqs, false);

      if (work < qWork)
      {
         q = d;
         qWork = work;
      }

      S[i].div = d;
      S[i].subseqs += subseqs;
   }

   seq->bestQ = q;
   //printf("best for %llu %u %.0lf\n", seq->k, q, qWork);
      
   return t;
}

// Express n as n = p_0^m_0 * p_1^m_1 * ... * p_(k-1)^m_(k-1) for prime p_i
// and positive m_j. Store p_i in P[i], m_j in M[j], and return k. Since n
// is a uint32_t, k cannot be greater than 9.
uint32_t    GenericSequenceHelper::FindMultiplicities(uint32_t n, uint32_t *P, uint32_t *M)
{
   uint32_t i, m, q;

   for (i = 0, q = 2; n > 1; q = (q + 1) | 1)
   {
      for (m = 0; n % q == 0; m++)
         n /= q;
      
      if (m > 0)
         P[i] = q, M[i] = m, i++;
   }

   return i;
}

#define BABY_WORK    1.0    // 1 mulmod, 1 insert
#define GIANT_WORK   1.0    // 1 mulmod, 1 lookup
#define EXP_WORK     0.7    // 1 mulmod
#define SUBSEQ_WORK  1.4    // 1 mulmod, 1 lookup (giant step 0)
                               
double    GenericSequenceHelper::EstimateWork(uint32_t Q, uint32_t s, bool print)
{
   uint32_t babySteps, giantSteps;
   double   work;
   
   ChooseSteps(Q, s, babySteps, giantSteps);

   work = babySteps*BABY_WORK + s*(giantSteps-1)*GIANT_WORK + Q*EXP_WORK + s*SUBSEQ_WORK;

   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) ip_App;
   
   if (print && srApp->ShowQEffort())
      srApp->WriteToConsole(COT_OTHER, "q = %4u with %7u subseq yields bs = %5u, gs = %5u, work = %8.0lf", Q, s, babySteps, giantSteps, work);
   
   return work;
}