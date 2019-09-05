/* NoLegendreWorker.cpp -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/sse-asm-x86.h"
#include "NoLegendreWorker.h"

NoLegendreWorker::NoLegendreWorker(uint32_t myId, App *theApp, AbstractSubsequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  NoLegendreWorker::CleanUp(void)
{
   uint32_t idx;

   xfree(ssHash);
   
   for (idx=0; idx<4; idx++)
      delete ip_HashTable[idx];
   
   for (idx=0; idx<ii_SubsequenceCount; idx++)
   {
      xfree(c32[idx]);
      xfree(d64[idx]);
   }

   xfree(c32);
   xfree(d64);
   
   for (idx=0; idx<ii_SequenceCount; idx++)
   {
      xfree(f32[idx]);
      xfree(g64[idx]);
   }

   xfree(f32);
   xfree(g64);
}

void  NoLegendreWorker::TestMegaPrimeChunk(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t primeList[4];
   
   vector<uint64_t>::iterator it = iv_Primes.begin();
   
   while (it != iv_Primes.end())
   {
      primeList[0] = *it;
      it++;
      
      primeList[1] = *it;
      it++;
      
      primeList[2] = *it;
      it++;
      
      primeList[3] = *it;
      it++;

      DiscreteLog(primeList);
      
      SetLargestPrimeTested(primeList[3], 4);
      
      if (primeList[3] >= maxPrime)
         break;
   }
}

void  NoLegendreWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("NoLegendreWorker::TestMiniPrimeChunk not implemented");
}

void  NoLegendreWorker::SetSequences(uint64_t largestPrimeTested, uint32_t bestQ, seq_t *sequences, uint32_t sequenceCount, subseq_t *subsequences, uint32_t subsequenceCount)
{
   if (ip_WorkerSequences)
      CleanUp();
   
   ii_BestQ = bestQ;
   
   ip_Sequences = sequences;
   ii_SequenceCount = sequenceCount;

   ip_Sequences = subsequences;
   ii_SubsequenceCount = subsequenceCount;
   
   InitializeDiscreteLog();
}

void  NoLegendreWorker::InitializeDiscreteLog(void)
{
   ii_BabySteps = 0;
      
   for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
      ii_BabySteps = MAX(ii_BabySteps, ip_SubSequenes[ssIdx].babySteps);

   ii_SieveLow = ii_MinN / ii_BestQ;

   f32 = (uint32_t **) xmalloc(ii_SequenceCount*sizeof(uint32_t *));
   g64 = (uint64_t **) xmalloc(ii_SequenceCount*sizeof(uint64_t *));
   
   c32 = (uint32_t **) xmalloc(ii_SubsequenceCount*sizeof(uint32_t *));
   d64 = (uint64_t **) xmalloc(ii_SubsequenceCount*sizeof(uint64_t *));
   
   for (idx=0; idx<ii_SequenceCount; idx++) 
   {
      f32[idx] = (uint32_t *) xmalloc(4*sizeof(uint32_t));
      g64[idx] = (uint64_t *) xmalloc(4*sizeof(uint64_t));
   }

   for (idx=0; idx<ii_SubsequenceCount; idx++) 
   {
      c32[idx] = (uint32_t *) xmalloc(4*sizeof(uint32_t));
      d64[idx] = (uint64_t *) xmalloc(4*sizeof(uint64_t));
   }
   
   for (idx=0; idx<4; idx++)
      ip_HashTable[idx] = new HashTable(ii_BabySteps);
}

void  NoLegendreWorker::DiscreteLog(uint64_t *primeList)
{
   uint32_t i, j, ssIdx;
   uint32_t pIdx;
   uint32_t orderOfB[4];   
   uint64_t b[4], bm64[4], bj0[4];
   double   invp[4];
      
   b[0] = ii_Base;
   b[1] = ii_Base;
   b[2] = ii_Base;
   b[3] = ii_Base;
   
   invp[0] = 1.0 / (double) primeList[0];
   invp[1] = 1.0 / (double) primeList[1];
   invp[2] = 1.0 / (double) primeList[2];
   invp[3] = 1.0 / (double) primeList[3];
   
  uint32_t i, j, kk, cc;
  uint64_t b, inv_b;
  uint64_t bQ;

  /* inv_b <-- 1/base (mod p) */
#if (BASE == 2)
  inv_b = (p+1)/2;
#elif (BASE == 0)
  inv_b = (b_term == 2) ? (p+1)/2 : invmod32_64(b_term,-1,p);
#else
  inv_b = invmod32_64(BASE,-1,p);
#endif

  vec_mulmod64_initp(p);

  /* Swap b and inv_b for dual sieve */
#if DUAL
  if (dual_opt)
  {
    b = inv_b;
# if (BASE == 0)
    inv_b = b_term;
# else
    inv_b = BASE;
# endif
  }
  else
#endif
# if (BASE == 0)
    b = b_term;
# else
    b = BASE;
# endif

  bQ = climb_ladder(BD64,b,p); /* bQ <-- b^Q. */

  if ((cc = setup64(inv_b,p)) > 0)
  {
    m = steps[cc].m;
    M = steps[cc].M;

    /* Baby steps. */
    if ((i = baby_steps(inv_b,p)) > 0) /* Unlikely */
    {
      /* i is the order of b (mod p). This is all the information we need to
         determine every solution for this p, so no giant steps are needed.
      */
      for (kk = 0; kk < cc; kk++)
        for (j = lookup64(D64[kk]); j < m*M; j += i)
          eliminate_term(C32[kk],sieve_low+j,p);
    }
    else
    {
      /* First giant step. */
      if (search_hashtable(D64,cc) < cc)
        eliminate_hashtable_terms(0,cc,p);

      /* Remaining giant steps. */
      if (M > 1)
#ifdef HAVE_new_giant_steps
        giant_steps(D64,cc,M,eliminate_1,powmod64(bQ,m,p),p);
#else
        giant_steps(cc,powmod64(bQ,m,p),p);
#endif
    }
  }

  vec_mulmod64_finip();
}

/* BD64 is a temporary array of Q+1 elements. For SSE2 operations BD64[1]
   must be at least 16-aligned, but 64-aligned is better.
*/
#define BD64 (BJ64+((POWER_RESIDUE_LCM+1)|7))

/* This function builds the list C32[] of subsequences (k*b^d)*(b^Q)^m+c for
   which p may be a factor (-ckb^d is a quadratic/cubic/quartic/quintic
   residue with respect to p) and initialises the table D64[] with the
   values -c/(k*b^d) (mod p). Returns the number of subsequences in C32[].
*/
static uint32_t attribute ((noinline)) setup64(uint64_t inv_b, uint64_t p)
{
  uint64_t neg_ck, bm, p_s;
#if NO_LOOKUP_OPT
  int64_t b_sym = 0, kc_sym;
#endif
  uint32_t div_ind;
  const uint32_t *list;
  uint32_t r, f, g, h, i, j;
  int16_t s;

  bm = p/2;

#if NO_LOOKUP_OPT
  if (no_lookup_opt)
  {
# if (BASE==0)
    b_sym = legendre32(b_term,p);
# else
    b_sym = legendre32(BASE,p);
# endif
  }
  else
#endif /* NO_LOOKUP_OPT */
  {
#ifdef HAVE_lookup_ind
    gen_lookup_ind(F32,lookup_ind_data,seq_count,bm);
#else
    for (i = 0; i < seq_count; i++)
      F32[i] = bm%SEQ[i].mod;
#endif /* HAVE_lookup_ind */
  }

#if SKIP_CUBIC_OPT
  if (skip_cubic_opt)
    s = 0;
  else
#endif
  s = div_shift[bm%(POWER_RESIDUE_LCM/2)];

  if (s == 0)
  {
    /* p = 1 (mod 2) is all we know, check for quadratic residues only.
     */
#if NO_LOOKUP_OPT
    if (no_lookup_opt)
      for (i = j = 0; i < seq_count; i++)
      {
        kc_sym = legendre64(SEQ[i].kc_core,p);
        switch (SEQ[i].parity)
        {
          default:
          case 0: /* even n */
            kc_sym = (kc_sym == 1);
            break;
          case 1: /* odd n */
            kc_sym = (kc_sym == b_sym);
            break;
#if ALLOW_MIXED_PARITY_SEQUENCES
          case 2: /* even and odd n */
            kc_sym = (kc_sym == 1 || kc_sym == b_sym);
            break;
#endif
        }
        if (kc_sym)  /* (-ckb^n/p) == +1 */
        {
          /* For each subsequence (k*b^d)*(b^Q)^(n/Q)+c, compute
             -c/(k*b^d) (mod p) and add the subsequence to the bsgs list.
          */
          neg_ck = (SEQ[i].c < 0) ? SEQ[i].k : p - SEQ[i].k;
          PRE2_MULMOD64_INIT(neg_ck);
          for (h = SEQ[i].first; h <= SEQ[i].last; h++)
          {
            C32[j] = h;
            D64[j] = PRE2_MULMOD64(BD64[subseq_d[h]],neg_ck,p); /* -c/(k*b^d) */
            j++;
          }
          PRE2_MULMOD64_FINI();
        }
      }
    else
#endif /* NO_LOOKUP_OPT */
      for (i = j = 0; i < seq_count; i++)
      {
        if (test_bit(SEQ[i].map,F32[i]))  /* (-ckb^n/p) == +1 */
        {
          /* For each subsequence (k*b^d)*(b^Q)^(n/Q)+c, compute
             -c/(k*b^d) (mod p) and add the subsequence to the bsgs list.
          */
          neg_ck = (SEQ[i].c < 0) ? SEQ[i].k : p - SEQ[i].k;
          PRE2_MULMOD64_INIT(neg_ck);
          for (h = SEQ[i].first; h <= SEQ[i].last; h++)
          {
            C32[j] = h;
            D64[j] = PRE2_MULMOD64(BD64[subseq_d[h]],neg_ck,p); /* -c/(k*b^d) */
            j++;
          }
          PRE2_MULMOD64_FINI();
        }
      }

    return j;
  }
  else if (s > 0)
  {
    /* p = 1 (mod s), where s is not a power of 2. Check for r-th power
       residues for each prime power divisor r of s.
    */
    p_s = p/s;
  }
  else /* s < 0 */
  {
    /* p = 1 (mod 2^s), where s > 1. Check for r-th power residues for each
       divisor r of s. We handle this case seperately to avoid computing p/s
       using plain division.
    */
    p_s = p >> (-s);
#ifndef NDEBUG
    s = 1 << (-s);
#endif
  }

  /* For 0 <= r < s, BJ64[r] <- 1/(b^r)^((p-1)/s) */
  BJ64[0] = 1;
  BJ64[1] = powmod64(inv_b,p_s,p);
  PRE2_MULMOD64_INIT(BJ64[1]);
  for (r = 1; BJ64[r] != 1; r++)
    BJ64[r+1] = PRE2_MULMOD64(BJ64[r],BJ64[1],p);
  PRE2_MULMOD64_FINI();
  assert(s%r == 0);
  /* 1/(b^r)^((p-1)/s)=1 (mod p) therefore (1/(b^r)^((p-1)/s))^y=1 (mod p)
     for 0 <= y < s/r. (Could we do more with this?)
  */
  div_ind = divisor_index[r];

#if USE_SETUP_HASHTABLE
  if (r > SMALL_HASH_THRESHOLD)
  {
    clear_hashtable(SMALL_HASH_SIZE);
    for (j = 0; j < r; j++)
      insert64_small(j);
  }
#endif

#if NO_LOOKUP_OPT
  if (no_lookup_opt)
    for (i = g = 0; i < seq_count; i++)
    {
      kc_sym = legendre64(SEQ[i].kc_core,p);
      switch (SEQ[i].parity)
      {
        default:
        case 0: /* even n */
          kc_sym = (kc_sym == 1);
          break;
        case 1: /* odd n */
          kc_sym = (kc_sym == b_sym);
          break;
#if ALLOW_MIXED_PARITY_SEQUENCES
        case 2: /* even and odd n */
          kc_sym = (kc_sym == 1 || kc_sym == b_sym);
          break;
#endif
      }
      if (kc_sym)  /* (-ckb^n/p) == +1 */
      {
        G64[g] = (SEQ[i].c < 0) ? SEQ[i].k : p - SEQ[i].k; /* -k/c (mod p) */
        F32[g] = i;
        g++;
      }
    }
  else
#endif /* NO_LOOKUP_OPT */
  for (i = g = 0; i < seq_count; i++)
  {
    if (test_bit(SEQ[i].map,F32[i]))  /* (-ckb^n/p) == +1 */
    {
      G64[g] = (SEQ[i].c < 0) ? SEQ[i].k : p - SEQ[i].k; /* -k/c (mod p) */
      F32[g] = i;
      g++;
    }
  }

#ifdef HAVE_vec_powmod64
  if (g > 1)
  {
    int remaining = g;
    // Because the underlying asm assumes that it can grab sizeof(uint64_t)*g bytes
    // from the stack, this can crash, so we reduce the size for each call to avoid it.
    for (i=0; i<g; i+=32*4)
    {
      if (remaining < 32*4)
        vec_powmod64(&G64[i], remaining, p_s, p);
      else
        vec_powmod64(&G64[i], 32*4, p_s, p);
     
      remaining -= 32*4;
    }
  }
  else
#endif
    for (i = 0; i < g; i++)
      G64[i] = powmod64(G64[i],p_s,p);

  for (i = j = 0; i < g; i++)
  {
    /* BJ64[r] <-- (-k/c)^((p-1)/s) */
    f = F32[i];
    BJ64[r] = G64[i];

    /* Find h such that BJ64[h]=BJ64[r], i.e. (-ckb^h)^((p-1)/r)=1 (mod p),
       or h >= r if not found. */
#if USE_SETUP_HASHTABLE
    if (r > SMALL_HASH_THRESHOLD)
      h = lookup64_small(BJ64[r]);
    else
#endif
    { /* Linear search */
      for (h = 0; BJ64[h] != BJ64[r]; h++)
        ;
    }

    if (h < r && (list = SCL[f].sc_lists[div_ind][h]) != NULL)
    {
      /* -c/(k*b^n) is an r-power residue for at least one term k*b^n+c
         of this sequence.
      */
      neg_ck = (SEQ[f].c < 0) ? SEQ[f].k : p - SEQ[f].k;
      PRE2_MULMOD64_INIT(neg_ck);
      while ((h = *list++) < SUBSEQ_MAX)
      {
        /* -ckb^d is an r-th power residue for at least one term
           (k*b^d)*(b^Q)^(n/Q)+c of this subsequence.
        */
        C32[j] = h;
        D64[j] = PRE2_MULMOD64(BD64[subseq_d[h]],neg_ck,p); /* -c/(k*b^d) */
        j++;
      }
      PRE2_MULMOD64_FINI();
    }
  }

  return j;
}
