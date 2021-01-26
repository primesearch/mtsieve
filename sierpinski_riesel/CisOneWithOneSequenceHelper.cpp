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
#include "GenericGpuWorker.h"
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
      theWorker = new GenericGpuWorker(id, ip_App, this);
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
bool   CisOneWithOneSequenceHelper::BuildLegendreTableForSequence(seq_t *seq, legendre_t *legendre)
{
   int64_t  r;
   uint64_t mod, ssqfb;
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
         
         legendre->oneParityMap = (uint8_t *) xmalloc(L_BYTES(mod) * sizeof(uint8_t));
         BuildLegendreMap(mod, r, "single", legendre->oneParityMap);
         break;

      // odd and even n, test for (-ck/p)==1 and (-bck/p)==1
      default:
         mod = mod*2*ssqfb;
      
         if ((2*mod+1) >= INT32_MAX)
         {
            ip_App->WriteToConsole(COT_OTHER, "Cannot use Legendre tables because square-free part of k is too large");
            return false;
         }
         
         if ((r*ssqfb) >= INT32_MAX)
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
               
         legendre->dualParityMapP1 = (uint8_t *) xmalloc(L_BYTES(mod) * sizeof(uint8_t));
         legendre->dualParityMapM1 = (uint8_t *) xmalloc(L_BYTES(mod) * sizeof(uint8_t));

         BuildLegendreMap(mod, r, "even", legendre->dualParityMapP1);
         BuildLegendreMap(mod, r*ssqfb, "odd", legendre->dualParityMapM1);
         break;
   }
   
   legendre->mod = mod;
   return true;
}

void    CisOneWithOneSequenceHelper::BuildLegendreMap(uint32_t size, int32_t r, const char *which, uint8_t *legendreMap)
{
   uint32_t i;
   uint32_t reset = 0;

   for (i=0; i<size; i++)
   {
      // Provide progress when building a large table
      if (reset == 1000000)
      {
         ip_App->WriteToConsole(COT_OTHER, "Building %s parity Legendre symbol table: %.1f%%", which, 100.0*i/size);
         reset = 0;
      }

      if (jacobi(r, 2*i+1) == 1)
         legendreMap[L_BYTE(i)] |= L_BIT(i);

      reset++;
   }
}
