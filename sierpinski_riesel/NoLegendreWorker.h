/* NoLegendreWorker.h -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _NoLegendreWorker_H
#define _NoLegendreWorker_H

#include "../x86_asm/avx-asm-x86.h"
#include "SierpinskiRieselApp.h"
#include "AbstractWorker.h"
#include "../core/HashTable.h"

using namespace std;

class NoLegendreWorker : public AbstractWorker
{
public:
   NoLegendreWorker(uint32_t myId, App *theApp, AbstractSubsequenceHelper *appHelper);

   ~NoLegendreWorker(void) {};

   void              SetSequences(uint64_t largestPrimeTested, uint32_t bestQ, seq_t *sequences, uint32_t sequenceCount, subseq_t *subsequences, uint32_t subsequenceCount);
   
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);
   
protected:

private:
   void              InitializeDiscreteLog(void);
   
   void              DiscreteLog(uint64_t *primeList);
   
   void              BabyStepsSSE(uint64_t *primeList, double *invp, uint64_t *b, uint64_t *bj0, uint32_t *orderOfB);
   void              BuildTablesSSE(uint64_t *baseToUse, uint64_t *primeList, uint64_t *bm64);

   void              BabyStepsAVX(double *dps, double *reciprocals, double *b, double *bj0, uint32_t *orderOfB);
   void              BuildTablesAVX(double *baseToUse, double *bm64d, uint64_t *primeList, double *dps, double *reciprocals);

   uint32_t          ii_SieveLow;
   
   HashTable        *ip_HashTable[4];
   
   uint32_t        **f32;     // there is one set of 4 per sequence   
   uint32_t        **g32;     // there is one set of 4 per sequence
   
   uint64_t        **c64;     // there is one set of 4 per subsequence
   uint64_t        **d64;     // there is one set of 4 per subsequence
};

#endif

