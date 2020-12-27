/* GenericWorker.h -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GenericWorker_H
#define _GenericWorker_H

#include "SierpinskiRieselApp.h"
#include "AbstractSubsequenceHelper.h"
#include "AbstractWorker.h"
#include "../core/HashTable.h"

using namespace std;

class GenericWorker : public AbstractWorker
{
public:
   GenericWorker(uint32_t myId, App *theApp, AbstractSubsequenceHelper *appHelper);

   ~GenericWorker(void) {};

   void              SetSequences(uint64_t largestPrimeTested, uint32_t bestQ, seq_t *sequences, uint32_t sequenceCount, subseq_t *subsequences, uint32_t subsequenceCount);
   
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);
   
protected:

private:
   void              InitializeDiscreteLog(void);
   
   uint64_t          ProcessSmallPrimes(void);
   void              DiscreteLogSmallPrimes(uint64_t *primeList, uint64_t *bases);
   void              DiscreteLogLargePrimes(uint64_t *primeList);
   
   void              BabyStepsSmallPrimes(uint64_t *primeList, uint64_t *b, uint64_t *bj0, uint32_t *orderOfB);
   void              BabyStepsBigPrimes(uint64_t *primeList, double *invp, uint64_t *b, uint64_t *bj0, uint32_t *orderOfB);
   void              BuildTables(uint64_t *baseToUse, uint64_t *primeList, double *invp, uint64_t *bm64);
      
   uint32_t          ii_SieveLow;
   uint32_t          ii_SieveRange;

   uint64_t          il_SmallPrimeSieveLimit;
   uint64_t          il_GenericSeveLimit;
   
   uint32_t          ii_BabySteps;
   uint32_t          ii_GiantSteps;

   HashTable        *ip_HashTable[4];
   
   uint32_t         *ssHash;        // there is one per subsequence
   
   uint64_t        **ck64;          // there is one set of 4 per sequence
   uint64_t        **bdck64;        // there is one set of 4 per subsequence
   uint64_t        **bd64;          // there is one set of 4 per Q
};

#endif

