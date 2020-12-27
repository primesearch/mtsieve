/* CarolKyneaWorker.h -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CarolKynea_H
#define _CarolKynea_H

#include "CarolKyneaApp.h"
#include "../core/Worker.h"
#include "../core/HashTable.h"

// There are 4 sequences.  2 for the Carol form, 2 for the Kynea form
// All of them will be sieved concurrently.
#define ROOT_COUNT	4

using namespace std;

typedef struct {
   uint64_t root;
   int32_t  c;
} seq_t;

class CarolKyneaWorker : public Worker
{
public:
   CarolKyneaWorker(uint32_t myId, App *theApp);

   ~CarolKyneaWorker(void) {};

   void              CleanUp(void);
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   bool              VerifyFactor(bool badFactorIsFatal, uint64_t p, uint32_t n, int32_t c);
   
protected:

private:
   CarolKyneaApp    *ip_CarolKyneaApp;
   HashTable        *ip_HashTable;
   
   void              DiscreteLog(uint64_t p);
   uint32_t          BabySteps(uint64_t b, uint64_t bj0, uint64_t p);
   uint64_t          FindRoot(uint64_t p);
   void              CheckFactor(uint64_t p, uint32_t n, int32_t c);

   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint32_t          ii_BabySteps;
   uint32_t          ii_GiantSteps;

   uint32_t          ii_SieveLow;
   uint32_t          ii_SieveRange;

   seq_t             io_Sequence[ROOT_COUNT];
      
   uint64_t          il_A[ROOT_COUNT];
};    

#endif

