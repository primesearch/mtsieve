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
#include "../core/MpArith.h"

// There are 4 sequences.  2 for the Carol form, 2 for the Kynea form
// All of them will be sieved concurrently.
#define ROOT_COUNT	4

using namespace std;

class CarolKyneaWorker : public Worker
{
public:
   CarolKyneaWorker(uint32_t myId, App *theApp);

   ~CarolKyneaWorker(void) {};

   void              CleanUp(void);
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
      
protected:

private:
   CarolKyneaApp    *ip_CarolKyneaApp;
   HashTable        *ip_HashTable;
   
   MpRes             FindRoot(uint64_t p, MpArith mp);
   void              DiscreteLog(uint64_t p, MpArith mp, uint64_t inv_pb);
   uint32_t          BabySteps(uint64_t p, MpArith mp, MpRes resBase);

   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint32_t          ii_BabySteps;
   uint32_t          ii_GiantSteps;

   uint32_t          ii_SieveLow;
   uint32_t          ii_SieveRange;

   MpRes             ir_Root1, ir_Root2;
};    

#endif

