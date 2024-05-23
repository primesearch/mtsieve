/* LifchitzWorker.h -- (C) Mark Rodenkirch, April 2024

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _LifchitzWorker_H
#define _LifchitzWorker_H

#include "LifchitzApp.h"
#include "../core/Worker.h"
#include "../core/MpArith.h"

using namespace std;

#define MAX_POWERS   50

typedef struct {
   uint32_t       x;
   uint32_t       next;
} hash_t;

class LifchitzWorker : public Worker
{
public:
   LifchitzWorker(uint32_t myId, App *theApp);

   ~LifchitzWorker() {};

   void           TestMegaPrimeChunk(void);
   void           TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void           CleanUp(void);
   
protected:
   void           NotifyPrimeListAllocated(uint32_t primesInList) {}
   
private:
   LifchitzApp   *ip_LifchitzApp;

   uint32_t       ii_MinX;
   uint32_t       ii_MaxX;
   uint32_t       ii_MinY;
   uint32_t       ii_MaxY;
   
   uint32_t       ii_HashSize;
   uint32_t       ii_Elements;
   uint32_t       ii_HashM1;
   
   hash_t        *ip_Hashes;
   MpRes         *ip_Residues;
};

#endif

