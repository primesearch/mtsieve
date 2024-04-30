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
#include "../core/MpArithVector.h"

using namespace std;

#define MAX_POWERS   50

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
   void           BuildTerms(void);

   uint32_t       ii_MinBase;
   uint32_t       ii_MaxBase;

   std::vector<bool>  iv_Bases;
   MpResVec      *ip_Remainders;
   uint64_t       il_NextTermsBuild;

   term_t        *ip_Terms;
};

#endif

