/* XYYXSparseWorker.h -- (C) Mark Rodenkirch, June 2014

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _XYYXSparseWorker_H
#define _XYYXSparseWorker_H

#include "XYYXApp.h"
#include "../core/Worker.h"

using namespace std;

#define MAX_POWERS   50

class XYYXSparseWorker : public Worker
{
public:
   XYYXSparseWorker(uint32_t myId, App *theApp);

   ~XYYXSparseWorker() {};

   void           TestMegaPrimeChunk(void);
   void           TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void           CleanUp(void);
   
protected:
   void           NotifyPrimeListAllocated(uint32_t primesInList) {}
   
private:  
   XYYXApp       *ip_XYYXApp;

   term_t        *ip_Terms;
   uint64_t       il_NextTermsBuild;

   bool           ib_IsPlus;
   bool           ib_IsMinus;
};

#endif

