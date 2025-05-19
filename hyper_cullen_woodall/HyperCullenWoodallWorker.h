/* HyperCullenWoodallWorker.h -- (C) Mark Rodenkirch

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _HyperCullenWoodallWorker_H
#define _HyperCullenWoodallWorker_H

#include "HyperCullenWoodallApp.h"
#include "../core/MpArithVector.h"
#include "../core/Worker.h"

using namespace std;

class HyperCullenWoodallWorker : public Worker
{
public:
   HyperCullenWoodallWorker(uint32_t myId, App *theApp);

   ~HyperCullenWoodallWorker() {};

   void           TestMegaPrimeChunk(void);
   void           TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void           CleanUp(void);
   
protected:
   void           NotifyPrimeListAllocated(uint32_t primesInList) {}
   
private:         
   void           FreeTerms(void);
   
   void           ComputeResidues(MpArithVec mp, base_t *terms, uint32_t minPower);
   
   HyperCullenWoodallApp       *ip_HyperCullenWoodallApp;

   uint32_t       ii_MinB;
   uint32_t       ii_MaxB;
   uint32_t       ii_MinN;
   uint32_t       ii_MaxN;
   uint32_t       ii_BCount;
   uint32_t       ii_NCount;

   bool           ib_IsPlus;
   bool           ib_IsMinus;

   uint64_t       il_NextTermsBuild;
   base_t        *ip_bTerms;
   base_t        *ip_nTerms;
};

#endif

