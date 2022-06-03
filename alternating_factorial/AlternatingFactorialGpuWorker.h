/* AlternatingFactorialGpuWorker.h -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _AlternatingFactorialGpuWorker_H
#define _AlternatingFactorialGpuWorker_H

#include "AlternatingFactorialApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

class AlternatingFactorialGpuWorker : public Worker
{
public:
   AlternatingFactorialGpuWorker(uint32_t myId, App *theApp);

   ~AlternatingFactorialGpuWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:   
   AlternatingFactorialApp      *ip_AlternatingFactorialApp;

   uint32_t          ii_MaxGpuFactors;
   
   uint64_t         *il_FactorialResiduals;
   uint64_t         *il_AltFactorialResiduals;
   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;
   uint32_t         *ii_Parameters;

   GpuKernel        *ip_Kernel;
};

#endif

