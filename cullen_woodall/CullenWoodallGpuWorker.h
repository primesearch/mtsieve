/* CullenWoodallGpuWorker.h -- (C) Mark Rodenkirch, May 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CullenWoodallGpuWorker_H
#define _CullenWoodallGpuWorker_H

#include "CullenWoodallApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

class CullenWoodallGpuWorker : public Worker
{
public:
   CullenWoodallGpuWorker(uint32_t myId, App *theApp);

   ~CullenWoodallGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   uint64_t          il_NextTermsBuild;   
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
   uint32_t         *ii_Terms;

   uint32_t          ii_Base;
   uint32_t          ii_Groups;
   uint32_t          ii_GroupSize;

   // Variables used by the kernel
   uint32_t         *ii_KernelTerms;
   uint32_t         *ii_FactorCount;
   int64_t          *il_FactorList;

   CullenWoodallApp *ip_CullenWoodallApp;
   
   GpuKernel        *ip_Kernel;
};

#endif

