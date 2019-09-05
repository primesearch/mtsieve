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
#include "../opencl/Kernel.h"
#include "../opencl/KernelArgument.h"

using namespace std;

class CullenWoodallGpuWorker : public Worker
{
public:
   CullenWoodallGpuWorker(uint32_t myId, App *theApp);

   ~CullenWoodallGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   void              CreateKernels(bool knownWorkSize);
   void              VerifyFactor(uint64_t p, uint32_t n, int32_t c);
   
   uint64_t         *il_MagicNumber;
   uint64_t         *il_MagicShift;
   int64_t          *il_FactorList;

   uint32_t          ii_Base;
   uint32_t         *ii_Terms;
   uint32_t          ii_Groups;
   uint32_t          ii_MinGroupSize;
   uint32_t          ii_MaxGroupSize;
   uint32_t          ii_MaxTermCount;
   uint64_t          il_NextTermsBuild;
   int32_t           ii_MaxFactorCount;
   int32_t           ii_FactorCount;

   CullenWoodallApp *ip_CullenWoodallApp;
   
   Kernel           *ip_MagicKernel;
   Kernel           *ip_GCWKernel;

   KernelArgument   *ip_KAPrime;
   KernelArgument   *ip_MKAMagicNumber;
   KernelArgument   *ip_MKAMagicShift;

   KernelArgument   *ip_GCWKAMagicNumber;
   KernelArgument   *ip_GCWKAMagicShift;
   KernelArgument   *ip_KATerms;
   KernelArgument   *ip_KAFactorCount;
   KernelArgument   *ip_KAFactorList;
};

#endif

