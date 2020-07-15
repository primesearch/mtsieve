/* XYYXGpuWorker.h -- (C) Mark Rodenkirch, June 2014

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _XYYXGpuWorker_H
#define _XYYXGpuWorker_H

#include "XYYXApp.h"
#include "../core/Worker.h"
#include "../opencl/Kernel.h"
#include "../opencl/KernelArgument.h"

using namespace std;

class XYYXGpuWorker : public Worker
{
public:
   XYYXGpuWorker(uint32_t myId, App *theApp);

   ~XYYXGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   void              VerifyFactor(uint64_t p, uint32_t x, uint32_t y, int32_t c);
   
   uint64_t         *il_MagicNumber;
   uint64_t         *il_MagicShift;
   uint64_t         *il_MinusFactorList;
   uint64_t         *il_PlusFactorList;
   
   uint32_t          ii_Groups;
   uint32_t          ii_GroupSize;
   uint32_t         *ii_GroupedTerms;
   uint64_t          il_NextTermsBuild;

   XYYXApp          *ip_XYYXApp;

   Kernel           *ip_MagicKernel;
   Kernel           *ip_XYYXKernel;

   KernelArgument   *ip_KAPrime;
   KernelArgument   *ip_MKAMagicNumber;
   KernelArgument   *ip_MKAMagicShift;

   KernelArgument   *ip_XYYXKAMagicNumber;
   KernelArgument   *ip_XYYXKAMagicShift;

   KernelArgument   *ip_FKAMagicNumber;
   KernelArgument   *ip_FKAMagicShift;
   KernelArgument   *ip_KATerms;
   KernelArgument   *ip_KAMinusFactors;
   KernelArgument   *ip_KAPlusFactors;
};

#endif

