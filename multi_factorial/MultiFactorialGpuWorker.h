/* MultiFactorialGpuWorker.h -- (C) Mark Rodenkirch, September 2016

   This is the header for the class sets up the call to the fsieve GPU function
   and parses the output from the GPU to determine if we have a factor.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _MultiFactorialGpuWorker_H
#define _MultiFactorialGpuWorker_H

#include "MultiFactorialApp.h"
#include "../core/Worker.h"
#include "../opencl/Kernel.h"
#include "../opencl/KernelArgument.h"

using namespace std;

class MultiFactorialGpuWorker : public Worker
{
public:
   MultiFactorialGpuWorker(uint32_t myId, App *theApp);

   ~MultiFactorialGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   void              CheckForPrimeOrFactor(int64_t p, uint32_t n, int32_t c);
   void              ValidateFactor(int64_t p, uint32_t n, int32_t c);
   
   MultiFactorialApp *ip_MultiFactorialApp;
   
   uint32_t          ii_StepN;
   uint32_t          ii_MultiFactorial;
   
   int64_t          *il_RemainderList;
   uint64_t         *il_MagicNumber;
   uint64_t         *il_MagicShift;
   uint32_t          ii_NRange[10];
   uint64_t         *il_MinusFactorList;
   uint64_t         *il_PlusFactorList;

   Kernel           *ip_MagicKernel;
   Kernel           *ip_FactorialKernel;

   KernelArgument   *ip_KAPrime;
   KernelArgument   *ip_MKAMagicNumber;
   KernelArgument   *ip_MKAMagicShift;

   KernelArgument   *ip_FKAMagicNumber;
   KernelArgument   *ip_FKAMagicShift;
   KernelArgument   *ip_KANRange;
   KernelArgument   *ip_KARemainder;
   KernelArgument   *ip_KAMinusFactor;
   KernelArgument   *ip_KAPlusFactor;
};

#endif

