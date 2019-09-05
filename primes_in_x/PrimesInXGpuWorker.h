/* PixSieve.h -- (C) Mark Rodenkirch, September 2016

   This is the header for the class sets up the call to the fsieve GPU function
   and parses the output from the GPU to determine if we have a factor.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _PrimesInXGpuWorker_H
#define _PrimesInXGpuWorker_H

#include "PrimesInXApp.h"
#include "../core/Worker.h"
#include "../opencl/Kernel.h"
#include "../opencl/KernelArgument.h"

using namespace std;

class PrimesInXGpuWorker : public Worker
{
public:
   PrimesInXGpuWorker(uint32_t myId, App *theApp);

   ~PrimesInXGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   void              VerifyFactor(uint64_t prime, uint32_t n);
   
   void              BuildMixedTerms(void);
   void              BuildSingleTerms(void);

   uint32_t          ii_MixedDigitListCount;
   uint32_t          ii_MixedDigitListUsed;
   uint32_t        **ii_MixedDigitList;

   uint32_t          ii_SingleDigitListCount;
   uint32_t          ii_SingleDigitListUsed;
   uint32_t        **ii_SingleDigitList;
   
   uint32_t         *ii_KernelDigitList;
   int64_t          *il_RemainderList;
   uint64_t         *il_FactorList;
   uint64_t         *il_MagicNumber;
   uint64_t         *il_MagicShift;
   uint32_t          ii_MaxSteps;

   PrimesInXApp     *ip_PrimesInXApp;

   Kernel           *ip_MagicKernel;
   Kernel           *ip_PrimesInXKernel;

   KernelArgument   *ip_KAPrime;
   KernelArgument   *ip_MKAMagicNumber;
   KernelArgument   *ip_MKAMagicShift;

   KernelArgument   *ip_SKAMagicNumber;
   KernelArgument   *ip_SKAMagicShift;
   KernelArgument   *ip_KARemainders;
   KernelArgument   *ip_KAFactors;
   KernelArgument   *ip_KADigits;
};

#endif

