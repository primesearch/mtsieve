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
#include "../opencl/Kernel.h"
#include "../opencl/KernelArgument.h"

using namespace std;

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
   void              ValidateFactor(int64_t p, uint32_t term);
   
   AlternatingFactorialApp      *ip_AlternatingFactorialApp;

   int64_t          *il_RemainderList;
   int64_t          *il_NM1TermList;
   uint64_t         *il_MagicNumber;
   uint64_t         *il_MagicShift;
   uint64_t         *il_FactorList;
   int32_t           ii_NRange[10];
   uint32_t          ii_MaxFactors;

   Kernel           *ip_MagicKernel;
   Kernel           *ip_AlternatingFactorialKernel;

   KernelArgument   *ip_KAPrime;
   KernelArgument   *ip_MKAMagicNumber;
   KernelArgument   *ip_MKAMagicShift;

   KernelArgument   *ip_FKAMagicNumber;
   KernelArgument   *ip_FKAMagicShift;
   KernelArgument   *ip_KANRange;
   KernelArgument   *ip_KARemainder;
   KernelArgument   *ip_KANM1Term;
   KernelArgument   *ip_KAFactor;
};

#endif

