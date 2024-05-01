/* LifchitzGpuWorker.h -- (C) Mark Rodenkirch, April 2024

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _LifchitzGpuWorker_H
#define _LifchitzGpuWorker_H

#include "LifchitzApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

// This maps to a uint4 in the GPU
typedef struct {
   uint32_t    x;
   uint32_t    y;
   uint32_t    sign;
   uint32_t    unused;
} gputerm_t;

// This maps to a ulong4 in the GPU
typedef struct {
   uint64_t    termIdx;
   uint64_t    sign;
   uint64_t    factor;
   uint64_t    unused;
} factor_t;


class LifchitzGpuWorker : public Worker
{
public:
   LifchitzGpuWorker(uint32_t myId, App *theApp);

   ~LifchitzGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList) {} 

   void              AllocateMemory(void);
   void              CreateTermGroups(void);

   uint64_t          il_NextTermsBuild;
   
   uint32_t          ii_MinX;
   uint32_t          ii_MaxX;
   uint32_t          ii_MinY;
   uint32_t          ii_MaxY;
   
   uint32_t          ii_XChunks;
   uint32_t          ii_YChunks;
   
   uint32_t          ii_MaxBases;
   
   uint32_t          ii_XBasesPerGroup;
   uint32_t          ii_YBasesPerGroup;
   uint32_t          ii_TermsPerGroup;
   
   uint32_t        **ip_XBases;
   uint32_t        **ip_YBases;
   uint32_t        **ip_TermsInGroup;
   gputerm_t      ***ip_GpuTerms;
   uint64_t       ***ip_TermIndexes;

   uint32_t          ii_MaxGpuFactors;
   
   uint32_t         *ii_KernelXBases;
   uint32_t         *ii_KernelYBases;
   gputerm_t        *ip_TermList;
   uint32_t         *ii_FactorCount;
   factor_t         *ip_FactorList;

   LifchitzApp      *ip_LifchitzApp;

   GpuKernel        *ip_Kernel;
};

#endif

