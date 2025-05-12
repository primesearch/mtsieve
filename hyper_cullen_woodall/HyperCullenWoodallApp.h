/* HyperCullenWoodallApp.h -- (C) Mark Rodenkirch
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _HyperCullenWoodallApp_H
#define _HyperCullenWoodallApp_H

#include <vector>
#include <set>

#include "../core/MpArithVector.h"
#include "../core/FactorApp.h"

typedef struct {
   uint32_t                base;
   uint32_t                powerCount;
   bool                    indexedByPower;
   std::vector<uint32_t>   powers;
   std::vector<MpResVec>   residues;
} base_t;

typedef struct {
   uint32_t   *terms;
   void       *nextGroup;
} gpugroup_t;

class HyperCullenWoodallApp : public FactorApp
{
public:
   HyperCullenWoodallApp(void);

   ~HyperCullenWoodallApp(void);

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength);
   
   bool              ReportFactor(uint64_t theFactor, uint32_t x, uint32_t y, int32_t sign);
   
   uint32_t          GetMinX(void) { return ii_MinX; };
   uint32_t          GetMaxX(void) { return ii_MaxX; };
   uint32_t          GetMinY(void) { return ii_MinY; };
   uint32_t          GetMaxY(void) { return ii_MaxY; };
   uint64_t          GetTermCount(void) { return il_TermCount; }
   uint32_t          GetXCount(void) { return (ii_MaxX - ii_MinX + 1); };
   uint32_t          GetYCount(void) { return (ii_MaxY - ii_MinY + 1); };

   void              BuildTerms(void);
   
   bool              IsPlus(void) { return ib_IsPlus; };
   bool              IsMinus(void) { return ib_IsMinus; };

   base_t           *GetXTerms(void) { return ip_xTerms; };
   base_t           *GetYTerms(void) { return ip_yTerms; };
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   void              FreeGpuGroupsOfTerms(gpugroup_t *gPtr);
   gpugroup_t       *GetGpuGroupsOfTerms(void);
   
   uint32_t          GetMaxGpuSteps(void) { return ii_MaxGpuSteps; };
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
#endif

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   void              OuptutAdditionalConsoleMessagesUponFinish(void) {};

private:
   void              SetInitialTerms(void);
   void              VerifyFactor(uint64_t theFactor, uint32_t x, uint32_t y, int32_t sign);
   
   std::vector<bool> iv_PlusTerms;
   std::vector<bool> iv_MinusTerms;
      
   bool              ib_IsPlus;
   bool              ib_IsMinus;
   uint32_t          ii_MinX;
   uint32_t          ii_MaxX;
   uint32_t          ii_MinY;
   uint32_t          ii_MaxY;

   base_t           *ip_xTerms;
   base_t           *ip_yTerms;

#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
#endif
};

#endif

