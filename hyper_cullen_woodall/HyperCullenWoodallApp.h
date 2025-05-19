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
   
   bool              ReportFactor(uint64_t theFactor, uint32_t b, uint32_t n, int32_t sign);
   
   uint32_t          GetMinB(void) { return ii_MinB; };
   uint32_t          GetMaxB(void) { return ii_MaxB; };
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };
   uint64_t          GetTermCount(void) { return il_TermCount; }
   uint32_t          GetBCount(void) { return (ii_MaxB - ii_MinB + 1); };
   uint32_t          GetNCount(void) { return (ii_MaxN - ii_MinN + 1); };

   void              BuildTerms(void);
   
   bool              IsPlus(void) { return ib_IsPlus; };
   bool              IsMinus(void) { return ib_IsMinus; };

   base_t           *GetBTerms(void) { return ip_bTerms; };
   base_t           *GetNTerms(void) { return ip_nTerms; };
   
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
   void              VerifyFactor(uint64_t theFactor, uint32_t b, uint32_t n, int32_t sign);
   
   std::vector<bool> iv_PlusTerms;
   std::vector<bool> iv_MinusTerms;
      
   bool              ib_IsPlus;
   bool              ib_IsMinus;
   uint32_t          ii_MinB;
   uint32_t          ii_MaxB;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;

   bool              ib_SplitTerms;
   uint32_t          ii_SplitB;
   uint32_t          ii_SplitN;
   uint64_t          il_SplitLargestPrimeTested;

   base_t           *ip_bTerms;
   base_t           *ip_nTerms;

#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
#endif
};

#endif

