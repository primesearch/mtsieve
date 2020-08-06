/* XYYXApp.h -- (C) Mark Rodenkirch, May 2014

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _XYYXApp_H
#define _XYYXApp_H

#include <vector>

#include "../core/FactorApp.h"

typedef struct {
   uint32_t   exponent;
   uint32_t   base;
   uint32_t   diffPowerIdx;
   uint64_t   fpuRemainders[4];
   double    *avxRemainders;
   void      *pairedPower;
} power_t;

typedef struct {
   uint32_t   base;
   uint32_t   powerCount;
   uint32_t   maxPowerDiff;
   uint32_t   distinctPowers;
   uint8_t   *neededPowers;
   power_t   *powers;
} base_t;

typedef struct {
   base_t    *xPowY;
   base_t    *yPowX;
} bases_t;

class XYYXApp : public FactorApp
{
public:
   XYYXApp(void);

   ~XYYXApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   bool              ReportFactor(uint64_t p, uint32_t x, uint32_t y, int32_t c);
   
   uint32_t          GetMinX(void) { return ii_MinX; };
   uint32_t          GetMaxX(void) { return ii_MaxX; };
   uint32_t          GetMinY(void) { return ii_MinY; };
   uint32_t          GetMaxY(void) { return ii_MaxY; };
   uint32_t          GetMaxPowerDiff(void) { return ii_MaxPowerDiff; };
   uint32_t          GetGpuSteps(void) { return ii_GpuSteps; }
   uint64_t          GetTermCount(void) { return il_TermCount; }
   uint32_t          GetXCount(void) { return (ii_MaxX - ii_MinX + 1); };
   uint32_t          GetYCount(void) { return (ii_MaxY - ii_MinY + 1); };

   bool              UseAvxIfAvailable(void) { return ib_UseAvx; };
   bool              IsPlus(void) { return ib_IsPlus; };
   bool              IsMinus(void) { return ib_IsMinus; };

   void              GetTerms(bool supportsAvx, uint32_t avxRemaindersCount, bases_t *bases);
   
#ifdef HAVE_GPU_WORKERS
   uint32_t          GetNumberOfGroups(void);
   uint32_t          GetGroupedTerms(uint32_t *terms);
#endif

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(void) {};
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);

private:
   void              SetInitialTerms(void);

   vector<bool>      iv_Terms;
   
   bool              ib_UseAvx;
   bool              ib_IsPlus;
   bool              ib_IsMinus;
   uint32_t          ii_MinX;
   uint32_t          ii_MaxX;
   uint32_t          ii_MinY;
   uint32_t          ii_MaxY;
   uint32_t          ii_GpuSteps;
   uint32_t          ii_MaxPowerDiff;
};

#endif

