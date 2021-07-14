/* AlgebraicFactorHelper.h -- (C) Mark Rodenkirch, December 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _AlgebraicFactorHelper_H
#define _AlgebraicFactorHelper_H

#include <stdio.h>
#include "SierpinskiRieselApp.h"

class AlgebraicFactorHelper
{  
public:
   AlgebraicFactorHelper(App *theApp, uint32_t base, uint32_t minN, uint32_t maxN);
   ~AlgebraicFactorHelper(void);

   // This can only search for all algebraic factors for a single k/c per call.
   uint64_t          RemoveTermsWithAlgebraicFactors(seq_t *seq);
   
private:
   App              *ip_App;
   
   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint32_t          ii_BRoot;
   uint32_t          ii_BPower;
   uint32_t          ii_KRoot;
   uint32_t          ii_KPower;
   
   // These variables are used by this class when looking for algebraic factors.
   vector<uint64_t>  iv_SmallPrimes;
   
   FILE             *ip_AlgebraicFactorFile;

   void              GetSmallPrimes(void);
   uint32_t          GetFactorList(uint64_t the_number, uint32_t *factor_list, uint32_t *power_list);
   
   void              CheckForSpecialForm(seq_t *seq);
   uint64_t          RemoveSimpleTerm(seq_t *seq);
   uint64_t          RemoveTermsWithKPowers(seq_t *seq);
   uint64_t          RemoveTermsWithKAndBPowers(seq_t *seq);
   uint64_t          RemoveComplexRoot(seq_t *seq);
   uint64_t          CheckBase2(seq_t *seq);
   uint64_t          CheckPower4(seq_t *seq);
   
#ifdef __MINGW_PRINTF_FORMAT
   uint32_t          CheckAndLogAlgebraicFactor(seq_t *seq, uint32_t n, const char *fmt, ...) __attribute__ ((format (__MINGW_PRINTF_FORMAT, 4, 5)));
#else
   uint32_t          CheckAndLogAlgebraicFactor(seq_t *seq, uint32_t n, const char *fmt, ...) __attribute__ ((format (printf, 4, 5)));
#endif

   void              GetRoot(uint64_t number, uint32_t *root, uint32_t *power);
   bool              IsGfnOrMersenneForm(seq_t *seq);
};

#endif
