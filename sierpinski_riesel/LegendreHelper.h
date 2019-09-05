/* LegendreHelper.h -- (C) Mark Rodenkirch, July 2019
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This class is used if any sequence has abs(c) = 1.
*/

#ifndef _LegendreHelper_H
#define _LegendreHelper_H

#include "AbstractSubsequenceHelper.h"

class LegendreHelper
{
public:
   LegendreHelper(uint32_t base, seq_t *sequences, uint32_t sequenceCount);

   ~LegendreHelper(void) {};

   bool              LoadLegendreTables(string legendreFileName);

   bool              BuildLegendreTables(void);
   
private:
   uint64_t          SquareFreeFactorization(uint64_t n);
   
   seq_t            *ip_Sequences;
   
   vector<uint64_t>  iv_Primes;
   
   uint32_t          ii_SequenceCount;
   uint32_t          ii_Base;
};

#endif

