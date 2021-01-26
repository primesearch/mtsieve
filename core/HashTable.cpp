/* HashTable.cpp -- (C) Mark Rodenkirch, January 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <inttypes.h>
#include <assert.h>
#include "main.h"
#include "HashTable.h"

#define HASH_MAX_DENSITY      0.60
#define HASH_MIN_DENSITY      0.10

#define HASH_MINIMUM_ELTS     8
#define L1_CACHE_SIZE         32

HashTable::HashTable(uint32_t elements)
{
   assert(elements <= HASH_MAX_ELTS);

   if (elements < HASH_MINIMUM_ELTS)
      elements = HASH_MINIMUM_ELTS;

   for (hsize = 1<<HASH_MINIMUM_SHIFT; hsize < elements/HASH_MAX_DENSITY; )
      hsize <<= 1;

   hsize_minus1 = hsize - 1;

   htable = (uint16_t *) xmalloc(hsize*sizeof(uint16_t));
   olist = (uint16_t *) xmalloc(elements*sizeof(uint16_t));

   // The j values are all in the range 0 <= j < M, so we can use M as an
   // empty slot marker as long as we fill BJ[M] with a value that will never
   // match a real b^j value. Since b^j is always in the range 0 <= b^j < p
   // for some prime p, any value larger than all 32/64 bit primes will do.
   empty_slot = elements;
   BJ64 = (uint64_t *) xmalloc((elements+1)*sizeof(uint64_t));
   BJ64[empty_slot] = UINT64_MAX;
}

HashTable::HashTable(uint32_t elements, uint32_t bestQ, uint32_t powerResidueLcm)
{
   uint32_t bjLen;
   
   assert(elements <= HASH_MAX_ELTS);

   if (elements < HASH_MINIMUM_ELTS)
      elements = HASH_MINIMUM_ELTS;

   for (hsize = 1<<HASH_MINIMUM_SHIFT; hsize*HASH_MAX_DENSITY < elements; )
      hsize *= 2;

   // Double size of hashtable if it will fit into L1 cache
   if (hsize*sizeof(uint16_t) < L1_CACHE_SIZE*1024 && hsize <= elements/HASH_MIN_DENSITY/2)
      hsize *= 2;
      
   if (hsize*sizeof(uint16_t)*2 < L1_CACHE_SIZE*1024 && hsize <= elements/HASH_MIN_DENSITY/2)
      hsize *= 2;

   // Don't exceed maximum hash
   if (hsize > HASH_MASK1)
      hsize = HASH_MASK1;

   hsize_minus1 = hsize - 1;

   htable = (uint16_t *) xmalloc(hsize*sizeof(uint16_t));
   olist = (uint16_t *) xmalloc(elements*sizeof(uint16_t));

   bjLen = bestQ+powerResidueLcm+9;

   // The j values are all in the range 0 <= j < M, so we can use M as an
   // empty slot marker as long as we fill BJ[M] with a value that will never
   // match a real b^j value. Since b^j is always in the range 0 <= b^j < p
   // for some prime p, any value larger than all 32/64 bit primes will do.
   empty_slot = elements+1;
   BJ64 = (uint64_t *) xmalloc(bjLen*sizeof(uint64_t));
   BJ64[empty_slot] = UINT64_MAX;
}

HashTable::~HashTable(void)
{
   xfree(BJ64);
   xfree(olist);
   xfree(htable);
}

