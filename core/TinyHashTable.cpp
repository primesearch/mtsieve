/* TinyHashTable.cpp -- (C) Mark Rodenkirch, May 2023

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <inttypes.h>
#include <assert.h>
#include "main.h"
#include "TinyHashTable.h"

#define HASH_MINIMUM_SHIFT    3

#define HASH_MAX_DENSITY      0.60
#define HASH_MINIMUM_ELTS     8

TinyHashTable::TinyHashTable(uint32_t elements) : HashTable()
{
   assert(elements <= TINY_HASH_MAX_ELTS);

   if (elements < HASH_MINIMUM_ELTS)
      elements = HASH_MINIMUM_ELTS;

   for (hsize = 1<<HASH_MINIMUM_SHIFT; hsize < elements/HASH_MAX_DENSITY; )
      hsize *= 2;

   hsize_minus1 = hsize - 1;

   htable = (uint8_t *) xmalloc(hsize*sizeof(uint8_t));
   olist = (uint8_t *) xmalloc(elements*sizeof(uint8_t));

   // The j values are all in the range 0 <= j < M, so we can use M as an
   // empty slot marker as long as we fill BJ[M] with a value that will never
   // match a real b^j value. Since b^j is always in the range 0 <= b^j < p
   // for some prime p, any value larger than all 32/64 bit primes will do.
   empty_slot = elements;
   BJ64 = (uint64_t *) xmalloc((elements+1)*sizeof(uint64_t));
   BJ64[empty_slot] = UINT64_MAX;
   
   Clear();
}

TinyHashTable::~TinyHashTable(void)
{
   xfree(BJ64);
   xfree(olist);
   xfree(htable);
}

