/* SmallHashTable.h -- (C) Mark Rodenkirch, May 2023

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SmallHashTable_H
#define _SmallHashTable_H

#include "HashTable.h"

#define SMALL_HASH_MASK1         (1<<15)
#define SMALL_HASH_MASK2         (SMALL_HASH_MASK1-1)
#define SMALL_HASH_MAX_ELTS      SMALL_HASH_MASK2-1

class SmallHashTable : public HashTable
{
public:
   SmallHashTable(uint32_t elements);
   
   ~SmallHashTable(void);

   inline void Clear(void)
   {
      // hsize is always a power of 2
      for (uint32_t i=0; i < hsize; i+=2)
      {
         htable[i] = empty_slot;
         htable[i+1] = empty_slot;
      }
      
      // will never match a Lookup value
      BJ64[empty_slot] = UINT64_MAX;
   }
   
   inline uint64_t get(uint32_t x) { return BJ64[x]; };
   
   inline void  Insert(uint64_t bj, uint32_t j)
   {
      uint32_t slot;
      il_Inserts++;

      BJ64[j] = bj;
      slot = bj & hsize_minus1;
      if (htable[slot] == empty_slot)
      {
         htable[slot] = j;
         return;
      }

      olist[j] = htable[slot];
      htable[slot] = (j | SMALL_HASH_MASK1);
      il_Conflicts++;

   };
   
   inline uint32_t Lookup(uint64_t bj)
   {
      uint32_t slot;
      uint16_t elt;
      uint16_t elt_low;

      slot = bj & hsize_minus1;
      elt = htable[slot];
      elt_low = elt & SMALL_HASH_MASK2;

      // Could check elt == empty_slot to avoid checking BJ64
      if (BJ64[elt_low] == bj)
         return elt_low;
      
      while (elt != elt_low) {
         elt = olist[elt & SMALL_HASH_MASK2];
         elt_low = elt & SMALL_HASH_MASK2;
         if (BJ64[elt_low] == bj)
            return elt_low;
      }
      
      return HASH_NOT_FOUND;
   };

private:
   uint16_t  empty_slot;
   uint16_t *htable;
   uint16_t *olist;
   uint32_t  hsize;
   uint32_t  hsize_minus1;
   uint64_t *BJ64;
};

#endif
