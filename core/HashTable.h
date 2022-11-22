/* HashTable.h -- (C) Mark Rodenkirch, January 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _HASHTABLE_H
#define _HASHTABLE_H

#define HASH_NOT_FOUND     UINT32_MAX
#define HASH_MASK1         (1<<15)
#define HASH_MASK2         (HASH_MASK1-1)
#define HASH_MAX_ELTS      HASH_MASK2-1
#define HASH_MINIMUM_SHIFT 11

class HashTable
{
public:
   HashTable(uint32_t elements);
   
   ~HashTable(void);

   void  Cleanup(void);
   
   uint64_t GetInserts(void) { return il_Inserts; };
   uint64_t GetConflicts(void) { return il_Conflicts; };
   
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
      htable[slot] = (j | HASH_MASK1);
      il_Conflicts++;

   };
   
   inline uint32_t Lookup(uint64_t bj)
   {
      uint32_t slot;
      uint16_t elt;
      uint16_t elt_low;

      slot = bj & hsize_minus1;
      elt = htable[slot];
      elt_low = elt & HASH_MASK2;

      // Could check elt == empty_slot to avoid checking BJ64
      if (BJ64[elt_low] == bj)
         return elt_low;
      
      while (elt != elt_low) {
         elt = olist[elt & HASH_MASK2];
         elt_low = elt & HASH_MASK2;
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
   
   uint64_t  il_Inserts;
   uint64_t  il_Conflicts;
};

#endif
