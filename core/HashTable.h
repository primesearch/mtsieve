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
#define HASH_MAX_ELTS      HASH_MASK2
#define HASH_MINIMUM_SHIFT 10

class HashTable
{
public:
   HashTable(uint32_t elements);
   ~HashTable(void);

   void  Cleanup(void);
   
   inline void Clear(void)
   {
      memset(htable, 0, hsize*sizeof(uint16_t));
   }
   
   inline void Insert(uint64_t bj, uint32_t j)
   {
      uint32_t slot;

      BJ64[j] = bj;
      slot = bj & hsize_minus1;
      if (htable[slot] == empty_slot)
         htable[slot] = j;
      else
      {
         olist[j] = (htable[slot] ^ HASH_MASK1);
         htable[slot] = (j | HASH_MASK1);
      }
   };
   
   inline uint32_t Lookup(uint64_t bj)
   {
      uint32_t slot;
      uint16_t elt;

      slot = bj & hsize_minus1;
      elt = htable[slot];

      if (BJ64[elt & HASH_MASK2] == bj)
         return elt & HASH_MASK2;
      
      if ((elt & HASH_MASK1) == 0)
         return HASH_NOT_FOUND;

      elt &= HASH_MASK2;
      do
      {
         if (BJ64[olist[elt] & HASH_MASK2] == bj)
            return olist[elt] & HASH_MASK2;
         
         elt = olist[elt];
      } while ((elt & HASH_MASK1) == 0);

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
