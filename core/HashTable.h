/* HashTable.h -- (C) Mark Rodenkirch, January 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _HASHTABLE_H
#define _HASHTABLE_H

#include <inttypes.h>

#define HASH_NOT_FOUND     UINT32_MAX

class HashTable
{
public:
   HashTable(void);
   
   virtual ~HashTable(void) {};
   
   uint64_t GetInserts(void) { return il_Inserts; };
   uint64_t GetConflicts(void) { return il_Conflicts; };
   
   virtual void Clear(void) = 0;
   
   virtual uint64_t get(uint32_t x) = 0;
   
   virtual void  Insert(uint64_t bj, uint32_t j) = 0;
   
   virtual uint32_t Lookup(uint64_t bj) = 0;

protected:  
   uint64_t  il_Inserts;
   uint64_t  il_Conflicts;
};

#endif
