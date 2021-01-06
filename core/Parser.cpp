/* Parser.cpp -- (C) Mark Rodenkirch, February 2012

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdlib.h>

#include "Parser.h"
#include "main.h"

Parser::Parser(void)
{
}

Parser::~Parser()
{
}

parse_t  Parser::Parse(const char *str, const char *valid, char &value)
{
   while (*valid)
   {
	   if (*str == *valid)
      {
         value = *str;
         return P_SUCCESS;
      }

      valid++;
   }

   return P_FAILURE;
}

parse_t  Parser::Parse(const char *str, int32_t lo, int32_t hi, int32_t &value)
{
   int64_t result64;
   parse_t status;

   status = Parser::Parse(str, lo, hi, result64);

   if (status == P_SUCCESS)
      value = (int32_t) result64;

   return status;
}

parse_t  Parser::Parse(const char *str, uint32_t lo, uint32_t hi, uint32_t &value)
{
   uint64_t result64;
   parse_t status;

   status = Parser::Parse(str, lo, hi, result64);

   if (status == P_SUCCESS)
      value = (uint32_t) result64;

   return status;
}

parse_t  Parser::Parse(const char *str, int64_t lo, int64_t hi, int64_t &value)
{
   uint64_t unsigned_result64;
   int64_t  signed_result64;
   parse_t status;

   bool negative = false;

   if (*str == '-')
   {
      negative = true;
      str++;
   }

   status = Parser::Parse(str, 0, INT64_MAX, unsigned_result64);

   if (status != P_SUCCESS)
      return status;

   signed_result64 = (int64_t) unsigned_result64;

   if (negative)
      signed_result64 *= -1;
   
   if (signed_result64 > hi)
      return P_OUT_OF_RANGE;

   if (signed_result64 < lo)
      return P_OUT_OF_RANGE;

   value = signed_result64;
   return P_SUCCESS;
}

parse_t  Parser::Parse(const char *str, uint64_t lo, uint64_t hi, uint64_t &value)
{
   uint64_t num;
   uint32_t expt;
   char *tail;

   if (*str == '-')
      return P_OUT_OF_RANGE;

   expt = 0;
   errno = 0;
   num = strtoull(str, &tail, 0);

   if (errno != 0 || num > hi)
      return P_OUT_OF_RANGE;

   switch (*tail)
   {
      case 'P': expt += 3;
      case 'T': expt += 3;
      case 'G': expt += 3;
      case 'M': expt += 3;
      case 'K': expt += 3;
         if (tail[1] != '\0')
            return P_FAILURE;
         for ( ; expt > 0; expt -= 3)
            if (num > hi/1000)
               return P_OUT_OF_RANGE;
            else
               num *= 1000;
         break;

      case 'e':
      case 'E':
         expt = strtoul(tail+1,&tail,0);
         if (errno != 0)
            return P_OUT_OF_RANGE;
         if (*tail != '\0')
            return P_FAILURE;
         while (expt-- > 0)
            if (num > hi/10)
               return P_OUT_OF_RANGE;
            else
               num *= 10;
         break;

      case 'p': expt += 10;
      case 't': expt += 10;
      case 'g': expt += 10;
      case 'm': expt += 10;
      case 'k': expt += 10;
         if (tail[1] != '\0')
            return P_FAILURE;
         if (num > (hi>>expt))
            return P_OUT_OF_RANGE;
         num <<= expt;
         break;

      case 'b':
      case 'B':
         expt = strtoul(tail+1,&tail,0);
         if (errno != 0)
            return P_OUT_OF_RANGE;
         if (*tail != '\0')
            return P_FAILURE;
         while (expt-- > 0)
            if (num > (hi>>1))
               return P_OUT_OF_RANGE;
            else
               num <<= 1;
         break;

      case '\0':
         break;

      default:
         return P_FAILURE;
   }

   if (num < lo)
      return P_OUT_OF_RANGE;

   value = num;
   return P_SUCCESS;
}

parse_t  Parser::Parse(const char *str, double lo, double hi, double &value)
{
   double num;

   if (*str == '-')
      return P_OUT_OF_RANGE;

   if (sscanf(str, "%lf", &num) != 1)
      return P_FAILURE;
   
   if (num < lo || num > hi)
      return P_OUT_OF_RANGE;
   
   value = num;
   return P_SUCCESS;
}