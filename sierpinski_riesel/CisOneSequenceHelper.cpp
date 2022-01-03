/* CisOneSequenceHelper.cpp -- (C) Mark Rodenkirch, January 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <assert.h>
#include <time.h>
#include "CisOneSequenceHelper.h"
#include "CisOneWithOneSequenceWorker.h"
#include "../core/inline.h"

#define NBIT(n)         ((n) - ii_MinN)

int sortByMapSize(const void *a, const void *b)
{
   legendre_t *aPtr = (legendre_t *) a;
   legendre_t *bPtr = (legendre_t *) b;

   return ( aPtr->bytesNeeded - bPtr->bytesNeeded );
}

int sortByMapSeqIdx(const void *a, const void *b)
{
   legendre_t *aPtr = (legendre_t *) a;
   legendre_t *bPtr = (legendre_t *) b;

   return ( aPtr->seqIdx - bPtr->seqIdx );
}

CisOneSequenceHelper::CisOneSequenceHelper(App *theApp, uint64_t largestPrimeTested) : AbstractSequenceHelper(theApp, largestPrimeTested)
{
}

void  CisOneSequenceHelper::BuildDivisorShifts(void)
{
   uint32_t shifts = ii_PowerResidueLcm / 2;
   uint32_t i, r;
   int32_t  shift, divide;

   ip_DivisorShifts = (int16_t *) xmalloc(shifts * sizeof(int16_t));

   for (i=0; i<shifts; i++)
   {
      r = gcd32(2*i, ii_PowerResidueLcm);

      divide = r;
      for (shift = 0; divide % 2 == 0; shift++)
         divide /= 2;

      if (divide > 1)
      {
         ip_DivisorShifts[i] = r;
         continue;
      }

      if (shift > 1)
         ip_DivisorShifts[i] = -shift;
      else
         ip_DivisorShifts[i] = 0;
   }
}

void  CisOneSequenceHelper::BuildPowerResidueIndices(void)
{
   uint32_t r, prIdx = 1;

   ii_UsedPowerResidueIndices = 0;
   ip_PowerResidueIndices = (uint16_t *) xmalloc((ii_PowerResidueLcm + 1) * sizeof(uint16_t));

   for (r=1; r<ii_PowerResidueLcm; r++)
   {
      if (ii_PowerResidueLcm % r != 0)
         continue;

      ip_PowerResidueIndices[r] = prIdx;
      prIdx++;
   }

   ii_UsedPowerResidueIndices = prIdx;

}

void   CisOneSequenceHelper::CleanUp(void)
{
   xfree(ip_DivisorShifts);
   xfree(ip_PowerResidueIndices);

   if (ip_Legendre != NULL)
      xfree(ip_Legendre);

   if (ip_LegendreTable != NULL)
      xfree(ip_LegendreTable);
}

uint32_t    CisOneSequenceHelper::FindBestQ(uint32_t &expectedSubsequences)
{
   uint32_t          i = 0, j, n;
   uint32_t          bit;
   vector<uint32_t>  S;
   vector<double>    W;
   vector<bool>      R;
   uint32_t          nDivisors = ii_LimitBase / ii_BaseMultiple;
   seq_t            *seqPtr;

   S.resize(nDivisors);
   W.resize(nDivisors);

   R.resize(ii_LimitBase);
      
   seqPtr = ip_FirstSequence;
   do
   {
      std::fill(R.begin(), R.end(), false);
      
      bit = NBIT(ii_MinN);

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seqPtr->nTerms[bit])
            R[n%ii_LimitBase] = true;

         bit++;
      }

      i = 0;
      for (j=0; j<nDivisors; j++)
      {
         if (nDivisors % (j+1) == 0)
            S[j] += CountResidueClasses((j+1)*ii_BaseMultiple, ii_LimitBase, R);
      }

      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
   
   for (i=0, j=0; j<nDivisors; j++)
      if (nDivisors % (j+1) == 0)
      {
         W[j] = RateQ((j+1)*ii_BaseMultiple, S[j]);

         if (W[j] < W[i])
            i = j;
      }

   expectedSubsequences = S[i];

   return (i+1)*ii_BaseMultiple;
}

// Return true iff subsequence h of k*b^n+c has any terms with n%a==b.
bool  CisOneSequenceHelper::HasCongruentTerms(uint32_t ssIdx, uint32_t a, uint32_t b)
{
   uint32_t g, m;

   g = gcd32(a, ii_BestQ);

   if (b % g == ip_Subsequences[ssIdx].q % g)
   {
      for (m=ii_MinM; m<=ii_MaxM; m++)
      {
         if (ip_Subsequences[ssIdx].mTerms[m-ii_MinM])
            if ((m*ii_BestQ + ip_Subsequences[ssIdx].q) % a == b)
               return true;
      }
   }

   return false;
}

void  CisOneSequenceHelper::LastChanceLogicBeforeSieving(void)
{
   uint32_t  ssIdx, babySteps, giantSteps;
   uint32_t  sieveLow = ii_MinN / ii_BestQ;
   uint32_t  sieveHigh = ii_MaxN / ii_BestQ;

   ii_MaxBabySteps = 0;

   for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
   {
      ChooseSteps(ii_BestQ, ssIdx+1, babySteps, giantSteps);

      if (sieveHigh >= sieveLow + babySteps*giantSteps)
         FatalError("LastChanceLogicBeforeSieving miscomputed the steps");

      ip_Subsequences[ssIdx].babySteps = babySteps;
      ip_Subsequences[ssIdx].giantSteps = giantSteps;

      if (ii_MaxBabySteps < babySteps)
         ii_MaxBabySteps = babySteps;
   }

   BuildLegendreTables();
   
   BuildCongruenceTables();
}

void   CisOneSequenceHelper::BuildLegendreTables()
{
   uint64_t     bytesNeeded, bytesUsed;
   uint64_t     stepsDone = 0;
   uint64_t     stepsToDo = 0;
   uint64_t     legendreTableBytes;
   uint32_t     seqNoLegendrePossible = 0;
   uint32_t     seqsWithLegendreMemory = 0;
   uint32_t     seqsWithLegendreFromFile = 0;
   bool         continueAllocating = true;
   seq_t       *seqPtr;
   time_t       startTime, stopTime;
   double       bytes;
   const char  *bytesPrecision;
   legendre_t  *legendrePtr;

   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) ip_App;
   
   legendreTableBytes = srApp->GetLegendreTableBytes();

   ip_Legendre = (legendre_t *) xmalloc(ii_SequenceCount * sizeof(legendre_t));

   seqPtr = ip_FirstSequence;
   bytesNeeded = 0;
   do
   {
      legendrePtr = &ip_Legendre[seqPtr->seqIdx];

      legendrePtr->seqIdx = seqPtr->seqIdx;
      legendrePtr->k = seqPtr->k;
      legendrePtr->c = seqPtr->c;
      legendrePtr->kcCore = seqPtr->kcCore;
      legendrePtr->nParity = seqPtr->nParity;
      legendrePtr->squareFreeK = seqPtr->squareFreeK;
      legendrePtr->canCreateMap = false;
      legendrePtr->haveMap = false;
      legendrePtr->loadedMapFromCache = false;

      ComputeLegendreMemoryToAllocate(legendrePtr, srApp->GetSquareFreeBase());

      if (legendrePtr->canCreateMap)
         bytesNeeded += legendrePtr->bytesNeeded;
      else
         seqNoLegendrePossible++;

      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);

   if (legendreTableBytes == 0)
      return;

   bytes = (double) bytesNeeded;
   bytesPrecision = "B";

   if (bytes >= 10.0 * 1024) { bytes /= 1024.0; bytesPrecision = "KB"; }
   if (bytes >= 10.0 * 1024) { bytes /= 1024.0; bytesPrecision = "MB"; }
   if (bytes >= 10.0 * 1024) { bytes /= 1024.0; bytesPrecision = "GB"; }
   if (bytes >= 10.0 * 1024) { bytes /= 1024.0; bytesPrecision = "TB"; }

   if (legendreTableBytes < bytesNeeded)
      ii_LegendreBytes = legendreTableBytes;
   else
      ii_LegendreBytes = bytesNeeded + 1;    // We are not using offset 0
   
   ip_LegendreTable = (uint8_t *) xmalloc(ii_LegendreBytes * sizeof(uint8_t));
   
   if (ip_LegendreTable == NULL)
   {
      ip_App->WriteToConsole(COT_OTHER, "Approximately %.0f %s needed for Legendre tables", bytes, bytesPrecision);
      FatalError("Not enough memory to allocate space for Legendre tables.  Adjust -l and try again");
   }

   // Sort so that we can allocate by increasing need so that we can build Legendre tables for the most sequences
   qsort(ip_Legendre, ii_SequenceCount, sizeof(legendre_t), sortByMapSize);

   // index = 0 is not used so that the GPU can use != 0 to know that a Legendre table exists
   bytesUsed = 1;
   for (uint32_t legIdx=0; legIdx<ii_SequenceCount; legIdx++)
   {
      legendrePtr = &ip_Legendre[legIdx];

      if (!legendrePtr->canCreateMap)
         continue;

      if (legendrePtr->bytesNeeded + bytesUsed > legendreTableBytes)
      {
         ip_App->WriteToConsole(COT_OTHER, "Stopped building Legendre tables due to limit specified by -l");
         continueAllocating = false;
      }

      if (!continueAllocating)
         continue;

      AssignMemoryToLegendreTable(legendrePtr, bytesUsed);

      seqsWithLegendreMemory++;
      stepsToDo += legendrePtr->mod;
      bytesUsed += legendrePtr->bytesNeeded;
      legendrePtr->haveMap = true;
   }
   
   // Restore the list to the original sequence
   qsort(ip_Legendre, ii_SequenceCount, sizeof(legendre_t), sortByMapSeqIdx);
      
   startTime = time(NULL);

   for (uint32_t legIdx=0; legIdx<ii_SequenceCount; legIdx++)
   {
      legendrePtr = &ip_Legendre[legIdx];

      // If we don't have memory for a Legendre table for this sequence, then there is nothing to do
      if (!legendrePtr->haveMap)
         continue;
      
      LoadLegendreTablesFromFile(legendrePtr);
      
      if (legendrePtr->loadedMapFromCache)
         seqsWithLegendreFromFile++;
      else
      {
         BuildLegendreTableForSequence(legendrePtr, srApp->GetSquareFreeBase(), stepsToDo, stepsDone, startTime);
                  
         WriteLegendreTableToFile(legendrePtr);
      }

      stepsDone += legendrePtr->mod;
   }

   stopTime = time(NULL);

   if (stopTime - startTime > 10)
      ip_App->WriteToConsole(COT_OTHER, "Took %u seconds to build Legendre tables", (uint32_t) (stopTime - startTime));
   
   uint32_t seqLegendreIsPossible = ii_SequenceCount - seqNoLegendrePossible;
      
   ip_App->WriteToConsole(COT_OTHER, "Legendre summary:  Approximately %.0f %s needed for Legendre tables", bytes, bytesPrecision);
   ip_App->WriteToConsole(COT_OTHER, "  %8u total sequences", ii_SequenceCount);
   ip_App->WriteToConsole(COT_OTHER, "  %8u are eligible for Legendre tables", seqLegendreIsPossible);   
   ip_App->WriteToConsole(COT_OTHER, "  %8u are not eligible for Legendre tables", seqNoLegendrePossible);
   ip_App->WriteToConsole(COT_OTHER, "  %8u have Legendre tables in memory", seqsWithLegendreMemory);   
   ip_App->WriteToConsole(COT_OTHER, "  %8u cannot have Legendre tables in memory", ii_SequenceCount - seqsWithLegendreMemory);   
   ip_App->WriteToConsole(COT_OTHER, "  %8u have Legendre tables loaded from files", seqsWithLegendreFromFile);
   ip_App->WriteToConsole(COT_OTHER, "  %8u required building of the Legendre tables", seqsWithLegendreMemory - seqsWithLegendreFromFile);
}

void   CisOneSequenceHelper::LoadLegendreTablesFromFile(legendre_t *legendrePtr)
{
   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) ip_App;
   
   string       directoryName = srApp->GetLegendreDirectoryName();
   char         fileName[500];
   v1_header_t  header;

   if (directoryName.size() == 0)
      return;

#ifdef WIN32
   sprintf(fileName, "%s\\b%u_k%" PRIu64"_c%" PRId64".leg", directoryName.c_str(), ii_Base, legendrePtr->k, legendrePtr->c);
#else
   sprintf(fileName, "%s/b%u_k%" PRIu64"_c%" PRId64".leg", directoryName.c_str(), ii_Base, legendrePtr->k, legendrePtr->c);
#endif
   
   FILE *fPtr = fopen(fileName, "rb");
   
   if (fPtr == NULL)
      return;
   
   if (fread(&header, sizeof(v1_header_t), 1, fPtr) < 1)
   {
      ip_App->WriteToConsole(COT_OTHER, "Could not read header from Legendre file %s", fileName);
      fclose(fPtr);
      return;
   }
   
   if (!ValidateLegendreFile(&header, legendrePtr))
   {
      fclose(fPtr);
      return;
   }

   if (!ReadLegendreTableFromFile(fPtr, legendrePtr->oneParityMap,    legendrePtr->mapSize, header.oneParityMapIndex) &&
       !ReadLegendreTableFromFile(fPtr, legendrePtr->dualParityMapM1, legendrePtr->mapSize, header.dualParityMapM1Index) &&
       !ReadLegendreTableFromFile(fPtr, legendrePtr->dualParityMapP1, legendrePtr->mapSize, header.dualParityMapP1Index))
   {
      fclose(fPtr);
      return;
   }
   
   legendrePtr->loadedMapFromCache = true;

   fclose(fPtr);
}

bool   CisOneSequenceHelper::ValidateLegendreFile(v1_header_t *headerPtr, legendre_t *legendrePtr)
{
   if (headerPtr->base != ii_Base)
   {
      ip_App->WriteToConsole(COT_OTHER, "base in Legendre file is not the expected base (%u != %u)", headerPtr->base, ii_Base);
      return false;
   }
   
   if (headerPtr->fileVersion != 1)
   {
      ip_App->WriteToConsole(COT_OTHER, "version %u in Legendre file is not supported", 1);
      return false;
   }

   if (headerPtr->structSize != sizeof(v1_header_t))
   {
      ip_App->WriteToConsole(COT_OTHER, "header size in Legendre file is not the expected header size (%u != %u)", headerPtr->structSize, (uint32_t) sizeof(v1_header_t));   
      return false;
   }
   
   if (headerPtr->k != legendrePtr->k)
   {
      ip_App->WriteToConsole(COT_OTHER, "k in Legendre file is not the expected k (%" PRIu64" != %" PRIu64")", headerPtr->k, legendrePtr->k);  
      return false;
   }
   
   if (headerPtr->c != legendrePtr->c)
   {
      ip_App->WriteToConsole(COT_OTHER, "c in Legendre file is not the expected c (%" PRId64" != %" PRId64")", headerPtr->c, legendrePtr->c);  
      return false;
   }
   
   if (headerPtr->mapSize != legendrePtr->mapSize)
   {
      ip_App->WriteToConsole(COT_OTHER, "mapSize for k=%" PRIu64" c=%" PRId64" is not the expected mapSize (%u != %u)",
               legendrePtr->k, legendrePtr->c, headerPtr->mapSize, legendrePtr->mapSize);
      return false;
   }

   if (headerPtr->oneParityMapIndex    == 0 &&
       headerPtr->dualParityMapM1Index == 0 &&
       headerPtr->dualParityMapP1Index == 0)
   {
      ip_App->WriteToConsole(COT_OTHER, "no maps are in the file for k=%" PRIu64" c=%" PRId64"", legendrePtr->k, legendrePtr->c);
      return false;
   }
   
   if ((headerPtr->oneParityMapIndex    == 0 && legendrePtr->oneParityMapIndex    != 0) ||
       (headerPtr->oneParityMapIndex    != 0 && legendrePtr->oneParityMapIndex    == 0) ||
       (headerPtr->dualParityMapM1Index == 0 && legendrePtr->dualParityMapM1Index != 0) ||
       (headerPtr->dualParityMapM1Index != 0 && legendrePtr->dualParityMapM1Index == 0) ||
       (headerPtr->dualParityMapP1Index == 0 && legendrePtr->dualParityMapP1Index != 0) ||
       (headerPtr->dualParityMapP1Index != 0 && legendrePtr->dualParityMapP1Index == 0))
   {
      ip_App->WriteToConsole(COT_OTHER, "wrong maps for k=%" PRIu64" c=%" PRId64"", legendrePtr->k, legendrePtr->c);
      return false;
   }

   return true;
}

bool   CisOneSequenceHelper::ReadLegendreTableFromFile(FILE *fPtr, uint8_t *map, uint32_t mapSize, uint64_t offset)
{
   // No map to read so assume success
   if (map == NULL)
      return true;
      
   if (fseek(fPtr, offset, SEEK_SET) != 0)
   {
      ip_App->WriteToConsole(COT_OTHER, "fseek to offset %" PRIu64" failed", offset);
      return false;
   }
   
   if (fread(map, sizeof(uint8_t), mapSize, fPtr) < mapSize)
   {
      ip_App->WriteToConsole(COT_OTHER, "could not read map from file");
      return false;
   }
   
   return true;
}

void   CisOneSequenceHelper::WriteLegendreTableToFile(legendre_t *legendrePtr)
{
   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) ip_App;
   
   string       directoryName = srApp->GetLegendreDirectoryName();
   char         fileName[500];
   uint32_t     offset;
   v1_header_t  header;

   if (directoryName.size() == 0)
      return;

#ifdef WIN32
   sprintf(fileName, "%s\\b%u_k%" PRIu64"_c%" PRId64".leg", directoryName.c_str(), ii_Base, legendrePtr->k, legendrePtr->c);
#else
   sprintf(fileName, "%s/b%u_k%" PRIu64"_c%" PRId64".leg", directoryName.c_str(), ii_Base, legendrePtr->k, legendrePtr->c);
#endif
      
   FILE *fPtr = fopen(fileName, "wb");
   
   if (fPtr == NULL)
      FatalError("Could not open Legendre file %s", fileName);
   
   offset = sizeof(v1_header_t);
   
   header.fileVersion = 1;
   header.structSize = sizeof(v1_header_t);
   header.base = ii_Base;
   header.k = legendrePtr->k;
   header.c = legendrePtr->c;
   header.mapSize = legendrePtr->mapSize;
   header.oneParityMapIndex = 0;
   header.dualParityMapM1Index = 0;
   header.dualParityMapP1Index = 0;
   
   if (legendrePtr->oneParityMap != NULL)
   {
      header.oneParityMapIndex = offset;
      offset += legendrePtr->mapSize;
   }
   
   if (legendrePtr->dualParityMapM1 != NULL)
   {
      header.dualParityMapM1Index = offset;
      offset += legendrePtr->mapSize;
   }
   
   if (legendrePtr->dualParityMapP1 != NULL)
   {
      header.dualParityMapP1Index = offset;
      offset += legendrePtr->mapSize;
   }
         
   if (fwrite(&header, sizeof(v1_header_t), 1, fPtr) < 1)
      FatalError("could not write header to file");
   
   if (legendrePtr->oneParityMap != NULL)
   {
      if (fwrite(legendrePtr->oneParityMap, sizeof(uint8_t), legendrePtr->mapSize, fPtr) < legendrePtr->mapSize)
         FatalError("could not write map to file");
   }

   if (legendrePtr->dualParityMapM1 != NULL)
   {
      if (fwrite(legendrePtr->dualParityMapM1, sizeof(uint8_t), legendrePtr->mapSize, fPtr) < legendrePtr->mapSize)
         FatalError("could not write map to file");
   }
   
   if (legendrePtr->dualParityMapP1 != NULL)
   {
      if (fwrite(legendrePtr->dualParityMapP1, sizeof(uint8_t), legendrePtr->mapSize, fPtr) < legendrePtr->mapSize)
         FatalError("could not write map to file");
   }

   fclose(fPtr);
}
