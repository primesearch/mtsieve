#ifdef __cplusplus
#ifndef _CISONEMULTIPLE_KERNEL_CL
#define _CISONEMULTIPLE_KERNEL_CL
const char *cisonemultiple_kernel= \
"#if !defined(USE_OPENCL) && !defined(USE_METAL)\n" \
"#define MEED_METALLIB\n" \
"#endif\n" \
"#ifdef USE_METAL\n" \
"#include <metal_stdlib>\n" \
"#include <metal_atomic>\n" \
"using namespace metal;\n" \
"#endif\n" \
"#ifdef MEED_METALLIB\n" \
"#define MAX_FACTORS           500\n" \
"#define BASE                    2\n" \
"#define BEST_Q                  5\n" \
"#define SIEVE_LOW               5\n" \
"#define SEQUENCES             100\n" \
"#define SUBSEQUENCES          100\n" \
"#define HASH_ELEMENTS         500\n" \
"#define HASH_SIZE             256\n" \
"#define HASH_SIZE_M1          255\n" \
"#define POWER_RESIDUE_LCM     720\n" \
"#define DIM2                  120\n" \
"#define DIM3                  120\n" \
"#endif\n" \
"#if defined(USE_OPENCL)\n" \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"#define MUL_HI          mul_hi\n" \
"void collectFactor(ulong   p,\n" \
"ulong    ssIdx,\n" \
"ulong    n,\n" \
"ulong    x,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors);\n" \
"#else\n" \
"#define MUL_HI          mulhi\n" \
"void collectFactor(ulong       p,\n" \
"uint        ssIdx,\n" \
"uint        n,\n" \
"volatile device atomic_int *factorCount,\n" \
"device ulong4     *factors);\n" \
"#endif\n" \
"#define SP_MIXED     0\n" \
"#define SP_EVEN      1\n" \
"#define SP_ODD       2\n" \
"#define SP_NO_PARITY 999\n" \
"#define MM_P(p)             (p.x)\n" \
"#define MM_Q(p)             (p.y)\n" \
"#define MM_ONE(p)           (p.z)\n" \
"#define MM_R2(p)            (p.w)\n" \
"#define SEQ_SEQIDX(a,b)     (a[b].x)\n" \
"#define SEQ_PARITY(a,b)     (a[b].y)\n" \
"#define SEQ_MINSS(a,b)      (a[b].z)\n" \
"#define SEQ_MAXSS(a,b)      (a[b].w)\n" \
"#define SUBSEQ_Q(a,b)       (a[b].x)\n" \
"#define SUBSEQ_BS(a,b)      (a[b].y)\n" \
"#define SUBSEQ_GS(a,b)      (a[b].z)\n" \
"#define USABLE_SEQIDX(a,b)  (a[b].x)\n" \
"#define USABLE_SSIDX(a,b)   (a[b].y)\n" \
"#define USABLE_Q(a,b)       (a[b].z)\n" \
"#define USABLE_BDCK(a,b)    (a[b].w)\n" \
"#define N_TERM(q, i, j)       ((SIEVE_LOW + (j) + (i)*babySteps)*BEST_Q + q)\n" \
"#define CSS_INDEX(a, b, c) ((a) * DIM2 + (b) * DIM3 + (c))\n" \
"#define L_BYTE(x)  ((x)>>3)\n" \
"#define L_BIT(x)   (1<<((x)&7))\n" \
"#define HASH_NOT_FOUND        UINT_MAX\n" \
"#define HASH_MASK1            (1<<15)\n" \
"#define HASH_MASK2            (HASH_MASK1-1)\n" \
"ulong  invmod(ulong a, ulong p);\n" \
"short  legendre(long a, ulong p);\n" \
"ulong  getNegCK(ulong k, int c, ulong p);\n" \
"uint      setupDiscreteLog(__global const ulong   *ks,\n" \
"__global const long    *kcCores,\n" \
"__global const uint4   *seqData,\n" \
"__global const uint4   *subseqData,\n" \
"__global const short   *divisorShifts,\n" \
"__global const ushort  *prlIndices,\n" \
"__global const uint    *cssIndices,\n" \
"__global const uint    *cssSubseqs,\n" \
"ulong   *resBD,\n" \
"ulong4   mmPrime,\n" \
"ulong    resInvBase,\n" \
"ulong   *resX,\n" \
"ulong4  *usableSubsequence);\n" \
"uint      getShift0Subsequences(__global const ulong   *ks,\n" \
"__global const long    *kcCores,\n" \
"__global const uint4   *seqData,\n" \
"__global const uint4   *subseqData,\n" \
"ulong   *resBD,\n" \
"ulong4   mmPrime,\n" \
"ulong4  *usableSubsequence);\n" \
"uint  getShiftXSubsequences(__global const ulong   *ks,\n" \
"__global const long    *kcCores,\n" \
"__global const uint4   *seqData,\n" \
"__global const uint4   *subseqData,\n" \
"__global const ushort   *prlIndices,\n" \
"__global const uint    *cssIndices,\n" \
"__global const uint    *cssSubseqs,\n" \
"ulong   *resBD,\n" \
"ulong4   mmPrime,\n" \
"ulong    pShift,\n" \
"uint     r,\n" \
"ulong   *resX,\n" \
"ulong4  *usableSubsequence);\n" \
"ulong climbLadder(__global const ushort  *ladders,\n" \
"ulong4   mmPrime,\n" \
"ulong    resBase,\n" \
"ulong   *resBD);\n" \
"ulong doBabySteps(ulong4 mmPrime,\n" \
"ulong resInvBase,\n" \
"uint bSteps,\n" \
"ushort *h_table,\n" \
"ushort *h_olist,\n" \
"ulong *h_BJ64);\n" \
"void  hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"uint  hashLookup(ulong bj, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"ulong mmmInvert(ulong p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong4 p);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong4 p);\n" \
"ulong mmmNToRes(ulong n, ulong4 p);\n" \
"ulong mmmPowmod(ulong resbase, ulong exp, ulong4 p);\n" \
"#if defined(USE_OPENCL)\n" \
"__kernel void cisonemultiple_kernel(__global const ulong   *primes,\n" \
"__global const ulong   *ks,\n" \
"__global const long    *kcCores,\n" \
"__global const uint4   *seqData,\n" \
"__global const uint4   *subseqData,\n" \
"__global const short   *divisorShifts,\n" \
"__global const ushort  *prlIndices,\n" \
"__global const uint    *cssIndices,\n" \
"__global const uint    *cssSubseqs,\n" \
"__global const ushort  *ladders,\n" \
"volatile __global       uint    *factorCount,\n" \
"__global       ulong4  *factors)\n" \
"#else\n" \
"kernel void cisonemultiple_kernel(device const ulong   *primes,\n" \
"device const ulong   *ks,\n" \
"device const long    *kcCores,\n" \
"device const uint4   *seqData,\n" \
"device const uint4   *subseqData,\n" \
"device const short   *divisorShifts,\n" \
"device const ushort  *prlIndices,\n" \
"device const uint    *cssIndices,\n" \
"device const uint    *cssSubseqs,\n" \
"device const ushort  *ladders,\n" \
"volatile device atomic_int   *factorCount,\n" \
"device       ulong4 *factors)\n" \
"uint    tpid [[thread_position_in_grid]])\n" \
"#endif\n" \
"{\n" \
"#if defined(USE_OPENCL)\n" \
"int    gid = get_global_id(0);\n" \
"#else\n" \
"int    gid = tpid;\n" \
"#endif\n" \
"ulong4 mmPrime;\n" \
"ulong  thePrime;\n" \
"MM_P(mmPrime)   = thePrime = primes[gid];\n" \
"MM_Q(mmPrime)   = mmmInvert(thePrime);\n" \
"MM_ONE(mmPrime) = mmmOne(thePrime);\n" \
"MM_R2(mmPrime)  = mmmR2(mmPrime);\n" \
"ulong   invBase = invmod(BASE, thePrime);\n" \
"ulong   resBase = mmmNToRes(BASE, mmPrime);\n" \
"ulong   resInvBase = mmmNToRes(invBase, mmPrime);\n" \
"ushort  h_table[HASH_SIZE];\n" \
"ushort  h_olist[HASH_ELEMENTS];\n" \
"ulong   h_BJ64[HASH_ELEMENTS+1];\n" \
"ulong4  usableSubsequences[SUBSEQUENCES];\n" \
"ulong   resBDCK[SUBSEQUENCES];\n" \
"ulong   resBD[BEST_Q];\n" \
"ulong   resX[POWER_RESIDUE_LCM];\n" \
"ulong   resBexpQ;\n" \
"ulong   ussIdx, ssCount;\n" \
"uint    babySteps, giantSteps;\n" \
"uint    idx, i, j;\n" \
"for (idx=0; idx<HASH_SIZE; idx+=2)\n" \
"{\n" \
"h_table[idx+0] = HASH_ELEMENTS;\n" \
"h_table[idx+1] = HASH_ELEMENTS;\n" \
"}\n" \
"h_BJ64[HASH_ELEMENTS] = ULONG_MAX;\n" \
"resBexpQ = climbLadder(ladders, mmPrime, resBase, resBD);\n" \
"ssCount = setupDiscreteLog(ks, kcCores, seqData, subseqData, divisorShifts, prlIndices, cssIndices, cssSubseqs, resBD, mmPrime, resInvBase, resX, usableSubsequences);\n" \
"if (ssCount == 0)\n" \
"return;\n" \
"babySteps  = SUBSEQ_BS(subseqData, ssCount-1);\n" \
"giantSteps = SUBSEQ_GS(subseqData, ssCount-1);\n" \
"ushort orderOfB = doBabySteps(mmPrime, resInvBase, babySteps, h_table, h_olist, h_BJ64);\n" \
"if (orderOfB > 0)\n" \
"{\n" \
"for (ussIdx=0; ussIdx<ssCount; ussIdx++)\n" \
"{\n" \
"j = hashLookup(USABLE_BDCK(usableSubsequences, ussIdx), h_table, h_olist, h_BJ64);\n" \
"while (j < babySteps * giantSteps)\n" \
"{\n" \
"uint q = USABLE_Q(usableSubsequences, ussIdx);\n" \
"collectFactor(thePrime, USABLE_SSIDX(usableSubsequences, ussIdx), N_TERM(q, 0, j), 0, factorCount, factors);\n" \
"j += orderOfB;\n" \
"}\n" \
"}\n" \
"}\n" \
"else\n" \
"{\n" \
"for (ussIdx=0; ussIdx<ssCount; ussIdx++)\n" \
"{\n" \
"j = hashLookup(USABLE_BDCK(usableSubsequences, ussIdx), h_table, h_olist, h_BJ64);\n" \
"if (j != HASH_NOT_FOUND)\n" \
"{\n" \
"uint q = USABLE_Q(usableSubsequences, ussIdx);\n" \
"collectFactor(thePrime, USABLE_SSIDX(usableSubsequences, ussIdx), N_TERM(q, 0, j), 0, factorCount, factors);\n" \
"}\n" \
"}\n" \
"if (giantSteps > 1)\n" \
"{\n" \
"ulong resBQM = mmmPowmod(resBexpQ, babySteps, mmPrime);\n" \
"for (i=1; i<giantSteps; i++)\n" \
"{\n" \
"for (ussIdx=0; ussIdx<ssCount; ussIdx++)\n" \
"{\n" \
"USABLE_BDCK(usableSubsequences, ussIdx) = mmmMulmod(USABLE_BDCK(usableSubsequences, ussIdx), resBQM, mmPrime);\n" \
"j = hashLookup(USABLE_BDCK(usableSubsequences, ussIdx), h_table, h_olist, h_BJ64);\n" \
"if (j != HASH_NOT_FOUND)\n" \
"{\n" \
"uint q = USABLE_Q(usableSubsequences, ussIdx);\n" \
"collectFactor(thePrime, USABLE_SSIDX(usableSubsequences, ussIdx), N_TERM(q, i, j), 0, factorCount, factors);\n" \
"}\n" \
"}\n" \
"}\n" \
"}\n" \
"}\n" \
"}\n" \
"ulong climbLadder(__global const ushort  *ladders,\n" \
"ulong4   mmPrime,\n" \
"ulong    resBase,\n" \
"ulong   *resBD)\n" \
"{\n" \
"uint  i, j, idx, lLen;\n" \
"lLen = *ladders;\n" \
"resBD[0] = MM_ONE(mmPrime);\n" \
"resBD[1] = resBase;\n" \
"resBD[2] = mmmMulmod(resBD[1], resBD[1], mmPrime);\n" \
"i = 2;\n" \
"for (j=0; j<lLen; j++)\n" \
"{\n" \
"idx = ladders[j+1];\n" \
"resBD[i+idx] = mmmMulmod(resBD[i], resBD[idx], mmPrime);\n" \
"i += idx;\n" \
"}\n" \
"return resBD[BEST_Q];\n" \
"}\n" \
"uint      setupDiscreteLog(__global const ulong   *ks,\n" \
"__global const long    *kcCores,\n" \
"__global const uint4   *seqData,\n" \
"__global const uint4   *subseqData,\n" \
"__global const short   *divisorShifts,\n" \
"__global const ushort  *prlIndices,\n" \
"__global const uint    *cssIndices,\n" \
"__global const uint    *cssSubseqs,\n" \
"ulong   *resBD,\n" \
"ulong4   mmPrime,\n" \
"ulong    resInvBase,\n" \
"ulong   *resX,\n" \
"ulong4  *usableSubsequence)\n" \
"{\n" \
"ulong      bm, pShift;\n" \
"uint       idx, r;\n" \
"int        shift;\n" \
"bm = MM_P(mmPrime) / 2;\n" \
"idx = bm % (POWER_RESIDUE_LCM/2);\n" \
"shift = divisorShifts[idx];\n" \
"if (shift == 0)\n" \
"return getShift0Subsequences(ks, kcCores, seqData, subseqData, resBD, mmPrime, usableSubsequence);\n" \
"if (shift > 0)\n" \
"{\n" \
"pShift = MM_P(mmPrime) / shift;\n" \
"}\n" \
"else\n" \
"{\n" \
"pShift = MM_P(mmPrime) >> (-shift);\n" \
"shift = 1 << (-shift);\n" \
"}\n" \
"resX[0] = MM_ONE(mmPrime);\n" \
"resX[1] = mmmPowmod(resInvBase, pShift, mmPrime);\n" \
"for (r=1; resX[r] != resX[0]; r++)\n" \
"resX[r+1] = mmmMulmod(resX[r], resX[1], mmPrime);\n" \
"if (shift % r != 0)\n" \
"return 0;\n" \
"return getShiftXSubsequences(ks, kcCores, seqData, subseqData, prlIndices, cssIndices, cssSubseqs, resBD, mmPrime, pShift, r, resX, usableSubsequence);\n" \
"}\n" \
"uint      getShift0Subsequences(__global const ulong   *ks,\n" \
"__global const long    *kcCores,\n" \
"__global const uint4   *seqData,\n" \
"__global const uint4   *subseqData,\n" \
"ulong   *resBD,\n" \
"ulong4   mmPrime,\n" \
"ulong4  *usableSubsequence)\n" \
"{\n" \
"short      kcLegendre;\n" \
"short      bLegendre = 0;\n" \
"short      usable = 0;\n" \
"uint       j = 0;\n" \
"bLegendre = legendre(BASE, MM_P(mmPrime));\n" \
"for (uint seqIdx=0; seqIdx<SEQUENCES; seqIdx++)\n" \
"{\n" \
"kcLegendre = legendre(kcCores[seqIdx], MM_P(mmPrime));\n" \
"switch (SEQ_PARITY(seqData, seqIdx))\n" \
"{\n" \
"case SP_EVEN:\n" \
"usable = (kcLegendre == 1 ? 1 : 0);\n" \
"break;\n" \
"case SP_ODD:\n" \
"usable = (kcLegendre == bLegendre ? 1 : 0);\n" \
"break;\n" \
"case SP_MIXED:\n" \
"usable = ((kcLegendre == 1 || kcLegendre == bLegendre) ? 1 : 0);\n" \
"break;\n" \
"}\n" \
"if (usable == 1)\n" \
"{\n" \
"ulong negCK = getNegCK(ks[seqIdx], (kcCores[seqIdx] < 0 ? 1 : -1), MM_P(mmPrime));\n" \
"ulong resNegCK = mmmNToRes(negCK, mmPrime);\n" \
"uint minss = SEQ_MINSS(seqData, seqIdx);\n" \
"uint maxss = SEQ_MAXSS(seqData, seqIdx);\n" \
"for (uint ssIdx=minss; ssIdx<maxss; ssIdx++)\n" \
"{\n" \
"uint q = SUBSEQ_Q(subseqData, ssIdx);\n" \
"USABLE_SEQIDX(usableSubsequence, j) = seqIdx;\n" \
"USABLE_SSIDX(usableSubsequence, j) = ssIdx;\n" \
"USABLE_Q(usableSubsequence, j) = q;\n" \
"USABLE_BDCK(usableSubsequence, j) = mmmMulmod(resBD[q], resNegCK, mmPrime);\n" \
"j++;\n" \
"}\n" \
"}\n" \
"}\n" \
"return j;\n" \
"}\n" \
"uint  getShiftXSubsequences(__global const ulong   *ks,\n" \
"__global const long    *kcCores,\n" \
"__global const uint4   *seqData,\n" \
"__global const uint4   *subseqData,\n" \
"__global const ushort  *prlIndices,\n" \
"__global const uint    *cssIndices,\n" \
"__global const uint    *cssSubseqs,\n" \
"ulong   *resBD,\n" \
"ulong4   mmPrime,\n" \
"ulong    pShift,\n" \
"uint     r,\n" \
"ulong   *resX,\n" \
"ulong4  *usableSubsequence)\n" \
"{\n" \
"short      kcLegendre;\n" \
"short      bLegendre = 0;\n" \
"short      usable = 0;\n" \
"uint       qr_mod, rIdx;\n" \
"uint       h, j = 0;\n" \
"uint       ssIdx, cssIdx;\n" \
"rIdx = prlIndices[r];\n" \
"bLegendre = legendre(BASE, MM_P(mmPrime));\n" \
"for (uint seqIdx=0; seqIdx<SEQUENCES; seqIdx++)\n" \
"{\n" \
"kcLegendre = legendre(kcCores[seqIdx], MM_P(mmPrime));\n" \
"switch (SEQ_PARITY(seqData, seqIdx))\n" \
"{\n" \
"case SP_EVEN:\n" \
"usable = (kcLegendre == 1 ? 1 : 0);\n" \
"break;\n" \
"case SP_ODD:\n" \
"usable = (kcLegendre == bLegendre ? 1 : 0);\n" \
"break;\n" \
"case SP_MIXED:\n" \
"usable = ((kcLegendre == 1 || kcLegendre == bLegendre) ? 1 : 0);\n" \
"break;\n" \
"}\n" \
"if (usable == 1)\n" \
"{\n" \
"ulong negCK = getNegCK(ks[seqIdx], (kcCores[seqIdx] < 0 ? 1 : -1), MM_P(mmPrime));\n" \
"ulong resNegCK = mmmNToRes(negCK, mmPrime);\n" \
"ulong resPowNegCK = mmmPowmod(resNegCK, pShift, mmPrime);\n" \
"resX[r] = resPowNegCK;\n" \
"for (h=0; resX[r] != resX[h]; h++)\n" \
";\n" \
"if (h < r)\n" \
"{\n" \
"qr_mod = cssIdx = CSS_INDEX(SEQ_SEQIDX(seqData, seqIdx), rIdx, h);\n" \
"cssIdx = cssIndices[cssIdx];\n" \
"if (cssIdx > 0 && cssSubseqs[cssIdx] > 0)\n" \
"{\n" \
"uint count = cssSubseqs[cssIdx];\n" \
"for (uint idx=1; idx<=count; idx++)\n" \
"{\n" \
"uint ssIdx = cssSubseqs[cssIdx+idx];\n" \
"uint q = SUBSEQ_Q(subseqData, ssIdx);\n" \
"USABLE_SEQIDX(usableSubsequence, j) = seqIdx;\n" \
"USABLE_SSIDX(usableSubsequence, j) = ssIdx;\n" \
"USABLE_Q(usableSubsequence, j) = q;\n" \
"USABLE_BDCK(usableSubsequence, j) = mmmMulmod(resBD[q], resNegCK, mmPrime);\n" \
"j++;\n" \
"}\n" \
"}\n" \
"}\n" \
"}\n" \
"}\n" \
"return j;\n" \
"}\n" \
"ulong doBabySteps(ulong4  mmPrime,\n" \
"ulong   resInvBase,\n" \
"uint    bSteps,\n" \
"ushort *h_table,\n" \
"ushort *h_olist,\n" \
"ulong  *h_BJ64)\n" \
"{\n" \
"ulong    resInvBaseExpQ;\n" \
"ulong    resBJ;\n" \
"ulong    firstResBJ;\n" \
"uint     j;\n" \
"resInvBaseExpQ = mmmPowmod(resInvBase, BEST_Q, mmPrime);\n" \
"firstResBJ = resBJ = mmmPowmod(resInvBaseExpQ, SIEVE_LOW, mmPrime);\n" \
"for (j=0; j<bSteps; j++)\n" \
"{\n" \
"hashInsert(resBJ, j, h_table, h_olist, h_BJ64);\n" \
"resBJ = mmmMulmod(resBJ, resInvBaseExpQ, mmPrime);\n" \
"if (resBJ == firstResBJ)\n" \
"return j + 1;\n" \
"}\n" \
"return 0;\n" \
"}\n" \
"ulong mmmInvert(ulong p)\n" \
"{\n" \
"ulong p_inv = 1;\n" \
"ulong prev = 0;\n" \
"while (p_inv != prev)\n" \
"{\n" \
"prev = p_inv;\n" \
"p_inv *= (2 - p * p_inv);\n" \
"}\n" \
"return p_inv;\n" \
"}\n" \
"inline ulong mmmOne(ulong _p)\n" \
"{\n" \
"return ((-_p) % _p);\n" \
"}\n" \
"inline ulong mmmR2(ulong4 mmPrime)\n" \
"{\n" \
"ulong t = mmmAdd(MM_ONE(mmPrime), MM_ONE(mmPrime), MM_P(mmPrime));\n" \
"t = mmmAdd(t, t, MM_P(mmPrime));\n" \
"for (size_t i=0; i<5; i++)\n" \
"t = mmmMulmod(t, t, mmPrime);\n" \
"return t;\n" \
"}\n" \
"inline ulong mmmAdd(ulong a, ulong b, ulong _p)\n" \
"{\n" \
"ulong c = (a >= _p - b) ? _p : 0;\n" \
"return a + b - c;\n" \
"}\n" \
"inline ulong mmmSub(ulong a, ulong b, ulong _p)\n" \
"{\n" \
"ulong c = (a < b) ? _p : 0;\n" \
"return a - b + c;\n" \
"}\n" \
"inline ulong mmmNToRes(ulong n, ulong4 mmPrime)\n" \
"{\n" \
"return mmmMulmod(n, MM_R2(mmPrime), mmPrime);\n" \
"}\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong4 mmPrime)\n" \
"{\n" \
"ulong lo = a * b;\n" \
"ulong hi = MUL_HI(a, b);\n" \
"ulong m = lo * MM_Q(mmPrime);\n" \
"ulong hi2 = MUL_HI(m, MM_P(mmPrime));\n" \
"long r = (long) hi - (long) hi2;\n" \
"if (r < 0)\n" \
"return (ulong) (r + MM_P(mmPrime));\n" \
"return (ulong) r;\n" \
"}\n" \
"ulong   mmmPowmod(ulong resbase, ulong exp, ulong4 mmPrime)\n" \
"{\n" \
"ulong x = resbase;\n" \
"ulong y = MM_ONE(mmPrime);\n" \
"while (true)\n" \
"{\n" \
"if (exp & 1)\n" \
"y = mmmMulmod(x, y, mmPrime);\n" \
"exp >>= 1;\n" \
"if (!exp)\n" \
"return y;\n" \
"x = mmmMulmod(x, x, mmPrime);\n" \
"}\n" \
"return 0;\n" \
"}\n" \
"ulong  invmod(ulong a, ulong p)\n" \
"{\n" \
"ulong ps1, ps2, q, r, t, dividend, divisor;\n" \
"uint parity;\n" \
"if (a < 3)\n" \
"return (a < 2) ? a : (p+1)/2;\n" \
"q = p / a;\n" \
"r = p % a;\n" \
"dividend = a;\n" \
"divisor = r;\n" \
"ps1 = q;\n" \
"ps2 = 1;\n" \
"parity = 0;\n" \
"while (divisor > 1)\n" \
"{\n" \
"r = dividend - divisor;\n" \
"t = r - divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t;\n" \
"if (r >= divisor) {\n" \
"q = dividend / divisor;\n" \
"r = dividend % divisor;\n" \
"q *= ps1;\n" \
"} } } } } } } } }\n" \
"q += ps2;\n" \
"parity = ~parity;\n" \
"dividend = divisor;\n" \
"divisor = r;\n" \
"ps2 = ps1;\n" \
"ps1 = q;\n" \
"}\n" \
"return (parity) ? ps1 : p - ps1;\n" \
"}\n" \
"short  legendre(long a, ulong p)\n" \
"{\n" \
"ulong x, y, t;\n" \
"short sign;\n" \
"if (a < 0)\n" \
"{\n" \
"a = -a;\n" \
"sign = (p % 4 == 1) ? 1 : -1;\n" \
"}\n" \
"else\n" \
"sign = 1;\n" \
"for (y = a; y % 2 == 0; y /= 2)\n" \
"if (p % 8 == 3 || p % 8 == 5)\n" \
"sign = -sign;\n" \
"if (p % 4 == 3 && y % 4 == 3)\n" \
"sign = -sign;\n" \
"for (x = p % y; x > 0; x %= y)\n" \
"{\n" \
"for ( ; x % 2 == 0; x /= 2)\n" \
"if (y % 8 == 3 || y % 8 == 5)\n" \
"sign = -sign;\n" \
"t = x, x = y, y = t;\n" \
"if (x % 4 == 3 && y % 4 == 3)\n" \
"sign = -sign;\n" \
"}\n" \
"return sign;\n" \
"}\n" \
"ulong getNegCK(ulong k, int c, ulong p)\n" \
"{\n" \
"ulong negCK;\n" \
"if (p < k)\n" \
"negCK = k % p;\n" \
"else\n" \
"negCK = k;\n" \
"if (c > 0)\n" \
"negCK = p - negCK;\n" \
"return negCK;\n" \
"}\n" \
"void hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64)\n" \
"{\n" \
"uint slot;\n" \
"h_BJ64[j] = bj;\n" \
"slot = bj & HASH_SIZE_M1;\n" \
"if (h_table[slot] == HASH_ELEMENTS)\n" \
"h_table[slot] = j;\n" \
"else\n" \
"{\n" \
"h_olist[j] = (h_table[slot] ^ HASH_MASK1);\n" \
"h_table[slot] = (j | HASH_MASK1);\n" \
"}\n" \
"}\n" \
"uint hashLookup(ulong bj, ushort *h_table, ushort *h_olist, ulong *h_BJ64)\n" \
"{\n" \
"uint   slot;\n" \
"ushort elt;\n" \
"slot = bj & HASH_SIZE_M1;\n" \
"elt = h_table[slot];\n" \
"if (h_BJ64[elt & HASH_MASK2] == bj)\n" \
"return elt & HASH_MASK2;\n" \
"if ((elt & HASH_MASK1) == 0)\n" \
"return HASH_NOT_FOUND;\n" \
"elt &= HASH_MASK2;\n" \
"do\n" \
"{\n" \
"if (h_BJ64[h_olist[elt] & HASH_MASK2] == bj)\n" \
"return h_olist[elt] & HASH_MASK2;\n" \
"elt = h_olist[elt];\n" \
"} while ((elt & HASH_MASK1) == 0);\n" \
"return HASH_NOT_FOUND;\n" \
"}\n" \
"#if defined(USE_OPENCL)\n" \
"void collectFactor(ulong   p,\n" \
"ulong   ssIdx,\n" \
"ulong   n,\n" \
"ulong   x,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors)\n" \
"#else\n" \
"void collectFactor(ulong       p,\n" \
"uint        ssIdx,\n" \
"uint        n,\n" \
"volatile device atomic_int *factorCount,\n" \
"device ulong4     *factors)\n" \
"#endif\n" \
"{\n" \
"#if defined(USE_OPENCL)\n" \
"int old = atomic_inc(factorCount);\n" \
"#else\n" \
"int old = atomic_fetch_add_explicit(factorCount, 1, memory_order_relaxed);\n" \
"#endif\n" \
"if (old >= MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = p;\n" \
"factors[old].y = ssIdx;\n" \
"factors[old].z = n;\n" \
"factors[old].w = x;\n" \
"}\n" \
;
#endif
#endif
