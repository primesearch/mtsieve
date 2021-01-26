#ifdef __cplusplus
#ifndef _SR_KERNEL_CL
#define _SR_KERNEL_CL
const char *sr_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n" \
"#define HASH_NOT_FOUND     UINT_MAX\n" \
"#define HASH_MASK1         (1<<15)\n" \
"#define HASH_MASK2         (HASH_MASK1-1)\n" \
"ulong  invmod(ulong a, ulong p);\n" \
"inline ulong  smod(long a, ulong p)\n" \
"{\n" \
"ulong ua;\n" \
"if (a >= 0)\n" \
"{\n" \
"ua = (ulong) a;\n" \
"return (ua < p) ? ua : ua%p;\n" \
"}\n" \
"else\n" \
"{\n" \
"ua = (ulong) -a;\n" \
"return (ua < p) ? p-ua : p-(ua%p);\n" \
"}\n" \
"}\n" \
"inline ulong  lmod(long a, ulong p)\n" \
"{\n" \
"if (a >= 0)\n" \
"return (a < p) ? a : a%p;\n" \
"else\n" \
"return (-a < p) ? p-(-a) : p-(-a%p);\n" \
"}\n" \
"inline ulong  umod(ulong a, ulong p)\n" \
"{\n" \
"return (a < p) ? a : a%p;\n" \
"}\n" \
"ulong buildTables(ulong thePrime, ulong _q, ulong _r2, ulong _one, ulong *resBDCK64,\n" \
"__global const ulong  *SEQ_K,\n" \
"__global const  long  *SEQ_C,\n" \
"__global const uint   *SUBSEQ_SEQ,\n" \
"__global const uint   *SUBSEQ_Q);\n" \
"ulong babySteps(ulong thePrime, ulong _q, ulong _r2, ulong _one, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"void  hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"uint  hashLookup(ulong bj, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"void collectFactor(ulong   p,\n" \
"uint    ssIdx,\n" \
"uint    n,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors);\n" \
"ulong mmmInvert(ulong p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one);\n" \
"__kernel void sr_kernel(__global const ulong  *primes,\n" \
"__global const ulong  *SEQ_K,\n" \
"__global const  long  *SEQ_C,\n" \
"__global const uint   *SUBSEQ_SEQ,\n" \
"__global const uint   *SUBSEQ_Q,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors)\n" \
"{\n" \
"int gid = get_global_id(0);\n" \
"ushort  h_table[HASH_SIZE];\n" \
"ushort  h_olist[HASH_ELEMENTS];\n" \
"ulong   h_BJ64[HASH_ELEMENTS+1];\n" \
"ulong   resBDCK64[SUBSEQUENCES];\n" \
"ulong   thePrime = primes[gid];\n" \
"uint    idx, j, ssIdx;\n" \
"for (idx=0; idx<HASH_SIZE; idx++)\n" \
"h_table[idx] = 0;\n" \
"for (idx=0; idx<HASH_ELEMENTS; idx++)\n" \
"h_BJ64[idx] = 0;\n" \
"h_BJ64[HASH_ELEMENTS] = ULONG_MAX;\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong _one = mmmOne(thePrime);\n" \
"ulong _r2 = mmmR2(thePrime, _q, _one);\n" \
"ulong resBM64 = buildTables(thePrime, _q, _r2, _one, resBDCK64, SEQ_K, SEQ_C, SUBSEQ_SEQ, SUBSEQ_Q);\n" \
"uint orderOfB = babySteps(thePrime, _q, _r2, _one, h_table, h_olist, h_BJ64);\n" \
"if (orderOfB > 0)\n" \
"{\n" \
"for (ssIdx=0; ssIdx<SUBSEQUENCES; ssIdx++)\n" \
"{\n" \
"j = hashLookup(resBDCK64[ssIdx], h_table, h_olist, h_BJ64);\n" \
"while (j < SIEVE_RANGE)\n" \
"{\n" \
"collectFactor(thePrime, ssIdx, N_TERM(ssIdx, 0, j), factorCount, factors);\n" \
"j += orderOfB;\n" \
"}\n" \
"}\n" \
"}\n" \
"else\n" \
"{\n" \
"for (ssIdx=0; ssIdx<SUBSEQUENCES; ssIdx++)\n" \
"{\n" \
"j = hashLookup(resBDCK64[ssIdx], h_table, h_olist, h_BJ64);\n" \
"if (j != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, ssIdx, N_TERM(ssIdx, 0, j), factorCount, factors);\n" \
"}\n" \
"}\n" \
"if (GIANT_STEPS < 2)\n" \
"return;\n" \
"resBM64 = mmmPowmod(resBM64, BABY_STEPS, thePrime, _q, _one);\n" \
"for (uint i=1; i<GIANT_STEPS; i++)\n" \
"{\n" \
"for (ssIdx=0; ssIdx<SUBSEQUENCES; ssIdx++)\n" \
"{\n" \
"resBDCK64[ssIdx] = mmmMulmod(resBDCK64[ssIdx], resBM64, thePrime, _q);\n" \
"uint j = hashLookup(resBDCK64[ssIdx], h_table, h_olist, h_BJ64);\n" \
"if (j != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, ssIdx, N_TERM(ssIdx, i, j), factorCount, factors);\n" \
"}\n" \
"}\n" \
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
"inline ulong mmmR2(ulong _p, ulong _q, ulong _one)\n" \
"{\n" \
"ulong t = mmmAdd(_one, _one, _p);\n" \
"t = mmmAdd(t, t, _p);\n" \
"for (size_t i=0; i<5; i++)\n" \
"t = mmmMulmod(t, t, _p, _q);\n" \
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
"inline ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2)\n" \
"{\n" \
"return mmmMulmod(n, _r2, _p, _q);\n" \
"}\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q)\n" \
"{\n" \
"ulong lo = a * b;\n" \
"ulong hi = mul_hi(a, b);\n" \
"ulong m = lo * _q;\n" \
"ulong hi2 = mul_hi(m, _p);\n" \
"long r = (long) hi - (long) hi2;\n" \
"if (r < 0)\n" \
"return (ulong) (r + _p);\n" \
"return (ulong) r;\n" \
"}\n" \
"ulong   mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one)\n" \
"{\n" \
"ulong x = resbase;\n" \
"ulong y = _one;\n" \
"while (true)\n" \
"{\n" \
"if (exp & 1)\n" \
"y = mmmMulmod(x, y, _p, _q);\n" \
"exp >>= 1;\n" \
"if (!exp)\n" \
"return y;\n" \
"x = mmmMulmod(x, x, _p, _q);\n" \
"}\n" \
"return 0;\n" \
"}\n" \
"ulong babySteps(ulong thePrime, ulong _q, ulong _r2, ulong _one, ushort *h_table, ushort *h_olist, ulong *h_BJ64)\n" \
"{\n" \
"uint    j;\n" \
"ulong   resBase = mmmN(BASE, thePrime, _q, _r2);\n" \
"ulong   resBexpQ = mmmPowmod(resBase, BESTQ, thePrime, _q, _one);\n" \
"ulong   firstResBJ = mmmPowmod(resBexpQ, SIEVE_LOW, thePrime, _q, _one);\n" \
"ulong   resBJ = firstResBJ;\n" \
"for (j=0; j<BABY_STEPS; j++)\n" \
"{\n" \
"hashInsert(resBJ, j, h_table, h_olist, h_BJ64);\n" \
"resBJ = mmmMulmod(resBJ, resBexpQ, thePrime, _q);\n" \
"if (resBJ == firstResBJ)\n" \
"return j + 1;\n" \
"}\n" \
"return 0;\n" \
"}\n" \
"ulong buildTables(ulong thePrime, ulong _q, ulong _r2, ulong _one, ulong *resBDCK64,\n" \
"__global const ulong  *SEQ_K,\n" \
"__global const  long  *SEQ_C,\n" \
"__global const uint   *SUBSEQ_SEQ,\n" \
"__global const uint   *SUBSEQ_Q)\n" \
"{\n" \
"uint  qIdx, seqIdx, ssIdx;\n" \
"ulong resBD64[BESTQ];\n" \
"ulong resCK64[SEQUENCES];\n" \
"ulong firstResBM64 = mmmN(invmod(BASE, thePrime), thePrime, _q, _r2);\n" \
"ulong resBM64 = firstResBM64;\n" \
"resBD64[0] = _one;\n" \
"for (qIdx=1; qIdx<BESTQ; qIdx++)\n" \
"{\n" \
"resBD64[qIdx] = resBM64;\n" \
"resBM64 = mmmMulmod(resBM64, firstResBM64, thePrime, _q);\n" \
"}\n" \
"for (seqIdx=0; seqIdx<SEQUENCES; seqIdx++)\n" \
"{\n" \
"ulong v1 = umod(SEQ_K[seqIdx], thePrime);\n" \
"ulong v2 = lmod(-SEQ_C[seqIdx], thePrime);\n" \
"ulong v3 = invmod(v1, thePrime);\n" \
"v2 = mmmN(v2, thePrime, _q, _r2);\n" \
"v3 = mmmN(v3, thePrime, _q, _r2);\n" \
"resCK64[seqIdx] = mmmMulmod(v2, v3, thePrime, _q);\n" \
"}\n" \
"for (ssIdx=0; ssIdx<SUBSEQUENCES; ssIdx++)\n" \
"{\n" \
"ulong ck = resCK64[SUBSEQ_SEQ[ssIdx]];\n" \
"ulong bd = resBD64[SUBSEQ_Q[ssIdx]];\n" \
"resBDCK64[ssIdx] = mmmMulmod(bd, ck, thePrime, _q);\n" \
"}\n" \
"return resBM64;\n" \
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
"void hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64)\n" \
"{\n" \
"uint slot;\n" \
"h_BJ64[j] = bj;\n" \
"slot = bj & (HASH_SIZE - 1);\n" \
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
"slot = bj & (HASH_SIZE - 1);\n" \
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
"void collectFactor(ulong   p,\n" \
"uint    ssIdx,\n" \
"uint    n,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors)\n" \
"{\n" \
"int old = atomic_inc(factorCount);\n" \
"if (old >= MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = ssIdx;\n" \
"factors[old].y = n;\n" \
"factors[old].z = p;\n" \
"}\n" \
;
#endif
#endif
