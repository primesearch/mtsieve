#ifdef __cplusplus
#ifndef _MF_KERNEL_CL
#define _MF_KERNEL_CL
const char *mf_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"void collect_factor(int  n,\n" \
"int  c,\n" \
"long p,\n" \
"volatile __global int *factorsCount,\n" \
"__global long4 *factors);\n" \
"ulong mmmInvert(ulong p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmN(ulong n, ulong _p);\n" \
"__kernel void mf_kernel(__global const  ulong  *primes,\n" \
"__global        ulong2 *rems,\n" \
"__global        int    *params,\n" \
"volatile __global        int    *factorsCount,\n" \
"__global        long4  *factors)\n" \
"{\n" \
"int    gid = get_global_id(0);\n" \
"int    steps;\n" \
"int    idx = 0;\n" \
"int    maxNFirstLoop = D_MIN_N - D_MULTIFACTORIAL;\n" \
"ulong  startN = params[0];\n" \
"ulong  currentN = params[1];\n" \
"ulong  n;\n" \
"ulong  thePrime = primes[gid];\n" \
"if (thePrime == 0) return;\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong pOne = mmmOne(thePrime);\n" \
"ulong mOne = mmmSub(0, pOne, thePrime);\n" \
"ulong mfrs = mmmN(D_MULTIFACTORIAL, thePrime);;\n" \
"ulong ri, rf;\n" \
"if (startN == currentN)\n" \
"{\n" \
"if (startN > thePrime)\n" \
"ri = mmmN(startN%thePrime, thePrime);\n" \
"else\n" \
"ri = mmmN(startN, thePrime);\n" \
"rf = ri;\n" \
"}\n" \
"else\n" \
"{\n" \
"ri = rems[gid].x;\n" \
"rf = rems[gid].y;\n" \
"}\n" \
"steps = 0;\n" \
"n = currentN + D_MULTIFACTORIAL;\n" \
"for (; n<maxNFirstLoop; n+=D_MULTIFACTORIAL, steps++)\n" \
"{\n" \
"if (steps >= D_MAX_STEPS)\n" \
"break;\n" \
"ri = mmmAdd(ri, mfrs, thePrime);\n" \
"rf = mmmMulmod(rf, ri, thePrime, _q);\n" \
"}\n" \
"for (; n <=D_MAX_N; n+=D_MULTIFACTORIAL, steps++)\n" \
"{\n" \
"if (steps >= D_MAX_STEPS)\n" \
"break;\n" \
"ri = mmmAdd(ri, mfrs, thePrime);\n" \
"rf = mmmMulmod(rf, ri, thePrime, _q);\n" \
"if (rf == pOne)\n" \
"collect_factor(n, -1, thePrime, factorsCount, factors);\n" \
"if (rf == mOne)\n" \
"collect_factor(n, +1, thePrime, factorsCount, factors);\n" \
"}\n" \
"rems[gid].x = ri;\n" \
"rems[gid].y = rf;\n" \
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
"ulong mmmOne(ulong _p)\n" \
"{\n" \
"return ((-_p) % _p);\n" \
"}\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p)\n" \
"{\n" \
"ulong c = (a >= _p - b) ? _p : 0;\n" \
"return a + b - c;\n" \
"}\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p)\n" \
"{\n" \
"ulong c = (a < b) ? _p : 0;\n" \
"return a - b + c;\n" \
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
"ulong   mmmN(ulong n, ulong _p)\n" \
"{\n" \
"if (n == 1)\n" \
"return mmmOne(_p);\n" \
"ulong list[64];\n" \
"ulong bit = 0x01;\n" \
"ulong value = 0;\n" \
"list[0] = mmmOne(_p);\n" \
"if (n & bit)\n" \
"{\n" \
"value = list[0];\n" \
"n &= ~bit;\n" \
"}\n" \
"for (uint idx=1; idx<64; idx++)\n" \
"{\n" \
"bit <<= 1;\n" \
"list[idx] = mmmAdd(list[idx-1], list[idx-1], _p);\n" \
"if (n & bit)\n" \
"{\n" \
"value = mmmAdd(list[idx], value, _p);\n" \
"n &= ~bit;\n" \
"}\n" \
"if (n == 0)\n" \
"break;\n" \
"}\n" \
"return value;\n" \
"}\n" \
"void collect_factor(int  n,\n" \
"int  c,\n" \
"long p,\n" \
"volatile __global int *factorsCount,\n" \
"__global long4 *factors)\n" \
"{\n" \
"if (n < D_MIN_N)\n" \
"return;\n" \
"if (n > D_MAX_N)\n" \
"return;\n" \
"int old = atomic_inc(factorsCount);\n" \
"if (old >= D_MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = n;\n" \
"factors[old].y = c;\n" \
"factors[old].z = p;\n" \
"}\n" \
;
#endif
#endif
