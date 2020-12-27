#ifdef __cplusplus
#ifndef _CW_KERNEL_CL
#define _CW_KERNEL_CL
const char *cw_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n" \
"void collect_factor(int  n,\n" \
"int  c,\n" \
"long p,\n" \
"volatile __global uint *factorsCount,\n" \
"__global long4 *factors);\n" \
"long compute_inverse(long a, long modulus);\n" \
"ulong mmmInvert(ulong p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmN(ulong n, ulong _p);\n" \
"ulong mmmPowmod(ulong base, ulong exp, ulong _p, ulong _q);\n" \
"#define MAX_POWERS 50\n" \
"__kernel void cw_kernel(__global const  long  *primes,\n" \
"__global const  int   *terms,\n" \
"volatile __global        uint  *factorsCount,\n" \
"__global        long4 *factors)\n" \
"{\n" \
"int gid = get_global_id(0);\n" \
"int   termIndex;\n" \
"int   theN, prevN;\n" \
"int   power, idx;\n" \
"long  thePrime = primes[gid];\n" \
"long  rem = 0;\n" \
"#ifdef CHECK_CULLEN\n" \
"ulong cullen;\n" \
"#endif\n" \
"#ifdef CHECK_WOODALL\n" \
"ulong woodall;\n" \
"#endif\n" \
"ulong inverse;\n" \
"ulong powers[MAX_POWERS+1];\n" \
"ulong resC[MAX_POWERS+1];\n" \
"ulong resW[MAX_POWERS+1];\n" \
"inverse = compute_inverse(BASE, thePrime);\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"theN = terms[0];\n" \
"rem = mmmPowmod(inverse, theN, thePrime, _q);\n" \
"#ifdef CHECK_CULLEN\n" \
"cullen  = mmmN(thePrime - theN, thePrime);\n" \
"if (rem == cullen)\n" \
"collect_factor(theN, +1, thePrime, factorsCount, factors);\n" \
"#endif\n" \
"#ifdef CHECK_WOODALL\n" \
"woodall = mmmN(theN, thePrime);\n" \
"if (rem == woodall)\n" \
"collect_factor(theN, -1, thePrime, factorsCount, factors);\n" \
"#endif\n" \
"powers[1] = mmmN(BASE, thePrime);\n" \
"powers[2] = mmmMulmod(powers[1], powers[1], thePrime, _q);\n" \
"resC[1] = mmmN(thePrime - 1, thePrime);\n" \
"resW[1] = mmmN(1, thePrime);\n" \
"resC[2] = mmmAdd(resC[1], resC[1], thePrime);\n" \
"resW[2] = mmmAdd(resW[1], resW[1], thePrime);\n" \
"if (BASE & 1)\n" \
"{\n" \
"for (idx=4; idx<MAX_POWERS+1; idx+=2)\n" \
"{\n" \
"powers[idx] = mmmMulmod(powers[idx-2], powers[2], thePrime, _q);\n" \
"resC[idx] = mmmAdd(resC[idx-2], resC[2], thePrime);\n" \
"resW[idx] = mmmAdd(resW[idx-2], resW[2], thePrime);\n" \
"}\n" \
"}\n" \
"else\n" \
"{\n" \
"powers[idx] = mmmMulmod(powers[idx-1], powers[1], thePrime, _q);\n" \
"resC[idx] = mmmAdd(resC[idx-1], resC[1], thePrime);\n" \
"resW[idx] = mmmAdd(resW[idx-1], resW[1], thePrime);\n" \
"}\n" \
"prevN = terms[0];\n" \
"termIndex = 1;\n" \
"while (terms[termIndex] > 0)\n" \
"{\n" \
"theN = terms[termIndex];\n" \
"power = prevN - theN;\n" \
"while (power > MAX_POWERS)\n" \
"{\n" \
"rem = mmmMulmod(rem, powers[MAX_POWERS], thePrime, _q);\n" \
"cullen = mmmSub(cullen, resC[MAX_POWERS], thePrime);\n" \
"woodall = mmmSub(woodall, resW[MAX_POWERS], thePrime);\n" \
"power -= MAX_POWERS;\n" \
"}\n" \
"rem = mmmMulmod(rem, powers[power], thePrime, _q);\n" \
"cullen = mmmSub(cullen, resC[power], thePrime);\n" \
"woodall = mmmSub(woodall, resW[power], thePrime);\n" \
"#ifdef CHECK_CULLEN\n" \
"if (rem == cullen)\n" \
"collect_factor(theN, +1, thePrime, factorsCount, factors);\n" \
"#endif\n" \
"#ifdef CHECK_WOODALL\n" \
"if (rem == woodall)\n" \
"collect_factor(theN, -1, thePrime, factorsCount, factors);\n" \
"#endif\n" \
"prevN = theN;\n" \
"termIndex++;\n" \
"};\n" \
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
"ulong   mmmPowmod(ulong base, ulong exp, ulong _p, ulong _q)\n" \
"{\n" \
"ulong x = mmmN(base, _p);\n" \
"ulong y = mmmN(1, _p);\n" \
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
"void collect_factor(int  n,\n" \
"int  c,\n" \
"long p,\n" \
"volatile __global uint *factorsCount,\n" \
"__global long4 *factors)\n" \
"{\n" \
"int old = atomic_inc(factorsCount);\n" \
"if (old >= MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = n;\n" \
"factors[old].y = c;\n" \
"factors[old].z = p;\n" \
"}\n" \
"long compute_inverse(long a, long p)\n" \
"{\n" \
"long     m = p;\n" \
"long     y = 0, x = 1;\n" \
"long     q, t;\n" \
"while (a > 1)\n" \
"{\n" \
"q = a / m;\n" \
"t = m;\n" \
"m = a % m;\n" \
"a = t;\n" \
"t = y;\n" \
"y = x - q * y;\n" \
"x = t;\n" \
"}\n" \
"if (x < 0)\n" \
"x += p;\n" \
"return x;\n" \
"}\n" \
;
#endif
#endif
