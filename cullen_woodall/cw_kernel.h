#ifdef __cplusplus
#ifndef _CW_KERNEL_CL
#define _CW_KERNEL_CL
const char *cw_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n" \
"void collect_factor(int  n,\n" \
"int  c,\n" \
"long p,\n" \
"volatile __global uint *factorCount,\n" \
"__global long4 *factors);\n" \
"long compute_inverse(long a, long p);\n" \
"ulong mmmInvert(ulong _p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2);\n" \
"#define MAX_POWERS 50\n" \
"__kernel void cw_kernel(__global const  long  *primes,\n" \
"__global const  int   *terms,\n" \
"volatile __global        uint  *factorCount,\n" \
"__global        long4 *factors)\n" \
"{\n" \
"int gid = get_global_id(0);\n" \
"ulong thePrime = primes[gid];\n" \
"ulong powers[MAX_POWERS+1];\n" \
"ulong resC[MAX_POWERS+1];\n" \
"ulong resW[MAX_POWERS+1];\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong _one = mmmOne(thePrime);\n" \
"ulong _r2 = mmmR2(thePrime, _q, _one);\n" \
"ulong inverse = compute_inverse(BASE, thePrime);\n" \
"int theN = terms[0];\n" \
"ulong resI = mmmN(inverse, thePrime, _q, _r2);\n" \
"ulong rem = mmmPowmod(resI, theN, thePrime, _q, _one, _r2);\n" \
"#ifdef CHECK_CULLEN\n" \
"ulong cullen  = mmmN(thePrime - theN, thePrime, _q, _r2);\n" \
"if (rem == cullen)\n" \
"collect_factor(theN, +1, thePrime, factorCount, factors);\n" \
"#endif\n" \
"#ifdef CHECK_WOODALL\n" \
"ulong woodall = mmmN(theN, thePrime, _q, _r2);\n" \
"if (rem == woodall)\n" \
"collect_factor(theN, -1, thePrime, factorCount, factors);\n" \
"#endif\n" \
"powers[1] = mmmN(BASE, thePrime, _q, _r2);\n" \
"powers[2] = mmmMulmod(powers[1], powers[1], thePrime, _q);\n" \
"resC[1] = mmmN(thePrime - 1, thePrime, _q, _r2);\n" \
"resW[1] = _one;\n" \
"resC[2] = mmmAdd(resC[1], resC[1], thePrime);\n" \
"resW[2] = mmmAdd(resW[1], resW[1], thePrime);\n" \
"if (BASE & 1)\n" \
"{\n" \
"for (size_t idx=4; idx<MAX_POWERS+1; idx+=2)\n" \
"{\n" \
"powers[idx] = mmmMulmod(powers[idx-2], powers[2], thePrime, _q);\n" \
"resC[idx] = mmmAdd(resC[idx-2], resC[2], thePrime);\n" \
"resW[idx] = mmmAdd(resW[idx-2], resW[2], thePrime);\n" \
"}\n" \
"}\n" \
"else\n" \
"{\n" \
"for (size_t idx=3; idx<MAX_POWERS+1; idx++)\n" \
"{\n" \
"powers[idx] = mmmMulmod(powers[idx-1], powers[1], thePrime, _q);\n" \
"resC[idx] = mmmAdd(resC[idx-1], resC[1], thePrime);\n" \
"resW[idx] = mmmAdd(resW[idx-1], resW[1], thePrime);\n" \
"}\n" \
"}\n" \
"int prevN = terms[0];\n" \
"size_t termIndex = 1;\n" \
"while (terms[termIndex] > 0)\n" \
"{\n" \
"theN = terms[termIndex];\n" \
"size_t power = (size_t)(prevN - theN);\n" \
"while (power > MAX_POWERS)\n" \
"{\n" \
"rem = mmmMulmod(rem, powers[MAX_POWERS], thePrime, _q);\n" \
"#ifdef CHECK_CULLEN\n" \
"cullen = mmmSub(cullen, resC[MAX_POWERS], thePrime);\n" \
"#endif\n" \
"#ifdef CHECK_WOODALL\n" \
"woodall = mmmSub(woodall, resW[MAX_POWERS], thePrime);\n" \
"#endif\n" \
"power -= MAX_POWERS;\n" \
"}\n" \
"rem = mmmMulmod(rem, powers[power], thePrime, _q);\n" \
"#ifdef CHECK_CULLEN\n" \
"cullen = mmmSub(cullen, resC[power], thePrime);\n" \
"#endif\n" \
"#ifdef CHECK_WOODALL\n" \
"woodall = mmmSub(woodall, resW[power], thePrime);\n" \
"#endif\n" \
"#ifdef CHECK_CULLEN\n" \
"if (rem == cullen)\n" \
"collect_factor(theN, +1, thePrime, factorCount, factors);\n" \
"#endif\n" \
"#ifdef CHECK_WOODALL\n" \
"if (rem == woodall)\n" \
"collect_factor(theN, -1, thePrime, factorCount, factors);\n" \
"#endif\n" \
"prevN = theN;\n" \
"termIndex++;\n" \
"};\n" \
"}\n" \
"ulong mmmInvert(ulong _p)\n" \
"{\n" \
"ulong p_inv = 1;\n" \
"ulong prev = 0;\n" \
"while (p_inv != prev)\n" \
"{\n" \
"prev = p_inv;\n" \
"p_inv *= (2 - _p * p_inv);\n" \
"}\n" \
"return p_inv;\n" \
"}\n" \
"ulong mmmOne(ulong _p)\n" \
"{\n" \
"return ((-_p) % _p);\n" \
"}\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one)\n" \
"{\n" \
"ulong t = mmmAdd(_one, _one, _p);\n" \
"t = mmmAdd(t, t, _p);\n" \
"for (size_t i=0; i<5; i++)\n" \
"t = mmmMulmod(t, t, _p, _q);\n" \
"return t;\n" \
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
"ulong lo, hi;\n" \
"#ifdef __NV_CL_C_VERSION\n" \
"const uint a0 = (uint)(a), a1 = (uint)(a >> 32);\n" \
"const uint b0 = (uint)(b), b1 = (uint)(b >> 32);\n" \
"uint c0 = a0 * b0, c1 = mul_hi(a0, b0), c2, c3;\n" \
"asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a0), \"r\" (b1), \"r\" (c1));\n" \
"asm volatile (\"madc.hi.u32 %0, %1, %2, 0;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b1));\n" \
"asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a1), \"r\" (b1), \"r\" (c2));\n" \
"asm volatile (\"madc.hi.u32 %0, %1, %2, 0;\" : \"=r\" (c3) : \"r\" (a1), \"r\" (b1));\n" \
"asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a1), \"r\" (b0), \"r\" (c1));\n" \
"asm volatile (\"madc.hi.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a1), \"r\" (b0), \"r\" (c2));\n" \
"asm volatile (\"addc.u32 %0, %1, 0;\" : \"=r\" (c3) : \"r\" (c3));\n" \
"lo = upsample(c1, c0); hi = upsample(c3, c2);\n" \
"#else\n" \
"lo = a * b; hi = mul_hi(a, b);\n" \
"#endif\n" \
"ulong m = lo * _q;\n" \
"ulong mp = mul_hi(m, _p);\n" \
"long r = (long)(hi - mp);\n" \
"return (r < 0) ? r + _p : r;\n" \
"}\n" \
"ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2)\n" \
"{\n" \
"return mmmMulmod(n, _r2, _p, _q);\n" \
"}\n" \
"ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2)\n" \
"{\n" \
"ulong x = resB;\n" \
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
"void collect_factor(int  n,\n" \
"int  c,\n" \
"long p,\n" \
"volatile __global uint *factorCount,\n" \
"__global long4 *factors)\n" \
"{\n" \
"int old = atomic_inc(factorCount);\n" \
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
