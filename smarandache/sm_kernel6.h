#ifdef __cplusplus
#ifndef _SM_KERNEL6_CL
#define _SM_KERNEL6_CL
const char *sm_kernel6= \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"ulong  invmod(ulong a, ulong p);\n" \
"void collect_factor(uint    n,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong2 *factors);\n" \
"ulong mmmInvert(ulong _p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q);\n" \
"__kernel void sm_kernel6(__global const  ulong  *primes,\n" \
"__global const  uint   *terms,\n" \
"volatile __global        uint   *factorCount,\n" \
"__global        ulong2 *factors)\n" \
"{\n" \
"int    gid = get_global_id(0);\n" \
"ulong  thePrime = primes[gid];\n" \
"if (thePrime == 0) return;\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong _one = mmmOne(thePrime);\n" \
"ulong _r2 = mmmR2(thePrime, _q, _one);\n" \
"ulong invmod2 = invmod(12321, thePrime);\n" \
"ulong invmod3 = invmod(((ulong) 1234321) % thePrime, thePrime);\n" \
"ulong invmod4 = invmod(((ulong) 123454321) % thePrime, thePrime);\n" \
"ulong resTenE1 = mmmNtoRes(10, thePrime, _q, _r2);\n" \
"ulong tempMul1 = mmmNtoRes(((ulong) 15208068915062105958) % thePrime, thePrime, _q, _r2);\n" \
"ulong tempSub1 = mmmNtoRes(((ulong) 11211123422) % thePrime, thePrime, _q, _r2);\n" \
"ulong tempSub2 = mmmNtoRes(((ulong) 1109890222) % thePrime, thePrime, _q, _r2);\n" \
"ulong resInvmod2 = mmmNtoRes(invmod2, thePrime, _q, _r2);\n" \
"ulong tempMul3 = mmmNtoRes(((ulong) 123454321) % thePrime, thePrime, _q, _r2);\n" \
"ulong tempSub3 = mmmNtoRes(((ulong) 1110988902222) % thePrime, thePrime, _q, _r2);\n" \
"ulong resInvmod3 = mmmNtoRes(invmod3, thePrime, _q, _r2);\n" \
"ulong tempMul4 = mmmNtoRes(((ulong) 12345654321) % thePrime, thePrime, _q, _r2);\n" \
"ulong tempSub4 = mmmNtoRes(((ulong) 1111098889022222) % thePrime, thePrime, _q, _r2);\n" \
"ulong resInvmod4 = mmmNtoRes(invmod4, thePrime, _q, _r2);\n" \
"ulong resTemp = mmmPowmod(resTenE1, 179, thePrime, _q, _one);\n" \
"ulong resC = mmmMulmod(resTemp, tempMul1, thePrime, _q);\n" \
"resC = mmmSub(resC, tempSub1, thePrime);\n" \
"resTemp = mmmPowmod(resTenE1, 2699, thePrime, _q, _one);\n" \
"resC = mmmMulmod(resC, resTemp, thePrime, _q);\n" \
"resC = mmmSub(resC, tempSub2, thePrime);\n" \
"resC = mmmMulmod(resC, resInvmod2, thePrime, _q);\n" \
"resTemp = mmmPowmod(resTenE1, 35999, thePrime, _q, _one);\n" \
"resC = mmmMulmod(resC, resTemp, thePrime, _q);\n" \
"resC = mmmMulmod(resC, tempMul3, thePrime, _q);\n" \
"resC = mmmSub(resC, tempSub3, thePrime);\n" \
"resC = mmmMulmod(resC, resInvmod3, thePrime, _q);\n" \
"resC = mmmMulmod(resC, tempMul4, thePrime, _q);\n" \
"resTemp = mmmPowmod(resTenE1, 449999, thePrime, _q, _one);\n" \
"resC = mmmMulmod(resC, resTemp, thePrime, _q);\n" \
"resC = mmmSub(resC, tempSub4, thePrime);\n" \
"resC = mmmMulmod(resC, resInvmod4, thePrime, _q);\n" \
"uint   exp = 6*terms[0] - 599989;\n" \
"ulong  m = ((ulong) terms[0] * 999999) + 1000000;\n" \
"resTemp = mmmPowmod(resTenE1, exp, thePrime, _q, _one);\n" \
"resC = mmmMulmod(resC, resTemp, thePrime, _q);\n" \
"ulong resM = mmmNtoRes(m, thePrime, _q, _r2);\n" \
"ulong res9s = mmmNtoRes(999999, thePrime, _q, _r2);\n" \
"ulong resR = mmmNtoRes(((ulong) 1000000000000000000) % thePrime, thePrime, _q, _r2);\n" \
"ulong resT[6];\n" \
"resT[1] = mmmMulmod(resR, resR, thePrime, _q);\n" \
"resT[2] = mmmMulmod(resT[1], resT[1], thePrime, _q);\n" \
"resT[3] = mmmMulmod(resT[1], resT[2], thePrime, _q);\n" \
"resT[4] = mmmMulmod(resT[2], resT[2], thePrime, _q);\n" \
"resT[5] = mmmMulmod(resT[2], resT[3], thePrime, _q);\n" \
"if (resC == resM)\n" \
"collect_factor(terms[0], thePrime, factorCount, factors);\n" \
"uint idx = 1;\n" \
"while (terms[idx] > 0)\n" \
"{\n" \
"uint dn = terms[idx] - terms[idx-1];\n" \
"if (dn <= 30)\n" \
"resC = mmmMulmod(resC, resT[dn/6], thePrime, _q);\n" \
"else\n" \
"{\n" \
"resTemp = mmmPowmod(resT[1], dn/6, thePrime, _q, _one);\n" \
"resC = mmmMulmod(resC, resTemp, thePrime, _q);\n" \
"}\n" \
"resTemp = mmmNtoRes(dn, thePrime, _q, _r2);\n" \
"resTemp = mmmMulmod(resTemp, res9s, thePrime, _q);\n" \
"resM = mmmAdd(resM, resTemp, thePrime);\n" \
"if (resC == resM)\n" \
"collect_factor(terms[idx], thePrime, factorCount, factors);\n" \
"idx++;\n" \
"}\n" \
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
"ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2)\n" \
"{\n" \
"return mmmMulmod(n, _r2, _p, _q);\n" \
"}\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q)\n" \
"{\n" \
"return mmmMulmod(res, 1, _p, _q);\n" \
"}\n" \
"ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one)\n" \
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
"void collect_factor(uint    n,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorsCount,\n" \
"__global ulong2 *factors)\n" \
"{\n" \
"uint old = atomic_inc(factorsCount);\n" \
"if (old >= D_MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = n;\n" \
"factors[old].y = p;\n" \
"}\n" \
;
#endif
#endif
