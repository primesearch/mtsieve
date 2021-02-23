#ifdef __cplusplus
#ifndef _GFN_KERNEL_CL
#define _GFN_KERNEL_CL
const char *gfn_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n" \
"void  collectFactor(ulong k, uint n, ulong p,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors);\n" \
"uint  getSmallDivisor(ulong k, uint n);\n" \
"ulong mmmInvert(ulong _p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2);\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q);\n" \
"#ifndef D_MULTI_PASS\n" \
"__kernel void gfn_kernel(__global const ulong  *primes,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors)\n" \
"{\n" \
"int    gid = get_global_id(0);\n" \
"uint   n;\n" \
"ulong  p = primes[gid];\n" \
"ulong  k = (1+p) >> 1;\n" \
"ushort bitsToShift;\n" \
"ulong _q = mmmInvert(p);\n" \
"ulong _one = mmmOne(p);\n" \
"ulong _r2 = mmmR2(p, _q, _one);\n" \
"ulong mpK = mmmNtoRes(k, p, _q, _r2);\n" \
"ulong mpRem = mmmPowmod(mpK, N_MIN, p, _q, _one, _r2);\n" \
"k = p - mmmResToN(mpRem, p, _q);\n" \
"n = N_MIN;\n" \
"while (n <= N_MAX)\n" \
"{\n" \
"bitsToShift = 63 - clz(k & -k);\n" \
"k >>= bitsToShift;\n" \
"n += bitsToShift;\n" \
"if (k <= K_MAX && k >= K_MIN && n <= N_MAX)\n" \
"collectFactor(k, n, p, factorCount, factors);\n" \
"k += p;\n" \
"}\n" \
"}\n" \
"#else\n" \
"__kernel void gfn_kernel(__global const ulong  *primes,\n" \
"__global const ulong  *params,\n" \
"__global       ulong  *rems,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors)\n" \
"{\n" \
"int    gid = get_global_id(0);\n" \
"uint   n, steps = 0;\n" \
"ulong  p = primes[gid];\n" \
"ulong  k = (1+p) >> 1;\n" \
"ushort bitsToShift;\n" \
"if (params[0] == N_MIN)\n" \
"{\n" \
"ulong _q = mmmInvert(p);\n" \
"ulong _one = mmmOne(p);\n" \
"ulong _r2 = mmmR2(p, _q, _one);\n" \
"ulong mpK = mmmNtoRes(k, p, _q, _r2);\n" \
"ulong mpRem = mmmPowmod(mpK, N_MIN, p, _q, _one, _r2);\n" \
"k = p - mmmResToN(mpRem, p, _q);\n" \
"n = N_MIN;\n" \
"}\n" \
"else\n" \
"{\n" \
"k = rems[gid];\n" \
"n = params[0];\n" \
"}\n" \
"while (n <= N_MAX - 64 && steps < D_MAX_STEPS - 64)\n" \
"{\n" \
"bitsToShift = 63 - clz(k & -k);\n" \
"k >>= bitsToShift;\n" \
"n += bitsToShift;\n" \
"steps += bitsToShift;\n" \
"if (k <= K_MAX && k >= K_MIN)\n" \
"collectFactor(k, n, p, factorCount, factors);\n" \
"k += p;\n" \
"}\n" \
"while (n <= N_MAX && steps < D_MAX_STEPS)\n" \
"{\n" \
"if (k & 1)\n" \
"{\n" \
"if (k <= K_MAX && k >= K_MIN && n <= N_MAX)\n" \
"collectFactor(k, n, p, factorCount, factors);\n" \
"k += p;\n" \
"}\n" \
"k >>= 1;\n" \
"n++;\n" \
"steps++;\n" \
"}\n" \
"rems[gid] = k;\n" \
"}\n" \
"#endif\n" \
"void  collectFactor(ulong k, uint n, ulong p,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors)\n" \
"{\n" \
"if (getSmallDivisor(k, n) > 0)\n" \
"return;\n" \
"int old = atomic_inc(factorCount);\n" \
"if (old >= D_MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = k;\n" \
"factors[old].y = n;\n" \
"factors[old].z = p;\n" \
"}\n" \
"uint  getSmallDivisor(ulong k, uint n)\n" \
"{\n" \
"uint  smallN;\n" \
"ulong smallK;\n" \
"smallN = n % (2);\n" \
"smallK = k % (3);\n" \
"if ((smallK << smallN) % (3) == 2) return 3;\n" \
"smallN = n % (4);\n" \
"smallK = k % (5);\n" \
"if ((smallK << smallN) % (5) == 4) return 5;\n" \
"smallN = n % (6);\n" \
"smallK = k % (7);\n" \
"if ((smallK << smallN) % (7) == 6) return 7;\n" \
"smallN = n % (10);\n" \
"smallK = k % (11);\n" \
"if ((smallK << smallN) % (11) == 10) return 11;\n" \
"smallN = n % (12);\n" \
"smallK = k % (13);\n" \
"if ((smallK << smallN) % (13) == 12) return 13;\n" \
"smallN = n % (16);\n" \
"smallK = k % (17);\n" \
"if ((smallK << smallN) % (17) == 16) return 17;\n" \
"smallN = n % (18);\n" \
"smallK = k % (19);\n" \
"if ((smallK << smallN) % (19) == 18) return 19;\n" \
"smallN = n % (22);\n" \
"smallK = k % (23);\n" \
"if ((smallK << smallN) % (23) == 22) return 23;\n" \
"smallN = n % (28);\n" \
"smallK = k % (29);\n" \
"if ((smallK << smallN) % (29) == 28) return 29;\n" \
"smallN = n % (30);\n" \
"smallK = k % (31);\n" \
"if ((smallK << smallN) % (31) == 30) return 31;\n" \
"smallN = n % (36);\n" \
"smallK = k % (37);\n" \
"if ((smallK << smallN) % (37) == 36) return 37;\n" \
"smallN = n % (40);\n" \
"smallK = k % (41);\n" \
"if ((smallK << smallN) % (41) == 40) return 41;\n" \
"smallN = n % (42);\n" \
"smallK = k % (43);\n" \
"if ((smallK << smallN) % (43) == 42) return 43;\n" \
"smallN = n % (46);\n" \
"smallK = k % (47);\n" \
"if ((smallK << smallN) % (47) == 46) return 47;\n" \
"return 0;\n" \
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
;
#endif
#endif
