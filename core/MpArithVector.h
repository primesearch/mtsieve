/* MpArithVector.h -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   Thanks to Yves Gallot for this implementation based upon 
   Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
*/

#ifndef _MpArithVector_H
#define _MpArithVector_H

#define		VECTOR_SIZE		4		// must be a power of two

// Arithmetic on vectors: hide the latency of the MUL instruction.

// Montgomery form: if 0 <= a < p then r is 2^64 * a mod p
template <size_t N>
class MpResVector
{
private:
	uint64_t _r[N];

public:
	uint64_t operator [](const size_t i) const { return _r[i]; }
	uint64_t & operator [](const size_t i) { return _r[i]; }
};

// Montgomery modular arithmetic in Z/pZ
template <size_t N>
class MpArithVector
{
private:
	uint64_t _p[N], _q[N];
	MpResVector<N> _one;		// 2^64 mod p
	MpResVector<N> _r2;		// (2^64)^2 mod p

private:
	// p * p_inv = 1 (mod 2^64) (Newton's method)
	constexpr uint64_t invert(const uint64_t p) const
	{
		uint64_t p_inv = 1, prev = 0;
      while (p_inv != prev)
      {
         prev = p_inv;
         p_inv *= 2 - p * p_inv;
      }
		return p_inv;
	}

	// The Montgomery REDC algorithm
	constexpr uint64_t REDC(const __uint128_t t, const uint64_t p, const uint64_t q) const
	{
		const uint64_t m = uint64_t(t) * q;
		const int64_t r = int64_t((t >> 64) - uint64_t((m * __uint128_t(p)) >> 64));
		return (r < 0) ? uint64_t(r + p) : uint64_t(r);
	}

public:
	MpArithVector(const uint64_t * const p)
	{
		for (size_t k = 0; k < N; ++k)
		{
			const uint64_t p_k = p[k];
			_p[k] = p_k;
			_q[k] = invert(p_k);
			_one[k] = (-p_k) % p_k;
		}
		
		MpResVector<N> t = add(_one, _one); t = add(t, t);	// 4
		for (size_t i = 0; i < 5; ++i) t = mul(t, t);	// 4^{2^5} = 2^64
		_r2 = t;
	}

	static MpResVector<N> zero()
	{
		MpResVector<N> r;
		for (size_t k = 0; k < N; ++k) r[k] = 0;
		return r;
	}

	MpResVector<N> one() const { return _one; }	// Montgomery form of 1

	uint64_t p(size_t k) const { return _p[k]; }
   
	static bool at_least_one_is_equal(const MpResVector<N> & a, const MpResVector<N> & b)
	{
		bool is_equal = false;
		for (size_t k = 0; k < N; ++k) is_equal |= (a[k] == b[k]);
		return is_equal;
	}

	MpResVector<N> add(const MpResVector<N> & a, const MpResVector<N> & b) const
	{
		MpResVector<N> r;
		for (size_t k = 0; k < N; ++k)
		{
			const uint64_t c = (a[k] >= _p[k] - b[k]) ? _p[k] : 0;
			r[k] = a[k] + b[k] - c;
		}
		return r;
	}

	MpResVector<N> sub(const MpResVector<N> & a, const MpResVector<N> & b) const
	{
		MpResVector<N> r;
		for (size_t k = 0; k < N; ++k)
		{
			const uint64_t c = (a[k] < b[k]) ? _p[k] : 0;
			r[k] = a[k] - b[k] + c;
		}
		return r;
	}
   
	uint64_t mul(const uint64_t a, const MpResVector<N> & b, size_t k) const
	{
	   return REDC(a * __uint128_t(b[k]), _p[k], _q[k]);
	}
  
   uint64_t mul(const MpResVector<N> & a, const MpResVector<N> & b, size_t k) const
	{
	   return REDC(a[k] * __uint128_t(b[k]), _p[k], _q[k]);
	}
   
	MpResVector<N> mul(const MpResVector<N> & a, const MpResVector<N> & b) const
	{
		MpResVector<N> r;
		for (size_t k = 0; k < N; ++k)
		{
			r[k] = REDC(a[k] * __uint128_t(b[k]), _p[k], _q[k]);
		}
		return r;
	}
   
	MpResVector<N> pow(const MpResVector<N> & a, size_t exp) const
	{
      MpResVector<N> x = a;
      MpResVector<N> y = _one;

      while (true)
      {
         if (exp & 1)
            y = mul(x, y);

         exp >>= 1;

         if (!exp)
            break;

         x = mul(x, x);
      }
      
      return y;
	}

	// Convert n to Montgomery representation
	MpResVector<N> nToRes(const uint64_t *n) const
	{
		// n * (2^64)^2 = (n * 2^64) * (1 * 2^64)
		MpResVector<N> r;
      
		for (size_t k = 0; k < N; ++k)
         r[k] = n[k];
         
		return mul(r, _r2);
	}
   
	// Convert n to Montgomery representation
	MpResVector<N> nToRes(uint64_t n) const
	{
		// n * (2^64)^2 = (n * 2^64) * (1 * 2^64)
		MpResVector<N> r;
      
		for (size_t k = 0; k < N; ++k)
         r[k] = n;
         
		return mul(r, _r2);
	}
   
	// Convert n to Montgomery representation
	MpResVector<N> nToRes(uint32_t *n) const
	{
		// n * (2^64)^2 = (n * 2^64) * (1 * 2^64)
		MpResVector<N> r;
      
		for (size_t k = 0; k < N; ++k)
         r[k] = n[k];
         
		return mul(r, _r2);
	}
   
   // Convert Montgomery representation to n
	MpResVector<N> resToN(const MpResVector<N> & a) const
	{
		MpResVector<N> r;
		for (size_t k = 0; k < N; ++k)
		{
			r[k] = REDC(a[k], _p[k], _q[k]);
		}
		return r;
	}
};


typedef MpResVector<VECTOR_SIZE> MpResVec;
typedef MpArithVector<VECTOR_SIZE> MpArithVec;

#endif