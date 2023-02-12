#include "scfdhmpi.h"
// #include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
// #include <cstring>
#include <cassert>

#include <x86intrin.h>

static inline __m256d rcp_full(const __m256d x){
	__m256d one = _mm256_set1_pd(1.0);
	__m256d two = _mm256_set1_pd(2.0);

	__m256d y = _mm256_cvtps_pd(
			_mm_rcp_ps(
				_mm256_cvtpd_ps(x)));

	__m256d h    = _mm256_fnmadd_pd(x, y, one);
	__m256d hp1  = _mm256_fnmadd_pd(x, y, two);
	__m256d yh   = _mm256_mul_pd(y, h);
	__m256d h2p1 = _mm256_fmadd_pd(h, h, one);
	__m256d poly = _mm256_mul_pd(hp1, h2p1);

	y = _mm256_fmadd_pd(poly, yh, y);  // y + yh * poly

#if 0 // for safety
	y = _mm256_div_pd(one, x);
#endif
	return y;
}

static inline __m256d rsqrt_full(const __m256d x){
	__m256d one = _mm256_set1_pd(1.0);
	__m256d a   = _mm256_set1_pd(1./2.);
	__m256d b   = _mm256_set1_pd(3./8.);
	__m256d c   = _mm256_set1_pd(5./16.);
	__m256d d   = _mm256_set1_pd(35./128.);


	__m256d y = _mm256_cvtps_pd(
			_mm_rsqrt_ps(
				_mm256_cvtpd_ps(x)));
	__m256d y2 = _mm256_mul_pd(y, y);
	__m256d h  = _mm256_fnmadd_pd(x, y2, one);
	__m256d yh = _mm256_mul_pd(y, h);

	__m256d poly = _mm256_fmadd_pd(h, d, c); // c + d*h
	poly = _mm256_fmadd_pd(h, poly, b); // b + h(c + dh)
	poly = _mm256_fmadd_pd(h, poly, a); // a + h(b + h(c + dh))

	y = _mm256_fmadd_pd(poly, yh, y);  // y + yh * poly

#if 0 // for safety
	y = _mm256_div_pd(one, _mm256_sqrt_pd(x));
#endif
	return y;
}

static inline void hadd_accum(double *acc, const __m256d s, const __m256d c){
	__m256d sum = _mm256_hadd_pd(s, c);
	__m128d xl = _mm256_extractf128_pd(sum, 0);
	__m128d xh = _mm256_extractf128_pd(sum, 1);
	__m128d xsum = _mm_add_pd(xl, xh);
	__m128d a0 = _mm_load_pd(acc);
	xsum = _mm_add_pd(xsum, a0);
	_mm_store_pd(acc, xsum);
}

static inline void hadd_accum2(double *acc, const __m256d a, const __m256d b, const __m256d c, const __m256d d){
	__m256d ab = _mm256_hadd_pd(a, b);
	__m256d cd = _mm256_hadd_pd(c, d);
	__m256d lo = _mm256_permute2f128_pd(ab, cd, 0x20);
	__m256d hi = _mm256_permute2f128_pd(ab, cd, 0x31);

	__m256d sum = _mm256_load_pd(acc);
	sum = _mm256_add_pd(sum, _mm256_add_pd(lo, hi));
	_mm256_store_pd(acc, sum);
}

template<int n, typename I>
auto align_up(const I i) -> I {
	return i%n ? i/n*n+n : i;
}

#define SCF_BASE_NAME CB
#define SCF_BASE_NAME_is_CB
#include "scfcoef-inc.h"
#undef SCF_BASE_NAME 
#undef SCF_BASE_NAME_is_CB
#define SCF_BASE_NAME LH
#define SCF_BASE_NAME_is_LH
#include "scfcoef-inc.h"
