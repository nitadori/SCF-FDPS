#include <stdio.h>
#include <x86intrin.h>

inline void hadd_accum2(double *acc, const __m256d a, const __m256d b, const __m256d c, const __m256d d){
	__m256d ab = _mm256_hadd_pd(a, b);
	__m256d cd = _mm256_hadd_pd(c, d);
	__m256d lo = _mm256_permute2f128_pd(ab, cd, 0x20);
	__m256d hi = _mm256_permute2f128_pd(ab, cd, 0x31);

	__m256d sum = _mm256_load_pd(acc);
	sum = _mm256_add_pd(sum, _mm256_add_pd(lo, hi));
	_mm256_store_pd(acc, sum);
}

int main(){
	double a[5][4] = {
		{ 1.,  2.,  3.,  4.},
		{ 5.,  6.,  7.,  8.},
		{ 9., 10., 11., 12.},
		{13., 14., 15., 16.},
		{1000, 2000, 3000, 4000}
	};
	hadd_accum2(a[4], *(__m256d*)a[0], *(__m256d*)a[1], *(__m256d*)a[2], *(__m256d*)a[3]); 
	printf("%f %f %f %f\n", a[4][0], a[4][1], a[4][2], a[4][3]);
}
