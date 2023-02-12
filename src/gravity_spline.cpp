#include <cmath>
#include "user_defined.hpp"
#include <x86intrin.h>

inline double pcut(double u0){
	using F = decltype(u0);
#if 0
	using std::min;
	using std::max;
#elif 1
	auto min = fmin;
	auto max = fmax;
#else
	auto min = [](auto a, auto b) { return a<b ? a : b; };
	auto max = [](auto a, auto b) { return a>b ? a : b; };
#endif

	double u  = min(F(2), u0);
	double u2 = u*u;

	double pl = fma(u, -3, 9);
	pl = fma(pl, u2, -20);
	pl = fma(pl, u2, 42);
	pl *= u;
	
	double pr = fma(u, 4, 2);

	double usub = max(F(0), u-1);
	double usub2 = usub  * usub;
	double usub4 = usub2 * usub2;
	double usub5 = usub  * usub4;

	double sum = fma(usub5, pr, pl);

	return F(1./30.) * sum;
} 

inline double acut(double u0){
	using F = decltype(u0);
#if 0
	using std::min;
	using std::max;
#elif 1
	auto min = fmin;
	auto max = fmax;
#else
	auto min = [](auto a, auto b) { return a<b ? a : b; };
	auto max = [](auto a, auto b) { return a>b ? a : b; };
#endif
	double u  = min(F(2), u0);
	double u2 = u*u;
	double u3 = u*u2;

	double pl = fma(u, 15, -36);
	pl = fma(pl, u2, 40);
	pl *= u3;

	double pr = fma(u, 20, 8);
	pr = fma(pr, u, 2);

	double usub = max(F(0), u-1);
	double usub2 = usub  * usub;
	double usub4 = usub2 * usub2;

	double sum = fma(usub4, -pr, pl);

	return F(1./30.) * sum;
}

#if 0
struct Body{
	double x, y, z, m;
};

struct Acceleration{
	double ax, ay, az, pot;
};

__attribute__((noinline))
void nbody_spline_ref(
	const int n,
	const double /* rcut2 */, // = 4 eps^2
	const double epsinv,
	const Body body[],
	Acceleration acc[])
{
	for(int i=0; i<n; i++){ 
		const double xi=body[i].x, yi=body[i].y, zi=body[i].z;
		double ax=0, ay=0, az=0, pot=0;

		for(int j=0; j<n; j++){
			if(j == i) continue;

			double dx = body[j].x - xi;
			double dy = body[j].y - yi;
			double dz = body[j].z - zi;

			double r2 = dx*dx + dy*dy + dz*dz;

			double ri = 1.0 / sqrt(r2);

			double mri = body[j].m * ri;
			double ri2 = ri * ri;

			double mri3 = mri * ri2;

			double u = r2 * ri * epsinv;

			mri3 *= acut(u);
			mri  *= pcut(u);

			ax += mri3 * dx;
			ay += mri3 * dy;
			az += mri3 * dz;
			pot -= mri;
		}
		assert(std::isfinite(pot));
		acc[i] = {ax, ay, az, pot};
	}
	return ;
}
#endif

#if 0
__attribute__((noinline))
void nbody_spline_list(
	const int n,
	const double rcut2, // = 4 eps^2
	const double epsinv,
	const Body body[],
	Acceleration acc[])
{
	enum{ LEN = 64, };
	int list[LEN];

	int len_max=0, len_sum=0;
	for(int i=0; i<n; i++){ 
		int num = 0;
		const double xi=body[i].x, yi=body[i].y, zi=body[i].z;
		double ax=0, ay=0, az=0, pot=0;

		for(int j=0; j<n; j++){
			double dx = body[j].x - xi;
			double dy = body[j].y - yi;
			double dz = body[j].z - zi;

			double r2 = dx*dx + dy*dy + dz*dz;

			double ri = 1.0 / sqrt(r2);

			double mri = body[j].m * ri;
			double ri2 = ri * ri;

			double mri3 = mri * ri2;

			if(r2 < rcut2){
				list[num++%LEN] = j;
				mri3 = 0.0;
			}

			ax += mri3 * dx;
			ay += mri3 * dy;
			az += mri3 * dz;
			pot -= mri3 * r2;
		}
		assert(std::isfinite(pot));
		// printf("%d : %d\n", i, num);
		assert(num <= LEN);

		len_max = std::max(len_max, num);
		len_sum += num;
		for(int jj=0; jj<num; jj++){
			int j = list[jj];

			double dx = body[j].x - xi;
			double dy = body[j].y - yi;
			double dz = body[j].z - zi;

			double r2 = dx*dx + dy*dy + dz*dz;

			if(0.0 == r2) continue;

			double ri = 1.0 / sqrt(r2);

			double mri = body[j].m * ri;
			double ri2 = ri * ri;

			double mri3 = mri * ri2;
			double u = r2 * ri * epsinv;

			mri3 *= acut(u);
			mri  *= pcut(u);

			ax += mri3 * dx;
			ay += mri3 * dy;
			az += mri3 * dz;
			pot -= mri;

		}
		assert(std::isfinite(pot));
		acc[i] = {ax, ay, az, pot};
	}
	printf("max: %d, avr: %f\n", len_max, len_sum * (1./n));
	return ;
}
#endif

inline void transpose_4ymm_pd(__m256d &r0, __m256d &r1, __m256d &r2, __m256d &r3){
	__m256d tmp0 = _mm256_unpacklo_pd(r0, r1); // | r12 | r02 | r10 | r00 |
	__m256d tmp1 = _mm256_unpackhi_pd(r0, r1); // | r13 | r03 | r11 | r01 |
	__m256d tmp2 = _mm256_unpacklo_pd(r2, r3); // | r32 | r22 | r30 | r20 |
	__m256d tmp3 = _mm256_unpackhi_pd(r2, r3); // | r33 | r23 | r31 | r21 |
	r0 = _mm256_permute2f128_pd(tmp0, tmp2, (0)+(2<<4));
	r1 = _mm256_permute2f128_pd(tmp1, tmp3, (0)+(2<<4));
	r2 = _mm256_permute2f128_pd(tmp0, tmp2, (1)+(3<<4));
	r3 = _mm256_permute2f128_pd(tmp1, tmp3, (1)+(3<<4));
}

inline __m256d rsqrtCubed_x5(
		const __m256d x,
		const __m256d m)
{
	const __m256d one  = _mm256_set1_pd(1.0);
	const __m256d a  = _mm256_set1_pd(  3./ 2.);
	const __m256d b  = _mm256_set1_pd( 15./ 8.);
	const __m256d c  = _mm256_set1_pd( 35./ 16.);
	const __m256d d  = _mm256_set1_pd(315./128.);

	__m256d y = _mm256_cvtps_pd(
			_mm_rsqrt_ps(
				_mm256_cvtpd_ps(x)));

	__m256d y2  = _mm256_mul_pd(y, y);
	__m256d my  = _mm256_mul_pd(m, y);

	__m256d z   = _mm256_mul_pd(my, y2);
	__m256d h   = _mm256_fnmadd_pd(x,  y2, one); // c - a * b

	__m256d poly = _mm256_fmadd_pd(d, h, c); // d*h + c
	poly = _mm256_fmadd_pd(poly, h, b); // (d*h + c)*h + b
	poly = _mm256_fmadd_pd(poly, h, a); // ((d*h + c)*h + b)*h + a

	__m256d zh  = _mm256_mul_pd(z, h);

	__m256d z1  = _mm256_fmadd_pd(zh, poly, z);

	return z1;
}

__attribute__((noinline))
void nbody_m256d_spline(
	const double eps2_sd, // = 4 eps^2
	const double epsinv,
	const EPI epi[],
	const int ni,
	const EPJ epj[],
	const int nj,
	Force force[])
{
	auto max = [](auto a, auto b) { return a>b ? a : b; };
	int max_length = 0;

	enum{ LEN = 256, };
	int list[4][LEN] __attribute__((aligned(64)));

	const __m256d rcut2 = _mm256_set1_pd(4.0 * eps2_sd);

	for(int i=0; i<ni; i+=4){
		__m256d xi = _mm256_loadu_pd(&epi[i+0].pos.x);
		__m256d yi = _mm256_loadu_pd(&epi[i+1].pos.x);
		__m256d zi = _mm256_loadu_pd(&epi[i+2].pos.x);
		__m256d mi = _mm256_loadu_pd(&epi[i+3].pos.x);
		transpose_4ymm_pd(xi, yi, zi, mi);

		__m256d ax, ay, az, pot;
		ax = ay = az = pot = _mm256_set1_pd(0);

		int n0, n1, n2, n3;
		n0 = n1 = n2 = n3 = 0;

		for(int j=0; j<nj; j+=1){
			__m256d dx = _mm256_sub_pd(xi, _mm256_set1_pd(epj[j+0].pos.x));
			__m256d dy = _mm256_sub_pd(yi, _mm256_set1_pd(epj[j+0].pos.y));
			__m256d dz = _mm256_sub_pd(zi, _mm256_set1_pd(epj[j+0].pos.z));
			__m256d mj = _mm256_set1_pd(epj[j+0].mass);

			__m256d r2  = 
				_mm256_fmadd_pd(dz, dz, 
					_mm256_fmadd_pd(dy, dy,
						_mm256_mul_pd(dx, dx)));

		       __m256d mask = _mm256_cmp_pd(r2, rcut2, _CMP_LT_OQ);

		       __m256d mri3 = rsqrtCubed_x5(r2, mj);

		       mri3 = _mm256_andnot_pd(mask, mri3);

		       ax  = _mm256_fnmadd_pd(mri3, dx, ax);
		       ay  = _mm256_fnmadd_pd(mri3, dy, ay);
		       az  = _mm256_fnmadd_pd(mri3, dz, az);
		       pot = _mm256_fnmadd_pd(mri3, r2, pot);

		       int m = _mm256_movemask_pd(mask);
		       if(m){
			       if(m&1){
				       list[0][n0++%LEN] = j;
			       }
			       if(m&2){
				       list[1][n1++%LEN] = j;
			       }
			       if(m&4){
				       list[2][n2++%LEN] = j;
			       }
			       if(m&8){
				       list[3][n3++%LEN] = j;
			       }
		       }
		}

		transpose_4ymm_pd(ax, ay, az, pot);
		if(i+0 < ni) _mm256_storeu_pd(&force[i+0].acc.x, ax);
		if(i+1 < ni) _mm256_storeu_pd(&force[i+1].acc.x, ay);
		if(i+2 < ni) _mm256_storeu_pd(&force[i+2].acc.x, az);
		if(i+3 < ni) _mm256_storeu_pd(&force[i+3].acc.x, pot);

		int ns[4] = {n0, n1, n2, n3};
		for(int ii=0; ii<4; ii++){
			if(i+ii >= ni) break;
			int num = ns[ii];
			// printf("%d  ", num);
			max_length = max(max_length, num);
			assert(num <= LEN);

			double xi = epi[i+ii].pos.x;
			double yi = epi[i+ii].pos.y;
			double zi = epi[i+ii].pos.z;

			double ax=0, ay=0, az=0, pot=0;

			for(int jj=0; jj<num; jj++){
				int j = list[ii][jj];

				double dx = epj[j].pos.x - xi;
				double dy = epj[j].pos.y - yi;
				double dz = epj[j].pos.z - zi;

				double r2 = dx*dx + dy*dy + dz*dz;

				if(0.0 == r2) continue;

				double ri = 1.0 / sqrt(r2);

				double mri = epj[j].mass * ri;
				double ri2 = ri * ri;

				double mri3 = mri * ri2;
				double u = r2 * ri * epsinv;

				mri3 *= acut(u);
				mri  *= pcut(u);

				ax += mri3 * dx;
				ay += mri3 * dy;
				az += mri3 * dz;
				pot -= mri;
			}
			force[i+ii].acc.x  += ax;
			force[i+ii].acc.y  += ay;
			force[i+ii].acc.z  += az;
			force[i+ii].pot    += pot;
		}
	}
	// puts("");
	TreeProf::max_list_length = max(TreeProf::max_list_length, max_length);
}



#if 0
void prlong(long l){
	l <<= 32;
	double d;
	memcpy(&d, &l, 8);
	printf("%lx : %f\n", l>>32, d);
}

int main(){
	enum { N = 20, };

	prlong(1072693248);
	prlong(1076101120);
	prlong(1077149696);
	prlong(1073741824);
	prlong(-1067941888);
	prlong(1080033280);

	for(int i=0; i<N+3; i++){
		double u = i * (2.0/N);
		double f = pcut(u);
		double g = acut(u);

		printf("%f\t%f\t%f\t%A\n", u, f, f*30, f*30);
		printf("%f\t%f\t%f\t%A\n", u, g, g*30, g*30);
	}
}
#endif

