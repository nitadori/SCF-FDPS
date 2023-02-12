#include "user_defined.hpp"
#include <x86intrin.h>

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

inline void rsqrt_nr_2ymm_pd(__m256d &x0, __m256d &x1){
	const __m256 half = _mm256_set1_ps(0.5f);
	const __m256 c1p5 = _mm256_set1_ps(1.5f);

	__m128 lo = _mm256_cvtpd_ps(x0);
	__m128 hi = _mm256_cvtpd_ps(x1);

	__m256 x = _mm256_set_m128(hi, lo);
	__m256 y = _mm256_rsqrt_ps(x);
	__m256 xh = _mm256_mul_ps(x, half); 
	__m256 y2 = _mm256_mul_ps(y, y); 
	__m256 p  = _mm256_fnmadd_ps(xh, y2, c1p5);

	y = _mm256_mul_ps(y, p); 

	lo = _mm256_extractf128_ps(y, 0);
	hi = _mm256_extractf128_ps(y, 1);
	
	x0 = _mm256_cvtps_pd(lo);
	x1 = _mm256_cvtps_pd(hi);
}

inline __m256d rsqrtCubed(
		const __m256d x,
		const __m256d y,
		const __m256d m)
{
	const __m256d one  = _mm256_set1_pd(1.0);
	const __m256d a    = _mm256_set1_pd(3./2.);
	const __m256d b    = _mm256_set1_pd(15./8.);

	__m256d y2  = _mm256_mul_pd(y, y);
	__m256d my  = _mm256_mul_pd(m, y);

	__m256d z   = _mm256_mul_pd(my, y2);
	__m256d h   = _mm256_fnmadd_pd(x,  y2, one); // c - a * b
	__m256d abh = _mm256_fmadd_pd(b, h, a); // a + b*h

	__m256d zh  = _mm256_mul_pd(z, h);

	// abh = a;  // force 2nd order

	__m256d z1  = _mm256_fmadd_pd(zh, abh, z);

	return z1;
}

// __attribute__((noinline))
template<typename Tepj>
void nbody_m256d(
	const double eps2_sd,
	const double epsinv_sd,
	const EPI epi[],
	const int ni,
	const Tepj epj[],
	const int nj,
	Force force[])
{
	const __m256d eps2 = _mm256_set1_pd(eps2_sd);
	const __m256d epsinv = _mm256_set1_pd(epsinv_sd);

	for(int i=0; i<ni; i+=4){
		__m256d xi = _mm256_loadu_pd(&epi[i+0].pos.x);
		__m256d yi = _mm256_loadu_pd(&epi[i+1].pos.x);
		__m256d zi = _mm256_loadu_pd(&epi[i+2].pos.x);
		__m256d mi = _mm256_loadu_pd(&epi[i+3].pos.x);
		transpose_4ymm_pd(xi, yi, zi, mi);

		__m256d ax, ay, az, pot;
		ax = ay = az = pot = _mm256_set1_pd(0);
		if constexpr ( std::is_same<Tepj, EPJ>{}){
			pot = mi * epsinv; // cancel self interaction
		}
		// pot = ax;

		/*
		if(nj%2){
			(*(EPJ*)&epj[nj]).pos  = 1.0;
			(*(EPJ*)&epj[nj]).mass = 0.0;
		}
		*/
		const int nj2 = nj/2*2;

		for(int j=0; j<nj2; j+=2){
			__m256d dx_0 = _mm256_sub_pd(xi, _mm256_set1_pd(epj[j+0].pos.x));
			__m256d dy_0 = _mm256_sub_pd(yi, _mm256_set1_pd(epj[j+0].pos.y));
			__m256d dz_0 = _mm256_sub_pd(zi, _mm256_set1_pd(epj[j+0].pos.z));
			__m256d mj_0 = _mm256_set1_pd(epj[j+0].mass);

			__m256d r2_0  = 
				_mm256_fmadd_pd(dz_0, dz_0, 
					_mm256_fmadd_pd(dy_0, dy_0,
						_mm256_fmadd_pd(dx_0, dx_0, eps2)));

			__m256d dx_1 = _mm256_sub_pd(xi, _mm256_set1_pd(epj[j+1].pos.x));
			__m256d dy_1 = _mm256_sub_pd(yi, _mm256_set1_pd(epj[j+1].pos.y));
			__m256d dz_1 = _mm256_sub_pd(zi, _mm256_set1_pd(epj[j+1].pos.z));
			__m256d mj_1 = _mm256_set1_pd(epj[j+1].mass);

			__m256d r2_1  = 
				_mm256_fmadd_pd(dz_1, dz_1, 
					_mm256_fmadd_pd(dy_1, dy_1,
						_mm256_fmadd_pd(dx_1, dx_1, eps2)));

			__m256d rinv_0 = r2_0;
			__m256d rinv_1 = r2_1;
			rsqrt_nr_2ymm_pd(rinv_0, rinv_1);

		       __m256d mri3_0 = rsqrtCubed(r2_0, rinv_0, mj_0);

		       ax  = _mm256_fnmadd_pd(mri3_0, dx_0, ax);
		       ay  = _mm256_fnmadd_pd(mri3_0, dy_0, ay);
		       az  = _mm256_fnmadd_pd(mri3_0, dz_0, az);
		       pot = _mm256_fnmadd_pd(mri3_0, r2_0, pot);

		       __m256d mri3_1 = rsqrtCubed(r2_1, rinv_1, mj_1);

		       ax  = _mm256_fnmadd_pd(mri3_1, dx_1, ax);
		       ay  = _mm256_fnmadd_pd(mri3_1, dy_1, ay);
		       az  = _mm256_fnmadd_pd(mri3_1, dz_1, az);
		       pot = _mm256_fnmadd_pd(mri3_1, r2_1, pot);
		}
		for(int j=nj2; j<nj; j+=2){
			__m256d dx_0 = _mm256_sub_pd(xi, _mm256_set1_pd(epj[j+0].pos.x));
			__m256d dy_0 = _mm256_sub_pd(yi, _mm256_set1_pd(epj[j+0].pos.y));
			__m256d dz_0 = _mm256_sub_pd(zi, _mm256_set1_pd(epj[j+0].pos.z));
			__m256d mj_0 = _mm256_set1_pd(epj[j+0].mass);

			__m256d r2_0  = 
				_mm256_fmadd_pd(dz_0, dz_0, 
					_mm256_fmadd_pd(dy_0, dy_0,
						_mm256_fmadd_pd(dx_0, dx_0, eps2)));

			__m256d rinv_0 = r2_0;
			__m256d rinv_1 = r2_0;
			rsqrt_nr_2ymm_pd(rinv_0, rinv_1);

		       __m256d mri3_0 = rsqrtCubed(r2_0, rinv_0, mj_0);

		       ax  = _mm256_fnmadd_pd(mri3_0, dx_0, ax);
		       ay  = _mm256_fnmadd_pd(mri3_0, dy_0, ay);
		       az  = _mm256_fnmadd_pd(mri3_0, dz_0, az);
		       pot = _mm256_fnmadd_pd(mri3_0, r2_0, pot);
		}

		transpose_4ymm_pd(ax, ay, az, pot);
		auto accumulate = [](double *dst, __m256d val){
			__m256d sum = _mm256_loadu_pd(dst);
			sum = _mm256_add_pd(sum, val);
			_mm256_storeu_pd(dst, sum);
		};
		if(i+0 < ni) accumulate(&force[i+0].acc.x, ax);
		if(i+1 < ni) accumulate(&force[i+1].acc.x, ay);
		if(i+2 < ni) accumulate(&force[i+2].acc.x, az);
		if(i+3 < ni) accumulate(&force[i+3].acc.x, pot);
	}
}

// explicit instanciations
template void nbody_m256d(
	const double eps2_sd,
	const double epsinv_sd,
	const EPI epi[],
	const int ni,
	const EPJ epj[],
	const int nj,
	Force force[]);

template void nbody_m256d(
	const double eps2_sd,
	const double epsinv_sd,
	const EPI epi[],
	const int ni,
	const SPJ epj[],
	const int nj,
	Force force[]);

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
void nbody_m256d_nj1(
	const double eps2_sd,
	const EPI epi[],
	const int ni,
	const ParticleSimulator::SPJQuadrupole spj[],
	const int nj,
	Force force[])
{
	const __m256d eps2 = _mm256_set1_pd(eps2_sd);

	for(int i=0; i<ni; i+=4){
		__m256d xi = _mm256_loadu_pd(&epi[i+0].pos.x);
		__m256d yi = _mm256_loadu_pd(&epi[i+1].pos.x);
		__m256d zi = _mm256_loadu_pd(&epi[i+2].pos.x);
		__m256d mi = _mm256_loadu_pd(&epi[i+3].pos.x);
		transpose_4ymm_pd(xi, yi, zi, mi);

		__m256d ax, ay, az, pot;
		ax = ay = az = pot = _mm256_set1_pd(0);

		for(int j=0; j<nj; j+=1){
			__m256d dx = _mm256_sub_pd(xi, _mm256_set1_pd(spj[j+0].pos.x));
			__m256d dy = _mm256_sub_pd(yi, _mm256_set1_pd(spj[j+0].pos.y));
			__m256d dz = _mm256_sub_pd(zi, _mm256_set1_pd(spj[j+0].pos.z));
			__m256d mj = _mm256_set1_pd(spj[j+0].mass);

			__m256d r2  = 
				_mm256_fmadd_pd(dz, dz, 
					_mm256_fmadd_pd(dy, dy,
						_mm256_fmadd_pd(dx, dx, eps2)));

		       __m256d mri3 = rsqrtCubed_x5(r2, mj);

		       ax  = _mm256_fnmadd_pd(mri3, dx, ax);
		       ay  = _mm256_fnmadd_pd(mri3, dy, ay);
		       az  = _mm256_fnmadd_pd(mri3, dz, az);
		       pot = _mm256_fnmadd_pd(mri3, r2, pot);
		}
		transpose_4ymm_pd(ax, ay, az, pot);

		auto accumulate = [](double *dst, __m256d val){
			__m256d sum = _mm256_loadu_pd(dst);
			sum = _mm256_add_pd(sum, val);
			_mm256_storeu_pd(dst, sum);
		};
		if(i+0 < ni) accumulate(&force[i+0].acc.x, ax);
		if(i+1 < ni) accumulate(&force[i+1].acc.x, ay);
		if(i+2 < ni) accumulate(&force[i+2].acc.x, az);
		if(i+3 < ni) accumulate(&force[i+3].acc.x, pot);
	}
}

inline __m256d rsqrt_pow5_x5(
		const __m256d x)
{
	const __m256d one  = _mm256_set1_pd(1.0);
	const __m256d a  = _mm256_set1_pd(   5./ 2.);
	const __m256d b  = _mm256_set1_pd(  35./ 8.);
	const __m256d c  = _mm256_set1_pd( 105./ 16.);
	const __m256d d  = _mm256_set1_pd(1155./128.);

	__m256d y = _mm256_cvtps_pd(
			_mm_rsqrt_ps(
				_mm256_cvtpd_ps(x)));

	__m256d y2  = _mm256_mul_pd(y, y);
	__m256d y4  = _mm256_mul_pd(y2, y2);

	__m256d z   = _mm256_mul_pd(y4, y);
	__m256d h   = _mm256_fnmadd_pd(x,  y2, one); // c - a * b

	__m256d poly = _mm256_fmadd_pd(d, h, c); // d*h + c
	poly = _mm256_fmadd_pd(poly, h, b); // (d*h + c)*h + b
	poly = _mm256_fmadd_pd(poly, h, a); // ((d*h + c)*h + b)*h + a

	__m256d zh  = _mm256_mul_pd(z, h);

	__m256d z1  = _mm256_fmadd_pd(zh, poly, z);

	return z1;
}

inline __m256d rsqrt_pow7_x5(
		const __m256d x)
{
	const __m256d one  = _mm256_set1_pd(1.0);
	const __m256d a  = _mm256_set1_pd(   7./ 2.);
	const __m256d b  = _mm256_set1_pd(  63./ 8.);
	const __m256d c  = _mm256_set1_pd( 231./ 16.);
	const __m256d d  = _mm256_set1_pd(3003./128.);

	__m256d y = _mm256_cvtps_pd(
			_mm_rsqrt_ps(
				_mm256_cvtpd_ps(x)));

	__m256d y2  = _mm256_mul_pd(y, y);
	__m256d y3  = _mm256_mul_pd(y, y2);
	__m256d y4  = _mm256_mul_pd(y2, y2);

	__m256d z   = _mm256_mul_pd(y3, y4);
	__m256d h   = _mm256_fnmadd_pd(x,  y2, one); // c - a * b

	__m256d poly = _mm256_fmadd_pd(d, h, c); // d*h + c
	poly = _mm256_fmadd_pd(poly, h, b); // (d*h + c)*h + b
	poly = _mm256_fmadd_pd(poly, h, a); // ((d*h + c)*h + b)*h + a

	__m256d zh  = _mm256_mul_pd(z, h);

	__m256d z1  = _mm256_fmadd_pd(zh, poly, z);

	return z1;
}


__attribute__((noinline))
void nbody_m256d_quad(
	const double eps2_sd,
	const EPI epi[],
	const int ni,
	const ParticleSimulator::SPJQuadrupole spj[],
	const int nj,
	Force force[])
{
	auto dup = [](double x)->__m256d { return _mm256_set1_pd(x); };

	const __m256d eps2 = dup(eps2_sd);

	double trqs[nj];
	for(int j=0; j<nj; j++){
		trqs[j] = 0.5 * spj[j].quad.getTrace();
	}

	for(int i=0; i<ni; i+=4){
		__m256d xi = _mm256_loadu_pd(&epi[i+0].pos.x);
		__m256d yi = _mm256_loadu_pd(&epi[i+1].pos.x);
		__m256d zi = _mm256_loadu_pd(&epi[i+2].pos.x);
		__m256d mi = _mm256_loadu_pd(&epi[i+3].pos.x);
		transpose_4ymm_pd(xi, yi, zi, mi);

		__m256d ax, ay, az, pot;
		ax = ay = az = pot = _mm256_set1_pd(0);

		for(int j=0; j<nj; j+=1){
			__m256d dx = _mm256_sub_pd(xi, dup(spj[j+0].pos.x));
			__m256d dy = _mm256_sub_pd(yi, dup(spj[j+0].pos.y));
			__m256d dz = _mm256_sub_pd(zi, dup(spj[j+0].pos.z));
			__m256d mj = _mm256_set1_pd(spj[j+0].mass);

			__m256d r2  = 
				_mm256_fmadd_pd(dz, dz, 
					_mm256_fmadd_pd(dy, dy,
						_mm256_fmadd_pd(dx, dx, eps2)));

		       __m256d ri7 = rsqrt_pow7_x5(r2);
		       __m256d ri5 = r2 * ri7;
		       __m256d ri3 = (r2 * r2) * ri7;

		       auto &qj = spj[j].quad;
		       __m256d qxx = dup(qj.xx), qyy = dup(qj.yy), qzz = dup(qj.zz);
		       __m256d qxy = dup(qj.xy), qyz = dup(qj.yz), qxz = dup(qj.xz);
		       __m256d trq = dup(trqs[j]); // 0.5*(qxx+qyy+qzz)

		       __m256d qr_x = qxx*dx + qxy*dy + qxz*dz;
		       __m256d qr_y = qxy*dx + qyy*dy + qyz*dz;
		       __m256d qr_z = qxz*dx + qyz*dy + qzz*dz;

		       __m256d qrr = qr_x*dx + qr_y*dy + qr_z*dz;

		       __m256d qrr_r7 = ri7 * qrr;
		       __m256d qrr_r5 = r2  * qrr_r7;
		       // ri5 *= dup(1.5);

		       __m256d B = dup(3.0) * ri5;
		       __m256d A = mj*ri3 - trq*B + dup(7.5)*qrr_r7;

#if 1
		       ax -= A*dx - B*qr_x;
		       ay -= A*dy - B*qr_y;
		       az -= A*dz - B*qr_z;
		       pot -= (mj*(r2*r2))*ri5 - trq*ri3 + dup(1.5)*qrr_r5;
#else
		       __m256d mri3 = mj * r2 * ri5;
		       ax  = _mm256_fnmadd_pd(mri3, dx, ax);
		       ay  = _mm256_fnmadd_pd(mri3, dy, ay);
		       az  = _mm256_fnmadd_pd(mri3, dz, az);
		       pot = _mm256_fnmadd_pd(mri3, r2, pot);
#endif
		}
		transpose_4ymm_pd(ax, ay, az, pot);

		auto accumulate = [](double *dst, __m256d val){
			__m256d sum = _mm256_loadu_pd(dst);
			sum = _mm256_add_pd(sum, val);
			_mm256_storeu_pd(dst, sum);
		};
		if(i+0 < ni) accumulate(&force[i+0].acc.x, ax);
		if(i+1 < ni) accumulate(&force[i+1].acc.x, ay);
		if(i+2 < ni) accumulate(&force[i+2].acc.x, az);
		if(i+3 < ni) accumulate(&force[i+3].acc.x, pot);
	}
}

#if 0
int main(){
	__attribute__((aligned(64)))
	double x[32], y[32];

	srand48(20220524);
	for(int i=0; i<32; i++){
		x[i] = drand48() * drand48();
	}
	for(int i=0; i<32; i+=4){
		*(__m256d *)&y[i] = rsqrt_pow7_x5(*(__m256d *)&x[i]);
	}

	for(int i=0; i<32; i++){
		double ref = pow(x[i], -3.5);
		double err = (y[i]-ref)/ref;
		printf("%a\t%a\t%e\n", ref, y[i], err);
	}
}
#endif
