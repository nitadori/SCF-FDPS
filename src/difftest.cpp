#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>
#include <ctime>

#include <algorithm>

#include <x86intrin.h>
//#include <omp.h>
//omp_get_wtime();

#include "scfdhmpi.h"

void accp_CB_mod0(struct Params *param,int ni,int ne,double *mass,
		double *x,double *y,double *z,double *ax,double *ay,
		double *az,double *pot,double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1])
{
	double sinmphi [LLMAX+1],
	       cosmphi [LLMAX+1];

	int nmax = param->nmax;
	int lmax = param->lmax;
	int lmin = param->lmin;
	int lskip = param->lskip;
	double G = param->G;
	double abasis = param->abasis;
	double abasis2 = param->abasis2;

	double scsum    [lmax+1][lmax+1][2][nmax+1]; // compact buffer, GNU ext.
	for(int l=lmin; l<=lmax; l+=lskip)
		for(int m=0; m<=l; m++)
			for(int n=0; n<=nmax; n++){
				scsum[l][m][0][n] = sinsum[n][l][m];
				scsum[l][m][1][n] = cossum[n][l][m];
			}
	double ultrasp  [lmax+1][nmax+1];  // transposed
	double ultrasp1 [lmax+1][nmax+1];  // transposed
	double lfact    [lmax+1];
	double plm      [lmax+1][lmax+1];
	double dplm     [lmax+1][lmax+1];

	for(int k=ni; k<ne; k++){
		double xi;
		double rho2   = x[k]*x[k] + y[k]*y[k];
		double rhoinv = 1.0/sqrt(rho2);
		double rho    = rho2 * rhoinv;
		double r2     = rho2 + z[k]*z[k];
		double rinv   = 1.0 / sqrt(r2);
		double r      = r2 * rinv;
		double costh  = z[k] * rinv;
		double sinth  = rho * rinv;

		cosmphi[0]=1.0;
		sinmphi[0]=0.0;
		double cosphi = cosmphi[1] = x[k] * rhoinv;
		double sinphi = sinmphi[1] = y[k] * rhoinv;

		for(int m=2; m<=lmax; m++){
			cosmphi[m]=cosphi*cosmphi[m-1] - sinphi*sinmphi[m-1];
			sinmphi[m]=cosphi*sinmphi[m-1] + sinphi*cosmphi[m-1];
		}

		{
#if 1
			xi = (r2-abasis2) / (r2+abasis2);
			xi = (r2-abasis2) / (r2+abasis2);
			double a2r2 = abasis2 + r2;
			double val  = abasis / sqrt(a2r2);
			double fac = abasis*r / a2r2;
#else
			xi = (r-abasis) / (r+abasis);
			double ar  = abasis + r;
			double val = mass[k] / ar;
			double fac = abasis*r / (ar*ar);
			(void)abasis2;
#endif
			lfact[0] = val;
			val *= fac;
			lfact[1] =  val;

			for(int l=2; l<=lmax; l++){
				val *= fac;
				lfact[l] = val;
			}
		}

		double potk = 0.0;
		double ar   = 0.0;
		double ath  = 0.0;
		double aphi = 0.0;

		for(int l=0; l<=lmax; l++){
			ultrasp [l][0] = 1.0;
			ultrasp [l][1] = twoalpha[l]*xi;
			ultrasp1[l][0] = 0.0;
			ultrasp1[l][1] = 1.0;
			
			double un   = ultrasp[l][1];
			double unm1 = 1.0;

			for(int n=1; n<=nmax-1; n++){
				ultrasp[l][n+1] = (c1[n][l]*xi*un - c2[n][l]*unm1) *c3[n];
				unm1 = un;
				un   = ultrasp[l][n+1];
				ultrasp1[l][n+1] = ((twoalpha[l]+n)*unm1 - (n+1)*xi*un) / (twoalpha[l] *(1.0-xi*xi));
			}
		}

		// sinth=sqrt(1.0-costh*costh);
		double ratio=1.0;
		double neginv_sin2th = -(1.0 / (sinth*sinth));

		plm [0][0] = 1.0;
		dplm[0][0] = 0.0;
		for(int m=0; m<=lmax; m++){
			if(m > 0){
				ratio *= -sinth;
				plm[m][m] = dblfact[m]*ratio;

				dplm[m][m]=m*costh*plm[m][m] * neginv_sin2th;
			}
			double plm1m = plm[m][m];
			double plm2m = 0.0;

			for(int l=m+1; l<=lmax; l++){
				plm[l][m] = (costh*(2*l-1)*plm1m - (l+m-1)*plm2m) * lmm[l][m];
				plm2m     = plm1m;
				plm1m     = plm[l][m];

				dplm[l][m]=(l*costh*plm[l][m] - (l+m)*plm[l-1][m]) * neginv_sin2th;
			}
		}

		for(int l=lmin; l<=lmax; l+=lskip){
			double temp3 = 0.0;
			double temp4 = 0.0;
			double temp5 = 0.0;
			double temp6 = 0.0;

			for(int m=0; m<=l; m++){
				double clm = 0.0;
				double dlm = 0.0;
				double elm = 0.0;
				double flm = 0.0;

				for(int n=0;n<=nmax;n++){
					clm += ultrasp [l][n] * scsum[l][m][1][n];
					dlm += ultrasp [l][n] * scsum[l][m][0][n];
					elm += ultrasp1[l][n] * scsum[l][m][1][n];
					flm += ultrasp1[l][n] * scsum[l][m][0][n];
				}

				temp3 +=   plm [l][m] * (clm*cosmphi[m] + dlm*sinmphi[m]);
				temp4 -=   plm [l][m] * (elm*cosmphi[m] + flm*sinmphi[m]);
				temp5 -=   dplm[l][m] * (clm*cosmphi[m] + dlm*sinmphi[m]);
				temp6 -= m*plm [l][m] * (dlm*cosmphi[m] - clm*sinmphi[m]);
			}

			double rnorm   = r/abasis;
			double rnorm2  = rnorm*rnorm;
			double rnorm21 = rnorm2+1.0;
			potk += lfact[l] *temp3;
			ar   += lfact[l] * (
				-temp3*(l/rnorm-(2*l+1)*rnorm/rnorm21)
				+temp4*4.0*rnorm*twoalpha[l]/(rnorm21*rnorm21) );
			ath  += lfact[l]*temp5;
			aphi += lfact[l]*temp6;
		}

#if 1
		ar    /=  abasis2;
		ath   *= -sinth/(abasis*r);
		aphi  /= (abasis*r*sinth);
		potk  /= abasis;
#else
		ar   /= abasis;
		ath  *= -sinth/r;
		aphi /= (r*sinth);
#endif
		ax [k] += G*(sinth*cosphi*ar + costh*cosphi*ath -sinphi*aphi);
		ay [k] += G*(sinth*sinphi*ar + costh*sinphi*ath +cosphi*aphi);
		az [k] += G*(costh*ar - sinth*ath);
		pot[k] += G*potk;
	}
}

#if 0
void test_rcp_rsqrt(){
	srand48(20220509);
	auto r = drand48; 
	double x[4] = {r()*r(), r()*r(), r()*r(), r()*r()}; 
	double rcp[4], rsq[4];
	double sum[2] = {-1.0, -1.0};

	*(__m256d *)rcp = rcp_full  (*(__m256d *)x);
	asm volatile ("nop");
	*(__m256d *)rsq = rsqrt_full(*(__m256d *)x);
	asm volatile ("nop");
	hadd_accum(sum, *(__m256d *)rcp, *(__m256d *)rsq);
	asm volatile ("nop");

	for(int i=0; i<4; i++){
		double h = 1.0 - x[i]*rcp[i];
		printf("%e %e %e\n", x[i], rcp[i], h);
	}

	for(int i=0; i<4; i++){
		double h = 1.0 - x[i]*rsq[i]*rsq[i];
		printf("%e %e %e\n", x[i], rsq[i], h);
	}

	printf("%24.16e %24.16e\n", sum[0], -1+rcp[0]+rcp[1]+rcp[2]+rcp[3]);
	printf("%24.16e %24.16e\n", sum[1], -1+rsq[0]+rsq[1]+rsq[2]+rsq[3]);
}
#endif

int main(int argc, char **argv)
{
	static double 
		mass[NMAX],
		x [NMAX], y [NMAX],z [NMAX],
		vx[NMAX], vy[NMAX],vz[NMAX] __attribute__((aligned(64)));
	static double ax[NMAX], ay[NMAX],az[NMAX], pot[NMAX];
	static double ax2[NMAX], ay2[NMAX],az2[NMAX], pot2[NMAX];
	struct Params param;
	double dblfact[LLMAX+1], twoalpha[LLMAX+1],
	       c1[LLMAX+1][LLMAX+1], c2[LLMAX+1][LLMAX+1],c3[NNMAX+1],
	       anltilde[NNMAX+1][LLMAX+1], lplusm[LLMAX+1][LLMAX+1], 
	       lmm[LLMAX+1][LLMAX+1], coeflm[LLMAX+1][LLMAX+1];

	MPI_Init(&argc,&argv);

	// test_rcp_rsqrt(); MPI_Finalize(); return 0;

	int me,n_pes;
	MPI_Comm_size(MPI_COMM_WORLD,&n_pes);
	MPI_Comm_rank(MPI_COMM_WORLD,&me);

	inparams(&param,me);

	int incont = param.incont;
	// int initsnap = param.initsnap;
	char *basisset = param.basisset;
	int nmax = param.nmax;
	int lmax = param.lmax;
	// int lmin = param.lmin;
	// int lskip = param.lskip;
	// int nsteps = param.nsteps;
	// int noutbod = param.noutbod;
	// int noutlog = param.noutlog;
	// double dtime = param.dtime;
	// double G = param.G;
	// double eps = param.eps;
	// double abasis = param.abasis;
	// double abasis2 = param.abasis2;
	// int zeroodd = param.zeroodd;
	// int zeroeven = param.zeroeven;
	// int fixacc = param.fixacc;
	// int ubodlog = param.ubodlog;

	int ndisk, nhalo;
	double tnow, etotal;
	inbods(incont, &ndisk, &nhalo, &tnow, mass, x, y, z, vx, vy, vz, &etotal, me, n_pes, 0);

	// initsys
	const double pi=4.0*atan(1.0);
	
	// accpot
	//const char *basisset = param.basisset;
	//int nmax = param.nmax;
	//int lmax = param.lmax;
	static double sinsum[NNMAX+1][LLMAX+1][LLMAX+1];
	static double cossum[NNMAX+1][LLMAX+1][LLMAX+1];
	static double sinsum2[NNMAX+1][LLMAX+1][LLMAX+1];
	static double cossum2[NNMAX+1][LLMAX+1][LLMAX+1];

	auto benchmark = [&](auto coef1, auto coef2, auto accp1, auto accp2, int nrep=1){
		for(int rep=0; rep<nrep; rep++){
			for(int n=0;n<=nmax;n++) for(int l=0;l<=lmax;l++) for(int m=0;m<=l;m++)
			{
				sinsum[n][l][m]=0.0;
				cossum[n][l][m]=0.0;
				sinsum2[n][l][m]=0.0;
				cossum2[n][l][m]=0.0;
			}

			for(int k=0; k<nhalo; k++){
				ax[k] = ay[k] = az[k] = pot[k] = 0.0;
				ax2[k] = ay2[k] = az2[k] = pot2[k] = 0.0;
			}
			double t0 = MPI_Wtime();

			coef1(&param,0,nhalo,sinsum,cossum,mass,x,y,z,
					c1,c2,c3,anltilde,dblfact,twoalpha,lplusm,lmm, coeflm);

			double t1 = MPI_Wtime();

			// scfcoef_CB_mod
			coef2(&param,0,nhalo,sinsum2,cossum2,mass,x,y,z,
					c1,c2,c3,anltilde,dblfact,twoalpha,lplusm,lmm, coeflm);

			double t2 = MPI_Wtime();

			accp1(&param,0,nhalo,mass,x,y,z,ax,ay,az,pot,dblfact,
					twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm, sinsum,cossum);

			double t3 = MPI_Wtime();

			accp2(&param,0,nhalo,mass,x,y,z,ax2,ay2,az2,pot2,dblfact,
					twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm, sinsum2,cossum2);

			double t4 = MPI_Wtime();

			printf("%e sec,  %e sec,  %e sec,  %e sec\n", t1-t0, t2-t1, t3-t2, t4-t3);

			int ncoef = (nmax+1)*(lmax+1)*(lmax+2)/2;
			double gops = 1.e-9 * nhalo * ncoef;
			printf("%e nsec, %e nsec, %e nsec, %e nsec\n", (t1-t0)/gops, (t2-t1)/gops, (t3-t2)/gops, (t4-t3)/gops);
			puts("");
			fflush(stdout);
		}
	};

	/*
	   if(strcmp(basisset,"CB") == 0){
	   puts("CB");
	// benchmark(scfcoef_CB, scfcoef_CB_mod, accp_CB, accp_CB_mod, 3);
	benchmark(scfcoef_CB, scfcoef_CB_avx, accp_CB, accp_CB_avx, 3);
	}else if(strcmp(basisset,"LH")==0){
	puts("LH");
	// benchmark(scfcoef_LH, scfcoef_LH_mod, accp_LH, accp_LH_mod, 3);
	benchmark(scfcoef_LH, scfcoef_LH_avx, accp_LH, accp_LH_avx, 3);
	}*/

	auto verify_coef = [&](){
		auto max = [](auto a, auto b) { return a>b ? a : b; };
		double err_max = 0.0;
		for(int n=0;n<=nmax;n++) for(int l=0;l<=lmax;l++) for(int m=0;m<=l;m++)
		{
			double c = cossum[n][l][m];
			double s = sinsum[n][l][m];
			double c2 = cossum2[n][l][m];
			double s2 = sinsum2[n][l][m];

			double dc = c2 - c;
			double ds = s2 - s;
			double err = sqrt( (dc*dc + ds*ds) / (c*c + s*s) );
			err_max = max(err, err_max);

#if 0
			printf("(%2d, %2d, %2d) : (%+e, %+e) (%+e, %+e) %e\n", 
					n, l, m, c, s, c2, s2, err);
#endif
		}
		printf("err_max : %e\n", err_max);
		fflush(stdout);
	};

	auto verify_accp = [&](){
		static double aerr[NMAX], perr[NMAX];
		for(int k=0; k<nhalo; k++){
			double x = ax[k];
			double y = ay[k];
			double z = az[k];
			double p = pot[k];
			double dx = ax2[k] - x;
			double dy = ay2[k] - y;
			double dz = az2[k] - z;
			double dp = pot2[k] - p;

			aerr[k] = sqrt( (dx*dx + dy*dy + dz*dz) / (x*x + y*y + z*z) );
			perr[k] = fabs(dp / p);
		}
		auto gt = [](auto a, auto b){ return a>b; };
		std::sort(aerr, aerr+nhalo, gt);
		std::sort(perr, perr+nhalo, gt);
		for(int k=0; k<5; k++){
			printf("%2d Δa/a : %e ΔΦ/Φ: %e\n", k, aerr[k], perr[k]);
		}
		puts("");
		fflush(stdout);
	};

	(void)basisset;
	getcoef_CB(nmax,lmax,pi,dblfact,anltilde,coeflm,twoalpha, lmm,lplusm,c1,c2,c3);
	puts("CB_mod");
	benchmark(scfcoef_CB, scfcoef_CB_mod, accp_CB, accp_CB_mod, 3);
	verify_coef(); verify_accp();
	puts("CB_avx");
	benchmark(scfcoef_CB, scfcoef_CB_avx, accp_CB, accp_CB_avx, 3);
	verify_coef();
	verify_accp();

	getcoef_LH(nmax,lmax,pi,dblfact,anltilde,coeflm,twoalpha, lmm,lplusm,c1,c2,c3);
	puts("LH_mod");
	benchmark(scfcoef_LH, scfcoef_LH_mod, accp_LH, accp_LH_mod, 3);
	verify_coef(); verify_accp();
	puts("LH_avx");
	benchmark(scfcoef_LH, scfcoef_LH_avx, accp_LH, accp_LH_avx, 3);
	verify_coef(); verify_accp();
	
	MPI_Finalize();

	return 0;
}
