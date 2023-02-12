#define CONN3_IMPL(x, y, z) x ## _ ## y ## _ ## z
#define CONN3(x, y, z) CONN3_IMPL(x, y, z)

void CONN3(scfcoef, SCF_BASE_NAME, mod)(
		struct Params *param,int ni,int ne,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double c1[][LLMAX+1],double c2[][LLMAX+1],
		double *c3,double anltilde[][LLMAX+1],double *dblfact,
		double *twoalpha,double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1]) 
{
	int nmax = param->nmax;
	int lmax = param->lmax;
	int lmin = param->lmin;
	int lskip = param->lskip;
	double abasis = param->abasis;
	double abasis2 = param->abasis2;

	double cosmphi  [lmax+1]; 
	double sinmphi  [lmax+1];
	// double lfact    [lmax+1];
	double plm      [lmax+1][lmax+1];

	double ultraspt [lmax+1][nmax+1];  // transposed
	double scsum    [lmax+1][lmax+1][2][nmax+1]; // compact buffer, GNU ext.

	for(int l=lmin; l<=lmax; l+=lskip){
		for(int m=0; m<=l; m++){
			for(int n=0; n<=nmax; n++){
				scsum[l][m][0][n]  = 0.0;
				scsum[l][m][1][n]  = 0.0;
			}
		}
	}

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

#if 1
		for(int m=2; m<=lmax; m++){
			cosmphi[m]=cosphi*cosmphi[m-1] - sinphi*sinmphi[m-1];
			sinmphi[m]=cosphi*sinmphi[m-1] + sinphi*cosmphi[m-1];
		}
#else
		double cos2phi = cosmphi[2] = cosphi*cosphi - sinphi*sinphi;
		double sin2phi = sinmphi[2] = 2.*cosphi*sinphi;
		for(int m=3; m<=lmax; m+=2){
			cosmphi[m]   = cos2phi*cosmphi[m-2] - sin2phi*sinmphi[m-2];
			sinmphi[m]   = cos2phi*sinmphi[m-2] + sin2phi*cosmphi[m-2];
			cosmphi[m+1] = cos2phi*cosmphi[m-1] - sin2phi*sinmphi[m-1];
			sinmphi[m+1] = cos2phi*sinmphi[m-1] + sin2phi*cosmphi[m-1];
		}
#endif

#if 0
		for(int m=2; m<=lmax; m++){
#endif

#if   defined SCF_BASE_NAME_is_CB
		xi = (r2-abasis2) / (r2+abasis2);
		double a2r2 = abasis2 + r2;
		double lc_val = mass[k] * abasis / sqrt(a2r2);
		double lc_fac = abasis*r / a2r2;
#elif defined SCF_BASE_NAME_is_LH
		xi = (r-abasis) / (r+abasis);
		double ar  = abasis + r;
		double lc_val = mass[k] / ar;
		double lc_fac = abasis*r / (ar*ar);
		(void)abasis2;
#endif

		for(int l=0; l<=lmax; l++){
			// double ultrasp_0l = 1.0;
			double ultrasp_1l = twoalpha[l]*xi;
			double un = ultrasp_1l * lc_val;
			double unm1 = lc_val;

			// lfact[l] = val;
			lc_val *= lc_fac;

			ultraspt[l][0] = unm1 * anltilde[0][l];
			ultraspt[l][1] = un   * anltilde[1][l];

			for(int n=1; n<=nmax-1; n++){
				double ultrasp_np1l = (c1[n][l]*xi*un - c2[n][l]*unm1) * c3[n];
				unm1 = un;
				un   = ultrasp_np1l;

				ultraspt[l][n+1] = un * anltilde[n+1][l];
			}
		}

		double ratio=1.0;
		for(int m=0; m<=lmax; m++){
			double plm_mm = 1.0;
			if(m > 0){
				ratio *= -sinth;
				plm_mm = dblfact[m]*ratio;
			}
			double plm1m = plm_mm;
			double plm2m = 0.0;
			// plm[m][m] = plm_mm  * coeflm[m][m] * lfact[m];
			plm[m][m] = plm_mm;//  * lfact[m];

			for(int l=m+1; l<=lmax; l++){
				double plm_lm = (costh*(2*l-1)*plm1m - lplusm[l][m]*plm2m) * lmm[l][m];
				plm2m = plm1m;
				plm1m = plm_lm;
				// plm[l][m] = plm_lm  * coeflm[l][m] * lfact[l];
				plm[l][m] = plm_lm;// * lfact[l];
			}
		}

		assert(1 == lskip);
		for(int l=lmin; l<=lmax; l+=lskip){
			for(int m=0; m<=l; m++){
				for(int n=0; n<=nmax; n++){
					//double ttemp5 = temp5*plm[l][m]*coeflm[l][m];
					double unl = ultraspt[l][n];
					double templm = unl * plm[l][m];
					double stemp  = templm * sinmphi[m];
					double ctemp  = templm * cosmphi[m];
					// sinsum[n][l][m] += stemp;
					// cossum[n][l][m] += ctemp;
					scsum[l][m][0][n] += stemp;
					scsum[l][m][1][n] += ctemp;
				}
			}
		}
	}

	scfmpi_sum_double(scsum[0][0][0], sizeof(scsum)/sizeof(double));

	for(int l=lmin; l<=lmax; l+=lskip)
		for(int m=0; m<=l; m++)
			for(int n=0; n<=nmax; n++){
				sinsum[n][l][m] = scsum[l][m][0][n] * coeflm[l][m];
				cossum[n][l][m] = scsum[l][m][1][n] * coeflm[l][m];
			}
#if 0
	icount=0;
	for(l=lmin;l<=lmax;l += lskip)
		for(m=0;m<=l;m++)
			for(n=0;n<=nmax;n++){
				icount += 2;	
				tmpsumsnd[icount-2]=sinsum[n][l][m];
				tmpsumsnd[icount-1]=cossum[n][l][m];
			}
	MPI_Allreduce(tmpsumsnd,tmpsumrcv,icount,MPI_DOUBLE,MPI_SUM,
			MPI_COMM_WORLD);
	icount=0;
	for(l=lmin;l<=lmax;l += lskip)
		for(m=0;m<=l;m++)
			for(n=0;n<=nmax;n++){
				icount += 2;	
				sinsum[n][l][m]=tmpsumrcv[icount-2];
				cossum[n][l][m]=tmpsumrcv[icount-1];
			}
#endif
	return;
}

void CONN3(scfcoef, SCF_BASE_NAME, avx)(
		struct Params *param,int ni,int ne,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double c1[][LLMAX+1],double c2[][LLMAX+1],
		double *c3,double anltilde[][LLMAX+1],double *dblfact,
		double *twoalpha,double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1]) 
{
#if 0
	using Ftype = double;
	const int VLEN = 1;
	auto fadd = [](Ftype a, Ftype b)->Ftype { return a + b; };
	auto fsub = [](Ftype a, Ftype b)->Ftype { return a - b; };
	auto fmul = [](Ftype a, Ftype b)->Ftype { return a * b; };
	auto fmul3  = [](Ftype a, Ftype b, Ftype c)->Ftype { return a * b * c; };
	auto finv   = [](Ftype x)->Ftype { return 1.0/x; };
	auto frsqrt = [](Ftype x)->Ftype { return 1.0/sqrt(x); };
	auto fneg   = [](Ftype x)->Ftype { return -x; };
	auto dup    = [](double x)->Ftype { return x; };
	auto load   = [](const double *x)->Ftype { return *x; };
	auto accumulate = [](double *sc, Ftype s, Ftype c){
		sc[0] += s; sc[1] += c;
	};
#else
	using Ftype = __m256d;
	const int VLEN = 4;
	auto fadd = [](Ftype a, Ftype b)->Ftype { return _mm256_add_pd(a, b); };
	auto fsub = [](Ftype a, Ftype b)->Ftype { return _mm256_sub_pd(a, b); };
	auto fmul = [](Ftype a, Ftype b)->Ftype { return _mm256_mul_pd(a, b); };
	auto fmul3  = [](Ftype a, Ftype b, Ftype c)->Ftype { 
		return _mm256_mul_pd(_mm256_mul_pd(a, b), c); };
	auto finv   = [](Ftype x)->Ftype { return rcp_full(x); };
	auto frsqrt = [](Ftype x)->Ftype { return rsqrt_full(x); };
	auto fneg   = [](Ftype x)->Ftype { return -x; };
	auto dup    = [](double x)->Ftype { return _mm256_set1_pd(x); };
	auto load   = [](const double *x)->Ftype { return _mm256_load_pd(x); };
	auto accumulate = [](double *sc, Ftype s, Ftype c){
		hadd_accum(sc, s, c);
	};
#endif
	if(0 == ni){
		// nhalo must be multiple of VLEN
		assert(0 == ne%VLEN);
	}else{
		assert(0 == ni%VLEN); // for safety
		// clean-up for ndisk
		auto end = align_up<VLEN>(ne);
		for(auto i=ne; i<end; i++){
			// need some value for coordinates, for safety
			x[i] = 1.0;
			y[i] = 2.0;
			z[i] = 3.0;
			mass[i] = 0.0;
		}
	}

	int nmax = param->nmax;
	int lmax = param->lmax;
	int lmin = param->lmin;
	int lskip = param->lskip;
	double abasis = param->abasis;
	double abasis2 = param->abasis2;

	Ftype cosmphi [lmax+1]; 
	Ftype sinmphi [lmax+1];
	// Ftype lfact   [lmax+1];
	Ftype plm     [lmax+1][lmax+1];

	Ftype ultraspt[lmax+1][nmax+1];
	int nsize = align_up<2>(nmax+1);
	double scsum  [lmax+1][lmax+1][nsize][2] __attribute__((aligned(64))); // compact buffer, GNU ext.
	(void)accumulate;

	for(int l=lmin; l<=lmax; l+=lskip){
		for(int m=0; m<=l; m++){
			for(int n=0; n<=nmax; n++){
				scsum[l][m][n][0]  = 0.0;
				scsum[l][m][n][1]  = 0.0;
			}
		}
	}

	for(int k=ni; k<ne; k+=VLEN){
		Ftype xk = load(x+k);
		Ftype yk = load(y+k);
		Ftype zk = load(z+k);

		Ftype rho2 = fadd(fmul(xk, xk), fmul(yk, yk));
		Ftype rhoinv = frsqrt(rho2);
		Ftype rho = fmul(rho2, rhoinv);

		Ftype r2   = fadd(rho2, fmul(zk, zk));
		Ftype rinv = frsqrt(r2);
		Ftype r    = fmul(r2, rinv);

		Ftype costh = fmul(zk, rinv);
		Ftype sinth = fmul(rho, rinv);

		Ftype a2 = dup(abasis2);
		Ftype a  = dup (abasis);
		Ftype xi;

		cosmphi[0] = dup(1.0);
		sinmphi[0] = dup(0.0);

		Ftype cosphi = cosmphi[1] = fmul(xk, rhoinv);
		Ftype sinphi = sinmphi[1] = fmul(yk, rhoinv);

		for(int m=2; m<=lmax; m++){
			cosmphi[m] = fsub(fmul(cosphi, cosmphi[m-1]), fmul(sinphi, sinmphi[m-1]));
			sinmphi[m] = fadd(fmul(cosphi, sinmphi[m-1]), fmul(sinphi, cosmphi[m-1]));
		}

		Ftype mk = load(mass+k);
#if   defined SCF_BASE_NAME_is_CB
		xi = fmul(fsub(r2, a2), finv(fadd(r2, a2)));
		Ftype a2r2_rsq = frsqrt(fadd(a2, r2));
		Ftype a1 = dup(abasis);
		Ftype lc_val = fmul3(mk, a1, a2r2_rsq);
		Ftype lc_fac = fmul3(a1, r, fmul(a2r2_rsq, a2r2_rsq));
		(void)a;
#elif defined SCF_BASE_NAME_is_LH
		xi = fmul(fsub(r, a), finv(fadd(r, a)));
		Ftype ar    = fadd(a, r);
		Ftype arinv = finv(ar);
		Ftype lc_val = fmul(mk, arinv);
		Ftype lc_fac = fmul(fmul(a, r), fmul(arinv, arinv));
		(void)abasis2;
		(void)a2;
#endif

		for(int l=0; l<=lmax; l++){
			Ftype ultrasp_1l = fmul(dup(twoalpha[l]), xi);
			Ftype un = fmul(lc_val, ultrasp_1l);
			Ftype unm1 = lc_val;

			lc_val = fmul(lc_val, lc_fac); // l-coeff value

			ultraspt[l][0] = fmul(unm1, dup(anltilde[0][l]));
			ultraspt[l][1] = fmul(un  , dup(anltilde[1][l]));

			for(int n=1; n<=nmax-1; n++){
				Ftype tmp1 = fmul3(dup(c1[n][l]), xi, un);
				Ftype tmp2 = fmul(dup(c2[n][l]), unm1);
				Ftype  ultrasp_np1l = fmul(dup(c3[n]), fsub(tmp1, tmp2));
				unm1 = un;
				un   = ultrasp_np1l;

				ultraspt[l][n+1] = fmul(un, dup(anltilde[n+1][l]));
			}
		}

		Ftype ratio = dup(1.0);
		for(int m=0; m<=lmax; m++){
			Ftype plm_mm = dup(1.0);
			if(m > 0){
				ratio = fmul(ratio, fneg(sinth));
				plm_mm = fmul(dup(dblfact[m]), ratio);
			}
			Ftype plm1m = plm_mm;
			Ftype plm2m = dup(0.0);
			plm[m][m] = plm_mm; // fmul(plm_mm, lfact[m]);

			for(int l=m+1; l<=lmax; l++){
				Ftype templ = fmul3(costh, plm1m, dup((double)(2*l-1)));
				Ftype tempr = fmul(plm2m, dup((double)(l+m-1))); // lplusm[k][m]=l+m-1.0
				Ftype plm_lm = fmul(fsub(templ, tempr), dup(lmm[l][m]));

				plm2m = plm1m;
				plm1m = plm_lm;
				plm[l][m] = plm_lm; // fmul(plm_lm, lfact[l]);
			}
		}

		// assert(1 == lskip);
		for(int l=lmin; l<=lmax; l+=lskip){
			for(int m=0; m<=l; m++){
				for(int n=0; n<nmax+1; n+=2){
					Ftype unl = ultraspt[l][n];
					Ftype templm = fmul(unl, plm[l][m]);
					Ftype stemp  = fmul(templm, sinmphi[m]);
					Ftype ctemp  = fmul(templm, cosmphi[m]);
					// accumulate(scsum[l][m][n], stemp, ctemp);

					Ftype unl1 = ultraspt[l][n+1];
					Ftype templm1 = fmul(unl1, plm[l][m]);
					Ftype stemp1  = fmul(templm1, sinmphi[m]);
					Ftype ctemp1  = fmul(templm1, cosmphi[m]);
					// accumulate(scsum[l][m][n+1], stemp1, ctemp1);

					hadd_accum2(scsum[l][m][n], stemp, ctemp, stemp1, ctemp1);
				}
				#if 0
				for(int n=0; n<nmax+1; n++){
					double s = scsum[l][m][n][0];
					double c = scsum[l][m][n][1];
					if(ext_isnan(s+c)){
						int r = scfmpi_rank();
						fprintf(stderr, "%d: NaN/inf at (%d, %d, %d), [%e, %e]\n",
								r, n, l, m,  s, c);
						int k4 = k/4*4;
						for(int kk=k4; kk<k4+4; kk++){
							fprintf(stderr, "k=%d, posm: (%e, %e, %e | %e)\n",
									kk, x[kk], y[kk], z[kk], mass[kk]);
						}
						scfmpi_abort(0);
					}
				}
				#endif
			}
		}
	}
	
	// MPI_Allreduce(MPI_IN_PLACE, scsum, sizeof(scsum)/sizeof(double), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	scfmpi_sum_double(scsum[0][0][0], sizeof(scsum)/sizeof(double));

	// int nnan = 0;
	for(int l=lmin; l<=lmax; l+=lskip){
		for(int m=0; m<=l; m++){
			for(int n=0; n<=nmax; n++){
				sinsum[n][l][m] = scsum[l][m][n][0] * coeflm[l][m];
				cossum[n][l][m] = scsum[l][m][n][1] * coeflm[l][m];

				#if 0
				double num = sinsum[n][l][m] + cossum[n][l][m];
				if( ext_isnan(num) || !ext_isfinite(num) ){
					int r = scfmpi_rank();
					fprintf(stderr, "%d: NaN/inf at (%d, %d, %d), [%e, %e]\n",
							r, n, l, m,  sinsum[n][l][m], cossum[n][l][m]);
					if(++nnan >= 10) exit(0);
				}
				#endif
			}
		}
	}
}

void CONN3(accp, SCF_BASE_NAME, mod)(
		struct Params *param,int ni,int ne,double *mass,
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

#if   defined SCF_BASE_NAME_is_CB
		xi = (r2-abasis2) / (r2+abasis2);
		double a2r2 = abasis2 + r2;
		double lc_val  = G * abasis / sqrt(a2r2);
		double lc_fac = abasis*r / a2r2;
#elif defined SCF_BASE_NAME_is_LH
		xi = (r-abasis) / (r+abasis);
		double apr = abasis + r;
		double lc_val = G * abasis / apr;
		double lc_fac = abasis*r / (apr*apr);
		(void)abasis2;
#endif

		double potk = 0.0;
		double ar   = 0.0;
		double ath  = 0.0;
		double aphi = 0.0;

		for(int l=0; l<=lmax; l++){
			lfact[l] = lc_val;
			lc_val *= lc_fac;

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

#if   defined SCF_BASE_NAME_is_CB
			double rnorm   = r/abasis;
			double rnorm2  = rnorm*rnorm;
			double rnorm21 = rnorm2+1.0;
			potk += lfact[l] *temp3;
			ar   += lfact[l] * (
				-temp3*(l/rnorm-(2*l+1)*rnorm/rnorm21)
				+temp4*4.0*rnorm*twoalpha[l]/(rnorm21*rnorm21) );
			ath  += lfact[l]*temp5;
			aphi += lfact[l]*temp6;
#elif defined SCF_BASE_NAME_is_LH
			double rnorm  = r/abasis;
			double rnorm1 = rnorm+1.0;
			potk += lfact[l] * temp3;
			ar   += lfact[l] * (-temp3*(l/rnorm-(2*l+1)/rnorm1) +temp4*(8*l+6)/(rnorm1*rnorm1));
			ath  += lfact[l] * temp5;
			aphi += lfact[l] * temp6;
#endif
		}

#if   defined SCF_BASE_NAME_is_CB
		ar    /=  abasis2;
		ath   *= -sinth/(abasis*r);
		aphi  /= (abasis*r*sinth);
		potk  /= abasis;
#elif defined SCF_BASE_NAME_is_LH
		ar   /= abasis;
		ath  *= -sinth/r;
		aphi /= (r*sinth);
#endif
		double arho = sinth*ar + costh*ath;
		ax [k] += (cosphi*arho - sinphi*aphi);
		ay [k] += (sinphi*arho + cosphi*aphi);
		az [k] += (costh*ar - sinth*ath);
		pot[k] += potk;
	}
}

void CONN3(accp, SCF_BASE_NAME, avx)(
		struct Params *param,int ni,int ne,double *mass,
		double *x,double *y,double *z,double *ax,double *ay,
		double *az,double *pot,double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1])
{
#if 0
	using Ftype = double;
	const int VLEN = 1;
	auto fadd = [](Ftype a, Ftype b)->Ftype { return a + b; };
	auto fsub = [](Ftype a, Ftype b)->Ftype { return a - b; };
	auto fmul = [](Ftype a, Ftype b)->Ftype { return a * b; };
	auto fmul3  = [](Ftype a, Ftype b, Ftype c)->Ftype { return a * b * c; };
	auto fmadd  = [](Ftype a, Ftype b, Ftype c)->Ftype { return c + a * b; };
	auto fmsub  = [](Ftype a, Ftype b, Ftype c)->Ftype { return c - a * b; };
	auto finv   = [](Ftype x)->Ftype { return 1.0/x; };
	auto frsqrt = [](Ftype x)->Ftype { return 1.0/sqrt(x); };
	auto fneg   = [](Ftype x)->Ftype { return -x; };
	auto dup    = [](double x)->Ftype { return x; };
	auto load   = [](const double *x)->Ftype { return *x; };
	auto store  = [](double *x, double v) { *x = v; };
	auto fdpad  = [](Ftype a, Ftype b, Ftype c, Ftype d)->Ftype { return a*b + c*d; };
	auto fdpsb  = [](Ftype a, Ftype b, Ftype c, Ftype d)->Ftype { return a*b - c*d; };
	// auto accumulate = [](double *sc, Ftype s, Ftype c){ sc[0] += s; sc[1] += c; };
#else
	using Ftype = __m256d;
	const int VLEN = 4;
	auto fadd = [](Ftype a, Ftype b)->Ftype { return _mm256_add_pd(a, b); };
	auto fsub = [](Ftype a, Ftype b)->Ftype { return _mm256_sub_pd(a, b); };
	auto fmul = [](Ftype a, Ftype b)->Ftype { return _mm256_mul_pd(a, b); };
	auto fmul3  = [](Ftype a, Ftype b, Ftype c)->Ftype { 
		return _mm256_mul_pd(_mm256_mul_pd(a, b), c); };
	auto fmadd = [](Ftype a, Ftype b, Ftype c)->Ftype { return _mm256_fmadd_pd(a, b, c); };
	auto fmsub = [](Ftype a, Ftype b, Ftype c)->Ftype { return _mm256_fnmadd_pd(a, b, c); };
	auto fdpad = [](Ftype a, Ftype b, Ftype c, Ftype d)->Ftype { 
		return _mm256_fmadd_pd(c, d, _mm256_mul_pd(a, b)); };
	auto fdpsb = [](Ftype a, Ftype b, Ftype c, Ftype d)->Ftype { 
		return _mm256_fnmadd_pd(c, d, _mm256_mul_pd(a, b)); };
	auto finv   = [](Ftype x)->Ftype { return rcp_full(x); };
	auto frsqrt = [](Ftype x)->Ftype { return rsqrt_full(x); };
	auto fneg   = [](Ftype x)->Ftype { return -x; };
	auto dup    = [](double x)->Ftype { return _mm256_set1_pd(x); };
	auto load   = [](const double *x)->Ftype { return _mm256_load_pd(x); };
	auto store  = [](double *p, Ftype v)->void { _mm256_store_pd(p, v); };
	// auto accumulate = [](double *sc, Ftype s, Ftype c){ hadd_accum(sc, s, c); };
#endif
	if(0 == ni){
		// nhalo must be multiple of VLEN
		assert(0 == ne%VLEN);
	}else{
		assert(0 == ni%VLEN); // for safety
	}

	Ftype sinmphi [LLMAX+1],
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
	Ftype ultrasp  [lmax+1][nmax+1];  // transposed
	Ftype ultrasp1 [lmax+1][nmax+1];  // transposed
	Ftype lfact    [lmax+1];
	Ftype plm      [lmax+1][lmax+1];
	Ftype dplm     [lmax+1][lmax+1];

	for(int k=ni; k<ne; k+=VLEN){
		Ftype xk = load(x+k);
		Ftype yk = load(y+k);
		Ftype zk = load(z+k);

		Ftype rho2 = fadd(fmul(xk, xk), fmul(yk, yk));
		Ftype rhoinv = frsqrt(rho2);
		Ftype rho = fmul(rho2, rhoinv);

		Ftype r2   = fadd(rho2, fmul(zk, zk));
		Ftype rinv = frsqrt(r2);
		Ftype r    = fmul(r2, rinv);

		Ftype costh = fmul(zk, rinv);
		Ftype sinth = fmul(rho, rinv);

		Ftype a2 = dup(abasis2);
		Ftype a  = dup (abasis);
		Ftype xi;

		cosmphi[0] = dup(1.0);
		sinmphi[0] = dup(0.0);

		Ftype cosphi = cosmphi[1] = fmul(xk, rhoinv);
		Ftype sinphi = sinmphi[1] = fmul(yk, rhoinv);

		for(int m=2; m<=lmax; m++){
			cosmphi[m] = fsub(fmul(cosphi, cosmphi[m-1]), fmul(sinphi, sinmphi[m-1]));
			sinmphi[m] = fadd(fmul(cosphi, sinmphi[m-1]), fmul(sinphi, cosmphi[m-1]));
		}

#if   defined SCF_BASE_NAME_is_CB
		xi = fmul(fsub(r2, a2), finv(fadd(r2, a2)));
		Ftype a2r2_rsq = frsqrt(fadd(a2, r2));
		Ftype a1 = a;
		Ftype lc_val = fmul3(dup(G), a1, a2r2_rsq);
		Ftype lc_fac = fmul3(a1, r, fmul(a2r2_rsq, a2r2_rsq));
		(void)a;
#elif defined SCF_BASE_NAME_is_LH
		xi = fmul(fsub(r, a), finv(fadd(r, a)));
		Ftype apr   = fadd(a, r);
		Ftype arinv = finv(apr);
		Ftype lc_val = fmul3(dup(G), a, arinv);
		Ftype lc_fac = fmul(fmul(a, r), fmul(arinv, arinv));
		(void)abasis2;
		(void)a2;
#endif

		Ftype potk = dup(0.0);
		Ftype ar   = dup(0.0);
		Ftype ath  = dup(0.0);
		Ftype aphi = dup(0.0);

		for(int l=0; l<=lmax; l++){
			lfact[l] = lc_val;
			lc_val = fmul(lc_val, lc_fac);

			ultrasp [l][0] = dup(1.0);
			ultrasp [l][1] = fmul(dup(twoalpha[l]), xi);
			ultrasp1[l][0] = dup(0.0);
			ultrasp1[l][1] = dup(1.0);
			
			Ftype un   = ultrasp[l][1];
			Ftype unm1 = dup(1.0);

			for(int n=1; n<=nmax-1; n++){
				//ultrasp[l][n+1] = (c1[n][l]*xi*un - c2[n][l]*unm1) *c3[n];
				Ftype tmp1 = fmul3(dup(c1[n][l]), xi, un);
				Ftype tmp2 = fmul(dup(c2[n][l]), unm1);
				ultrasp[l][n+1] = fmul(dup(c3[n]), fsub(tmp1, tmp2));

				unm1 = un;
				un   = ultrasp[l][n+1];
				// ultrasp1[l][n+1] = ((twoalpha[l]+n)*unm1 - (n+1)*xi*un) / (twoalpha[l] * (1.0-xi*xi));
				Ftype tmp3 = fsub(fmul(dup(twoalpha[l]+n), unm1), fmul3(dup(n+1), xi, un));
				Ftype tmp4 = fmul(dup(twoalpha[l]), fmsub(xi, xi, dup(1.0))); 
				ultrasp1[l][n+1] = fmul(tmp3, finv(tmp4));
			}
		}

		Ftype ratio = dup(1.0);
		Ftype neginv_sin2th = finv(fmul(fneg(sinth) , sinth));

		plm [0][0] = dup(1.0);
		dplm[0][0] = dup(0.0);
		for(int m=0; m<=lmax; m++){
			if(m > 0){
				ratio = fmul(ratio, fneg(sinth));
				plm[m][m] = fmul(dup(dblfact[m]), ratio);

				dplm[m][m] = fmul(fmul(dup(m), costh), fmul(plm[m][m], neginv_sin2th));
			}
			Ftype plm1m = plm[m][m];
			Ftype plm2m = dup(0.0);

			for(int l=m+1; l<=lmax; l++){
				// plm[l][m] = (costh*(2*l-1)*plm1m - (l+m-1)*plm2m) * lmm[l][m];
				plm[l][m] = fmul(fsub(fmul3(costh, dup(2*l-1), plm1m), fmul(dup(l+m-1), plm2m)) ,dup(lmm[l][m]));
				plm2m     = plm1m;
				plm1m     = plm[l][m];

				dplm[l][m]=fmul(fsub(fmul3(dup(l), costh, plm[l][m]),  fmul(dup(l+m), plm[l-1][m])), neginv_sin2th);
			}
		}

		for(int l=lmin; l<=lmax; l+=lskip){
			Ftype temp3 = dup(0.0);
			Ftype temp4 = dup(0.0);
			Ftype temp5 = dup(0.0);
			Ftype temp6 = dup(0.0);

			for(int m=0; m<=l; m++){
				Ftype clm = dup(0.0);
				Ftype dlm = dup(0.0);
				Ftype elm = dup(0.0);
				Ftype flm = dup(0.0);

				for(int n=0;n<=nmax;n++){
					Ftype s = dup(scsum[l][m][0][n]);
					Ftype c = dup(scsum[l][m][1][n]);
					Ftype u = ultrasp [l][n];
					Ftype v = ultrasp1[l][n];
					clm = fmadd(u, c, clm);
					dlm = fmadd(u, s, dlm);
					elm = fmadd(v, c, elm);
					flm = fmadd(v, s, flm);
				}

				temp3 = fmadd(plm [l][m], fdpad(clm, cosmphi[m], dlm, sinmphi[m]), temp3);
				temp4 = fmsub(plm [l][m], fdpad(elm, cosmphi[m], flm, sinmphi[m]), temp4);
				temp5 = fmsub(dplm[l][m], fdpad(clm, cosmphi[m], dlm, sinmphi[m]), temp5);
				Ftype mplm = dup(m) * plm[l][m];
				temp6 = fmsub(mplm, fdpsb(dlm, cosmphi[m], clm, sinmphi[m]), temp6);
			}

#if   defined SCF_BASE_NAME_is_CB
			Ftype rnorm   = fmul(r, dup(1.0/abasis));
			Ftype rnorm2  = fmul(rnorm, rnorm);
			Ftype rnorm21 = fadd(rnorm2, dup(1.0));
			Ftype invr21  = finv(rnorm21);
			/* ar   += lfact[l] * (
			 	-temp3*(l/rnorm - (2*l+1)*rnorm/rnorm21)
			 	+temp4*4.0*rnorm*twoalpha[l]/(rnorm21*rnorm21) ); */
			Ftype tt4 = fmul(fmul3(temp4, rnorm, dup(4.0*twoalpha[l])), fmul(invr21, invr21));
			Ftype tt3 = fmul(temp3, fsub(fmul(dup(l), finv(rnorm)), fmul3(dup(2*l+1), rnorm, invr21)));
#elif defined SCF_BASE_NAME_is_LH
			Ftype rnorm  = fmul(r, dup(1.0/abasis));
			Ftype rnorm1 = fadd(rnorm, dup(1.0));
			Ftype invr1  = finv(rnorm1);
			// ar   += lfact[l] * (-temp3*(l/rnorm-(2*l+1)/rnorm1) +temp4*(8*l+6)/(rnorm1*rnorm1));
			Ftype tt3 = fmul(temp3, fsub(fmul(dup(l), finv(rnorm)), fmul(dup(2*l+1), invr1)));
			Ftype tt4 = fmul(fmul(temp4, dup(8*l+6)), fmul(invr1, invr1));
#endif
			potk = fmadd(lfact[l], temp3, potk);
			ar   = fmadd(lfact[l], fsub(tt4, tt3), ar);
			ath  = fmadd(lfact[l], temp5, ath);
			aphi = fmadd(lfact[l], temp6, aphi);
		}

#if   defined SCF_BASE_NAME_is_CB
		ar   = fmul(ar, dup(1.0/abasis2));
		ath  = fmul3(ath, fneg(sinth), fmul(rinv, dup(1.0/abasis)));
		aphi = fmul(aphi, finv(fmul3(dup(abasis), r, sinth)));
		potk = fmul(potk, dup(1.0/abasis));
#elif defined SCF_BASE_NAME_is_LH
		ar   = fmul(ar, dup(1.0/abasis));
		ath  = fmul3(ath, fneg(sinth), rinv);
		aphi = fmul(aphi, finv(r*sinth));
#endif
		// ax [k] += G*(sinth*cosphi*ar + costh*cosphi*ath -sinphi*aphi);
		// ay [k] += G*(sinth*sinphi*ar + costh*sinphi*ath +cosphi*aphi);
		// az [k] += G*(costh*ar - sinth*ath);
		// pot[k] += G*potk;
		Ftype arho = fdpad(sinth, ar, costh, ath);
		Ftype tx = fdpsb(cosphi, arho, sinphi, aphi);
		Ftype ty = fdpad(sinphi, arho, cosphi, aphi);
		Ftype tz = fdpsb(costh, ar, sinth, ath);

		Ftype axk = load(ax+k);
		Ftype ayk = load(ay+k);
		Ftype azk = load(az+k);
		Ftype pok = load(pot+k);

		axk = fadd(axk, tx);
		ayk = fadd(ayk, ty);
		azk = fadd(azk, tz);
		pok = fadd(pok, potk);

		store(ax+k, axk);
		store(ay+k, ayk);
		store(az+k, azk);
		store(pot+k, pok);
	}
}
