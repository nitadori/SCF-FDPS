// #include <mpi.h>
#include <cmath>
#include "scfdhmpi.h"

void accp_CB(struct Params *param,int ni,int ne,double *mass,
		double *x,double *y,double *z,double *ax,double *ay,
		double *az,double *pot,double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1])
{
	int k,l,m,n;
	double sinth,costh,sinmphi[LLMAX+1],cosmphi[LLMAX+1],phinltil,
	       ar,ath,aphi,r,xi,temp3,temp4,temp5,temp6,clm,dlm,elm,flm,
	       plm[LLMAX+1][LLMAX+1],dplm[LLMAX+1][LLMAX+1],
	       ultrasp[NNMAX+1][LLMAX+1],ultrasp1[NNMAX+1][LLMAX+1],
	       potk,un,unm1,plm1m,plm2m,cosmphi1,sinmphi1,rnorm,rnorm2,
	       rnorm21,ratio,rxy2d,rxy2dinv,r2,r2a2;

	int nmax = param->nmax;
	int lmax = param->lmax;
	int lmin = param->lmin;
	int lskip = param->lskip;
	double G = param->G;
	double abasis = param->abasis;
	double abasis2 = param->abasis2;

	for(k=ni;k<ne;k++){
		rxy2d=x[k]*x[k]+y[k]*y[k];
		rxy2dinv=1.0/sqrt(rxy2d);
		r2=rxy2d+z[k]*z[k];
		r=sqrt(r2);
		costh=z[k]/r;
		xi=(r2-abasis2)/(r2+abasis2);

		cosmphi1=x[k]*rxy2dinv;
		sinmphi1=y[k]*rxy2dinv;
		cosmphi[0]=1.0;
		sinmphi[0]=0.0;
		cosmphi[1]=cosmphi1;
		sinmphi[1]=sinmphi1;

		for(m=2;m<=lmax;m++){
			cosmphi[m]=cosmphi1*cosmphi[m-1]-sinmphi1*sinmphi[m-1];
			sinmphi[m]=cosmphi1*sinmphi[m-1]+sinmphi1*cosmphi[m-1];
		}

		potk=0.0;
		ar=0.0;
		ath=0.0;
		aphi=0.0;

		for(l=0;l<=lmax;l++){
			ultrasp[0][l]=1.0;
			ultrasp[1][l]=twoalpha[l]*xi;
			ultrasp1[0][l]=0.0;
			ultrasp1[1][l]=1.0;
			
			un=ultrasp[1][l];
			unm1=1.0;

			for(n=1;n<=nmax-1;n++){
				ultrasp[n+1][l]=(c1[n][l]*xi*un-c2[n][l]*unm1)
					*c3[n];
				unm1=un;
				un=ultrasp[n+1][l];
				ultrasp1[n+1][l]=((twoalpha[l]+n)*unm1-(n+1)
						*xi*un)/(twoalpha[l]
						*(1.0-xi*xi));
			}
		}

		sinth=sqrt(1.0-costh*costh);
		ratio=1.0;

		for(m=0;m<=lmax;m++){
			plm[m][m]=1.0;
			if(m > 0){
				ratio *= -sinth;
				plm[m][m]=dblfact[m]*ratio;
			}
			plm1m=plm[m][m];
			plm2m=0.0;

			for(l=m+1;l<=lmax;l++){
				plm[l][m]=(costh*(2*l-1)*plm1m-lplusm[l][m]
						*plm2m)*lmm[l][m];
				plm2m=plm1m;
				plm1m=plm[l][m];
			}
		}

		dplm[0][0]=0.0;

		for(l=1;l<=lmax;l++)
			for(m=0;m<=l;m++)
				if(l == m)
					dplm[l][m]=l*costh*plm[l][m]
						/(costh*costh-1.0);
				else
					dplm[l][m]=(l*costh*plm[l][m]-(l+m)
						*plm[l-1][m])/(costh*costh-1.0);
		r2a2=abasis2+r2;
		if(lmin == 0)
			phinltil=abasis/sqrt(r2a2);
		else
			phinltil=abasis2*r/(r2a2*sqrt(r2a2));


		ratio=abasis*r/r2a2;
		if(lskip == 2) ratio *= ratio;

		for(l=lmin;l<=lmax;l += lskip){
			temp3=0.0;
			temp4=0.0;
			temp5=0.0;
			temp6=0.0;

			for(m=0;m<=l;m++){
				clm=0.0;
				dlm=0.0;
				elm=0.0;
				flm=0.0;

				for(n=0;n<=nmax;n++){
					clm += ultrasp[n][l]*cossum[n][l][m];
					dlm += ultrasp[n][l]*sinsum[n][l][m];
					elm += ultrasp1[n][l]*cossum[n][l][m];
					flm += ultrasp1[n][l]*sinsum[n][l][m];
				}

				temp3 += plm[l][m]*(clm*cosmphi[m]
						+dlm*sinmphi[m]);
				temp4 -= plm[l][m]*(elm*cosmphi[m]
						+flm*sinmphi[m]);
				temp5 -= dplm[l][m]*(clm*cosmphi[m]
						+dlm*sinmphi[m]);
				temp6 -= m*plm[l][m]*(dlm*cosmphi[m]
						-clm*sinmphi[m]);
			}

			rnorm=r/abasis;
			rnorm2=rnorm*rnorm;
			rnorm21=rnorm2+1.0;
			potk += temp3*phinltil;
			ar += phinltil*(-temp3*(l/rnorm-(2*l+1)*rnorm/rnorm21)
				+temp4*4.0*rnorm*twoalpha[l]/(rnorm21*rnorm21));
			ath += temp5*phinltil;
			aphi += temp6*phinltil;
			phinltil *= ratio;
		}

		ar /= abasis2;
		ath *= -sinth/(abasis*r);
		aphi /= (abasis*r*sinth);
		ax[k] += G*(sinth*cosmphi1*ar+costh*cosmphi1*ath
				-sinmphi1*aphi);
		ay[k] += G*(sinth*sinmphi1*ar+costh*sinmphi1*ath
				+cosmphi1*aphi);
		az[k] += G*(costh*ar-sinth*ath);
		pot[k] += G*potk/abasis;
	}
}

void accp_LH(struct Params *param,int ni,int ne,double *mass,
		double *x,double *y,double *z,double *ax,double *ay,
		double *az,double *pot,double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1])
{
	int k,l,m,n;
	double sinth,costh,sinmphi[LLMAX+1],cosmphi[LLMAX+1],phinltil,
	       ar,ath,aphi,r,xi,temp3,temp4,temp5,temp6,clm,dlm,elm,flm,
	       plm[LLMAX+1][LLMAX+1],dplm[LLMAX+1][LLMAX+1],
	       ultrasp[NNMAX+1][LLMAX+1],ultrasp1[NNMAX+1][LLMAX+1],
	       potk,un,unm1,plm1m,plm2m,cosmphi1,sinmphi1,rnorm,rnorm1,
	       ratio,rxy2d,rxy2dinv,ra;

	int nmax = param->nmax;
	int lmax = param->lmax;
	int lmin = param->lmin;
	int lskip = param->lskip;
	double G = param->G;
	double abasis = param->abasis;
	double abasis2 = param->abasis2;

	for(k=ni;k<ne;k++){
		rxy2d=x[k]*x[k]+y[k]*y[k];
		rxy2dinv=1.0/sqrt(rxy2d);
		r=sqrt(rxy2d+z[k]*z[k]);
		costh=z[k]/r;
		sinth=sqrt(rxy2d)/r;
		xi=(r-abasis)/(r+abasis);

		cosmphi1=x[k]*rxy2dinv;
		sinmphi1=y[k]*rxy2dinv;
		cosmphi[0]=1.0;
		sinmphi[0]=0.0;
		cosmphi[1]=cosmphi1;
		sinmphi[1]=sinmphi1;

		for(m=2;m<=lmax;m++){
			cosmphi[m]=cosmphi1*cosmphi[m-1]-sinmphi1*sinmphi[m-1];
			sinmphi[m]=cosmphi1*sinmphi[m-1]+sinmphi1*cosmphi[m-1];
		}

		potk=0.0;
		ar=0.0;
		ath=0.0;
		aphi=0.0;

		for(l=0;l<=lmax;l++){
			ultrasp[0][l]=1.0;
			ultrasp[1][l]=twoalpha[l]*xi;
			ultrasp1[0][l]=0.0;
			ultrasp1[1][l]=1.0;
			
			un=ultrasp[1][l];
			unm1=1.0;

			for(n=1;n<=nmax-1;n++){
				ultrasp[n+1][l]=(c1[n][l]*xi*un-c2[n][l]*unm1)
						*c3[n];
				unm1=un;
				un=ultrasp[n+1][l];
				ultrasp1[n+1][l]=((twoalpha[l]+n)*unm1
						-(n+1.0)*xi*un)
						/(twoalpha[l]*(1.0-xi*xi));
			}
		}

		ratio=1.0;

		for(m=0;m<=lmax;m++){
			plm[m][m]=1.0;
			if(m > 0){
				ratio *= -sinth;
				plm[m][m]=dblfact[m]*ratio;
			}
			plm1m=plm[m][m];
			plm2m=0.0;

			for(l=m+1;l<=lmax;l++){
				plm[l][m]=(costh*(2*l-1)*plm1m
					-lplusm[l][m]*plm2m)*lmm[l][m];
				plm2m=plm1m;
				plm1m=plm[l][m];
			}
		}

		dplm[0][0]=0.0;

		for(l=1;l<=lmax;l++)
			for(m=0;m<=l;m++)
				if(l == m)
					dplm[l][m]=l*costh*plm[l][m]
						/(costh*costh-1.0);
				else
					dplm[l][m]=(l*costh*plm[l][m]
						-(l+m)*plm[l-1][m])
						/(costh*costh-1.0);
		ra=abasis+r;
		if(lmin == 0)
			phinltil=abasis/ra;
		else
			phinltil=abasis2*r/(ra*ra*ra);


		ratio=abasis*r/(ra*ra);
		if(lskip == 2) ratio *= ratio;

		for(l=lmin;l<=lmax;l += lskip){
			temp3=0.0;
			temp4=0.0;
			temp5=0.0;
			temp6=0.0;

			for(m=0;m<=l;m++){
				clm=0.0;
				dlm=0.0;
				elm=0.0;
				flm=0.0;

				for(n=0;n<=nmax;n++){
					clm += ultrasp[n][l]*cossum[n][l][m];
					dlm += ultrasp[n][l]*sinsum[n][l][m];
					elm += ultrasp1[n][l]*cossum[n][l][m];
					flm += ultrasp1[n][l]*sinsum[n][l][m];
				}

				temp3 += plm[l][m]*(clm*cosmphi[m]
					+dlm*sinmphi[m]);
				temp4 -= plm[l][m]*(elm*cosmphi[m]
					+flm*sinmphi[m]);
				temp5 -= dplm[l][m]*(clm*cosmphi[m]
					+dlm*sinmphi[m]);
				temp6 -= m*plm[l][m]*(dlm*cosmphi[m]
					-clm*sinmphi[m]);
			}

			rnorm=r/abasis;
			rnorm1=rnorm+1.0;
			potk += temp3*phinltil;
			ar += phinltil*(-temp3*(l/rnorm-(2*l+1)/rnorm1)
				+temp4*4.0*(2*l+1.5)/(rnorm1*rnorm1));
			ath += temp5*phinltil;
			aphi += temp6*phinltil;
			phinltil *= ratio;
		}

		ar /= abasis;
		ath *= -sinth/r;
		aphi /= (r*sinth);
		ax[k] += G*(sinth*cosmphi1*ar+costh*cosmphi1*ath
				-sinmphi1*aphi);
		ay[k] += G*(sinth*sinmphi1*ar+costh*sinmphi1*ath
				+cosmphi1*aphi);
		az[k] += G*(costh*ar-sinth*ath);
		pot[k] += G*potk;
	}
}

void scfcoef_CB(struct Params *param,int ni,int ne,
		// double sinsum[][LLMAX+1][LLMAX+1],
		// double cossum[][LLMAX+1][LLMAX+1],
		double (* __restrict sinsum)[LLMAX+1][LLMAX+1],
		double (* __restrict cossum)[LLMAX+1][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double c1[][LLMAX+1],double c2[][LLMAX+1],
		double *c3,double anltilde[][LLMAX+1],double *dblfact,
		double *twoalpha,double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1])
{
	int k,l,m,n,icount;
	double rxy2d,rxy2dinv,r2,r,costh,sinth,xi,un,unm1,ratio,a2r2,
		plm1m,plm2m,temp3,temp4,temp5,ttemp5;
	double cosmphi[LLMAX+1], sinmphi[LLMAX+1], plm[LLMAX+1][LLMAX+1];
	double ultrasp[NNMAX+1][LLMAX+1],ultraspt[NNMAX+1][LLMAX+1];

	int nmax = param->nmax;
	int lmax = param->lmax;
	int lmin = param->lmin;
	int lskip = param->lskip;
	double abasis = param->abasis;
	double abasis2 = param->abasis2;
	// static double tmpsumsnd[(NNMAX+1)*(LLMAX+1)*(LLMAX+2)];
	// static double tmpsumrcv[(NNMAX+1)*(LLMAX+1)*(LLMAX+2)];
	static double tmpsum[(NNMAX+1)*(LLMAX+1)*(LLMAX+2)];

	for(k=ni;k<ne;k++){
		rxy2d=x[k]*x[k]+y[k]*y[k];
		rxy2dinv=1.0/sqrt(rxy2d);
		r2=rxy2d+z[k]*z[k];
		r=sqrt(r2);
		costh=z[k]/r;
		sinth=sqrt(1.0-costh*costh);
		xi=(r2-abasis2)/(r2+abasis2);
		cosmphi[0]=1.0;
		sinmphi[0]=0.0;
		cosmphi[1]=x[k]*rxy2dinv;
		sinmphi[1]=y[k]*rxy2dinv;

		for(m=2;m<=lmax;m++){
			cosmphi[m]=cosmphi[1]*cosmphi[m-1]
				-sinmphi[1]*sinmphi[m-1];
			sinmphi[m]=cosmphi[1]*sinmphi[m-1]
				+sinmphi[1]*cosmphi[m-1];
		}

		for(l=0;l<=lmax;l++){
			ultrasp[0][l]=1.0;
			ultrasp[1][l]=twoalpha[l]*xi;
			un=ultrasp[1][l];
			unm1=1.0;

			for(n=1;n<=nmax-1;n++){
				ultrasp[n+1][l]=(c1[n][l]*xi*un-c2[n][l]
						*unm1)*c3[n];
				unm1=un;
				un=ultrasp[n+1][l];
			}

			for(n=0;n<=nmax;n++)
				ultraspt[n][l]=ultrasp[n][l]*anltilde[n][l];
		}

		ratio=1.0;
		for(m=0;m<=lmax;m++){
			plm[m][m]=1.0;
			if(m > 0){
				ratio *= -sinth;
				plm[m][m]=dblfact[m]*ratio;
			}
			plm1m=plm[m][m];
			plm2m=0.0;

			for(l=m+1;l<=lmax;l++){
				plm[l][m]=(costh*(2*l-1)*plm1m
						-lplusm[l][m]*plm2m)*lmm[l][m];
				plm2m=plm1m;
				plm1m=plm[l][m];
			}
		}

		a2r2 = abasis2+r2;
		if(lmin == 0)
			temp5=mass[k]*abasis/sqrt(a2r2);
		else
			temp5=mass[k]*abasis2*r/((a2r2)*sqrt(a2r2));

		ratio=abasis*r/a2r2;

		if(lskip == 2) ratio *= ratio;

		for(l=lmin;l<=lmax;l += lskip){
			for(m=0;m<=l;m++){
				ttemp5=temp5*plm[l][m]*coeflm[l][m];
				temp3=ttemp5*sinmphi[m];
				temp4=ttemp5*cosmphi[m];
				for(n=0;n<=nmax;n++){
					sinsum[n][l][m] += temp3*ultraspt[n][l];
					cossum[n][l][m] += temp4*ultraspt[n][l];
				}
			}
			temp5 *= ratio;
		}
	}
	icount=0;
	for(l=lmin;l<=lmax;l += lskip)
		for(m=0;m<=l;m++)
			for(n=0;n<=nmax;n++){
				icount += 2;	
				tmpsum[icount-2]=sinsum[n][l][m];
				tmpsum[icount-1]=cossum[n][l][m];
			}
	// MPI_Allreduce(tmpsumsnd,tmpsumrcv,icount,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
	scfmpi_sum_double(tmpsum, icount);

	icount=0;
	for(l=lmin;l<=lmax;l += lskip)
		for(m=0;m<=l;m++)
			for(n=0;n<=nmax;n++){
				icount += 2;	
				sinsum[n][l][m]=tmpsum[icount-2];
				cossum[n][l][m]=tmpsum[icount-1];
			}
}

void scfcoef_LH(struct Params *param,int ni,int ne,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double c1[][LLMAX+1],double c2[][LLMAX+1],
		double *c3,double anltilde[][LLMAX+1],double *dblfact,
		double *twoalpha,double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1])
{
	int k,l,m,n,icount;
	double rxy2d,rxy2dinv,r,costh,sinth,xi,cosmphi[LLMAX+1],
		sinmphi[LLMAX+1],ultrasp[NNMAX+1][LLMAX+1],
		ultraspt[NNMAX+1][LLMAX+1],un,unm1,ratio,ar,
		plm[LLMAX+1][LLMAX+1],plm1m,plm2m,temp3,temp4,temp5,ttemp5;

	int nmax = param->nmax;
	int lmax = param->lmax;
	int lmin = param->lmin;
	int lskip = param->lskip;
	double abasis = param->abasis;
	// static double tmpsumsnd[(NNMAX+2)*(LLMAX+1)*(LLMAX+2)];
	// static double tmpsumrcv[(NNMAX+1)*(LLMAX+1)*(LLMAX+2)];
	static double tmpsum[(NNMAX+1)*(LLMAX+1)*(LLMAX+2)];

	for(k=ni;k<ne;k++){
		rxy2d=x[k]*x[k]+y[k]*y[k];
		rxy2dinv=1.0/sqrt(rxy2d);
		r=sqrt(rxy2d+z[k]*z[k]);
		costh=z[k]/r;
		sinth=sqrt(1.0-costh*costh);
		xi=(r-abasis)/(r+abasis);
		cosmphi[0]=1.0;
		sinmphi[0]=0.0;
		cosmphi[1]=x[k]*rxy2dinv;
		sinmphi[1]=y[k]*rxy2dinv;

		for(m=2;m<=lmax;m++){
			cosmphi[m]=cosmphi[1]*cosmphi[m-1]
				-sinmphi[1]*sinmphi[m-1];
			sinmphi[m]=cosmphi[1]*sinmphi[m-1]
				+sinmphi[1]*cosmphi[m-1];
		}

		for(l=0;l<=lmax;l++){
			ultrasp[0][l]=1.0;
			ultrasp[1][l]=twoalpha[l]*xi;
			un=ultrasp[1][l];
			unm1=1.0;

			for(n=1;n<=nmax-1;n++){
				ultrasp[n+1][l]=(c1[n][l]*xi*un
						-c2[n][l]*unm1)*c3[n];
				unm1=un;
				un=ultrasp[n+1][l];
			}

			for(n=0;n<=nmax;n++)
				ultraspt[n][l]=ultrasp[n][l]*anltilde[n][l];
		}

		ratio=1.0;
		for(m=0;m<=lmax;m++){
			plm[m][m]=1.0;
			if(m > 0){
				ratio *= -sinth;
				plm[m][m]=dblfact[m]*ratio;
			}
			plm1m=plm[m][m];
			plm2m=0.0;

			for(l=m+1;l<=lmax;l++){
				plm[l][m]=(costh*(2*l-1)*plm1m-lplusm[l][m]
						*plm2m)*lmm[l][m];
				plm2m=plm1m;
				plm1m=plm[l][m];
			}
		}

		ar = abasis+r;
		if(lmin == 0)
			temp5=mass[k]/ar;
		else
			temp5=mass[k]*abasis*r/(ar*ar*ar);

		ratio=abasis*r/(ar*ar);

		if(lskip == 2) ratio *= ratio;

		for(l=lmin;l<=lmax;l += lskip){
			for(m=0;m<=l;m++){
				ttemp5=temp5*plm[l][m]*coeflm[l][m];
				temp3=ttemp5*sinmphi[m];
				temp4=ttemp5*cosmphi[m];
				for(n=0;n<=nmax;n++){
					sinsum[n][l][m] += temp3*ultraspt[n][l];
					cossum[n][l][m] += temp4*ultraspt[n][l];
				}
			}
			temp5 *= ratio;
		}

	}
	icount=0;
	for(l=lmin;l<=lmax;l += lskip)
		for(m=0;m<=l;m++)
			for(n=0;n<=nmax;n++){
				icount += +2;	
				tmpsum[icount-2]=sinsum[n][l][m];
				tmpsum[icount-1]=cossum[n][l][m];
			}
	// MPI_Allreduce(tmpsumsnd,tmpsumrcv,icount,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
	scfmpi_sum_double(tmpsum, icount);

	icount=0;
	for(l=lmin;l<=lmax;l += lskip)
		for(m=0;m<=l;m++)
			for(n=0;n<=nmax;n++){
				icount += 2;	
				sinsum[n][l][m]=tmpsum[icount-2];
				cossum[n][l][m]=tmpsum[icount-1];
			}
}

