#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <omp.h>
//omp_get_wtime();

#define NMAX 60000
#define NNMAX 64
#define LLMAX 32

struct Params{
        enum {TEXT_LEN = 8} len;
        int incont;
        int initsnap;
        char basisset[TEXT_LEN];
        int nmax;
        int lmax;
        int lmin;
        int lskip;
        int nsteps;
        int noutbod;
        int noutlog;
        double dtime;
        double G;
        double abasis;
        double abasis2;
        int outpcoef;
        int zeroodd;
        int zeroeven;
        int fixacc;
        int ubodlog;
    
        char softening_type[TEXT_LEN];
        double theta;
        int n_leaf_limit;
        int n_group_limit;
    
};

void accp_CB(struct Params *param,int nbodies,double *mass,
		double *x,double *y,double *z,double *ax,double *ay,
		double *az,double *pot,double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1]);

void accp_LH(struct Params *param,int nbodies,double *mass,
		double *x,double *y,double *z,double *ax,double *ay,
		double *az,double *pot,double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1]);

void accpot(struct Params *param,int nbodies,
		double *tnow,double *mass,double *x,double *y,double *z,
		double *ax,double *ay,double *az,double *pot,
		double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],
		double lplusm[][LLMAX+1],double lmm[][LLMAX+1],
		double coeflm[][LLMAX+1],int me);

void corracc(int nbodies,double *mass,double *ax,double *ay,double *az);

void corrvel(double rcsign,int nbodies,double *tvel,double *tnow,double dtime,
		double *vx,double *vy,double *vz,double *ax,double *ay,
		double *az);

void endrun(double totime,int me);

double factrl(int n);

double gammln(double xx);

void getcoef_CB(int nmax,int lmax,double pi,double *dblfact,
	double anltilde[][LLMAX+1],double coeflm[][LLMAX+1],double *twoalpha,
	double lmm[][LLMAX+1],double lplusm[][LLMAX+1],double c1[][LLMAX+1],
	double c2[][LLMAX+1],double *c3);

void getcoef_LH(int nmax,int lmax,double pi,double *dblfact,
	double anltilde[][LLMAX+1],double coeflm[][LLMAX+1],double *twoalpha,
	double lmm[][LLMAX+1],double lplusm[][LLMAX+1],double c1[][LLMAX+1],
	double c2[][LLMAX+1],double *c3);

void inbods(int incont,int *nbodies,double *tnow,double *mass,
		double *x,double *y,double *z,double *vx,double *vy,
		double *vz,double *etotal,int me,int n_pes);

void initpars(double *tpos,double *tvel,double *tnow);

void initsys(struct Params *param,int nbodies,
		double *tnow,double *tpos,double *tvel,double *dblfact,
		double *twoalpha,double c1[][LLMAX+1],
		double c2[][LLMAX+1],double *c3,double anltilde[][LLMAX+1],
		double lplusm[][LLMAX+1],double lmm[][LLMAX+1],
		double coeflm[][LLMAX+1],
		double *mass,double *x,double *y, double *z,
		double *vx,double *vy,double *vz,
		double *ax,double *ay,double *az,double *pot,
		double *etotal,double cputime,int me,int n_pes);

void initvel(int nbodies,double *vx,double *vy,double *vz,double *ax,
		double *ay,double *az,double dtime,double *tnow,double *tvel);

void inparams(struct Params *param,int me);

void mpiget(double *rbuf,double *sbuf,int n,int idest,int isrc,int me);

void writecoef(double *tnow,int nmax,int lmax,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],int me);

void outbods(int nbodies,int initsnap,double *tnow,double *mass,double *x,
		double *y,double *z,double *vx,double *vy,double *vz,
		double *pot,double *etotal,int ubodlog,int me,int n_pes);

void outlog(struct Params *param,int nbodies,double *tnow,
		double *mass,double *x,double *y,double *z,double *vx,
		double *vy,double *vz,double *ax,double *ay,double *az,
		double *pot,double *etotal,double cputime0,int me,int n_pes);

void outstate(struct Params *param,int kstep,int nbodies,double *tvel,
		double *tnow,double *mass,double *x,double *y,double *z,
		double *vx,double *vy,double *vz,double *ax,double *ay,
		double *az,double *pot,double *etotal,double cputime,
		int me,int n_pes);

void outterm(char *message,int n,int me);

void scfcoef_CB(struct Params *param,int nbodies,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double c1[][LLMAX+1],double c2[][LLMAX+1],
		double *c3,double anltilde[][LLMAX+1],double *dblfact,
		double *twoalpha,double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1]);

void scfcoef_LH(struct Params *param,int nbodies,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double c1[][LLMAX+1],double c2[][LLMAX+1],
		double *c3,double anltilde[][LLMAX+1],double *dblfact,
		double *twoalpha,double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1]);

void share(double *x,int n);

void shareint(int *x,int n);

void startout(int me);

void steppos(int nbodies,double *x,double *y,double *z,double *vx,
		double *vy,double *vz,double dtime,double *tnow,
		double *tpos);

void stepsys(struct Params *param,int kstep,int nbodies,
		double *tnow,double *tpos,double *tvel,
		double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double *vx,double *vy,double *vz,
		double *ax,double *ay,double *az,double *pot,
		double *etotal,double cputime,int me,int n_pes);

void stepvel(int nbodies,double *vx,double *vy,double *vz,double *ax,
		double *ay,double *az,double dtime,double *tnow,
		double *tvel);

int main(int argc,char **argv)
{
	struct Params param;
	double t0mpi,t1mpi,t2mpi,tmpi;
	int n,nbodies;
	int me,n_pes;
	double tnow,tpos,tvel,etotal;
	double mass[NMAX],x[NMAX],y[NMAX],z[NMAX],vx[NMAX],
		vy[NMAX],vz[NMAX],ax[NMAX],ay[NMAX],az[NMAX],
		pot[NMAX];
	double dblfact[LLMAX+1],twoalpha[LLMAX+1],c1[LLMAX+1][LLMAX+1],
		c2[LLMAX+1][LLMAX+1],c3[NNMAX+1],anltilde[NNMAX+1][LLMAX+1],
		lplusm[LLMAX+1][LLMAX+1],lmm[LLMAX+1][LLMAX+1],
		coeflm[LLMAX+1][LLMAX+1];

	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD,&n_pes);

	MPI_Comm_rank(MPI_COMM_WORLD,&me);

	t0mpi = MPI_Wtime();

	startout(me);

	inparams(&param,me);

	int incont = param.incont;
	int initsnap = param.initsnap;
	char *basisset = param.basisset;
	int nmax = param.nmax;
	int lmax = param.lmax;
	int lmin = param.lmin;
	int lskip = param.lskip;
	int nsteps = param.nsteps;
	int noutbod = param.noutbod;
	int noutlog = param.noutlog;
	double dtime = param.dtime;
	double G = param.G;
	double abasis = param.abasis;
	int outpcoef = param.outpcoef;
	int zeroodd = param.zeroodd;
	int zeroeven = param.zeroeven;
	int fixacc = param.fixacc;
	int ubodlog = param.ubodlog;

	(void)incont;
	(void)initsnap;
	(void)basisset;
	(void)nmax;
	(void)lmax;
	(void)lmin;
	(void)lskip;
	(void)nsteps;
	(void)noutbod;
	(void)noutlog;
	(void)dtime;
	(void)G;
	(void)abasis;
	(void)outpcoef;
	(void)zeroodd;
	(void)zeroeven;
	(void)fixacc;
	(void)ubodlog;

	inbods(incont,&nbodies,&tnow,mass,x,y,z,vx,vy,vz,&etotal,me,n_pes);

	initsys(&param,nbodies,&tnow,&tpos,&tvel,dblfact,twoalpha,
		c1,c2,c3,anltilde,lplusm,lmm,coeflm,
		mass,x,y,z,vx,vy,vz,ax,ay,az,pot,&etotal,t0mpi,me,n_pes);

	t1mpi = MPI_Wtime();
	tmpi = t1mpi - t0mpi;

	for(n=1;n<=nsteps;n++)
		stepsys(&param,n,nbodies,&tnow,&tpos,&tvel,dblfact,
			twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm,
			mass,x,y,z,vx,vy,vz,ax,ay,az,pot,&etotal,
			t0mpi,me,n_pes);

	t2mpi = MPI_Wtime();
	tmpi = t2mpi-t1mpi;

	endrun(tmpi,me);

	MPI_Finalize();

	return 0;
}

double factrl(int n)
{
	int ntop,j;
	double a[34],arggam,ftrl;

	ntop=0;
	a[1]=1.0;

	if(n < 0){
		printf("negative factorial\n");
		exit(1);
	}else if(n <= ntop)
		ftrl=a[n+1];
	else if(n <= 32){
		for(j=ntop+1;j<=n;j++)
			a[j+1]=j*a[j];
		ntop=n;
		ftrl=a[n+1];
	}else{
		arggam=n+1.0;
		ftrl=exp(gammln(arggam));
	}

	return ftrl;
}	

void accp_CB(struct Params *param,int nbodies,double *mass,
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

	for(k=0;k<nbodies;k++){
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

			for(n=1;n<nmax;n++){
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

void accp_LH(struct Params *param,int nbodies,double *mass,
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

	for(k=0;k<nbodies;k++){
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

			for(n=1;n<nmax;n++){
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

void accpot(struct Params *param,int nbodies,
		double *tnow,double *mass,double *x,double *y,double *z,
		double *ax,double *ay,double *az,double *pot,
		double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],
		double lplusm[][LLMAX+1],double lmm[][LLMAX+1],
		double coeflm[][LLMAX+1],int me)
{
	int k,l,m,n;

	char *basisset = param->basisset;
	int nmax = param->nmax;
	int lmax = param->lmax;
	int outpcoef = param->outpcoef;
	double sinsum[NNMAX+1][LLMAX+1][LLMAX+1];
	double cossum[NNMAX+1][LLMAX+1][LLMAX+1];

	for(l=0;l<=lmax;l++)
		for(m=0;m<=l;m++)
			for(n=0;n<=nmax;n++){
				sinsum[n][l][m]=0.0;
				cossum[n][l][m]=0.0;
			}

	if(strcmp(basisset,"CB") == 0)
		scfcoef_CB(param,nbodies,sinsum,cossum,mass,x,y,z,c1,c2,c3,
				anltilde,dblfact,twoalpha,lplusm,lmm,coeflm);
	else if(strcmp(basisset,"LH")==0)
		scfcoef_LH(param,nbodies,sinsum,cossum,mass,x,y,z,c1,c2,c3,
				anltilde,dblfact,twoalpha,lplusm,lmm,coeflm);

	if(outpcoef == 1)
		writecoef(tnow,nmax,lmax,sinsum,cossum,me);
	for(k=0;k<=nbodies;k++){
		ax[k]=0.0;
		ay[k]=0.0;
		az[k]=0.0;
		pot[k]=0.0;
	}

	if(strcmp(basisset,"CB") == 0)
		accp_CB(param,nbodies,mass,x,y,z,ax,ay,az,pot,dblfact,
			twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm,
			sinsum,cossum);
	else if(strcmp(basisset,"LH") == 0)
		accp_LH(param,nbodies,mass,x,y,z,ax,ay,az,pot,dblfact,
			twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm,
			sinsum,cossum);
}

void corracc(int nbodies,double *mass,double *ax,double *ay,double *az)
{
	int k;
	double mtot=0.0,axcm=0.0,aycm=0.0,azcm=0.0;
	double tmpsumsnd[4],tmpsumrcv[4];

	for(k=0;k<nbodies;k++){
		mtot += mass[k];
		axcm += mass[k]*ax[k];
		aycm += mass[k]*ay[k];
		azcm += mass[k]*az[k];
	}

	tmpsumsnd[0]=mtot;
	tmpsumsnd[1]=axcm;
	tmpsumsnd[2]=aycm;
	tmpsumsnd[3]=azcm;

	MPI_Allreduce(tmpsumsnd,tmpsumrcv,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	mtot = tmpsumrcv[0];
	axcm = tmpsumrcv[1];
	aycm = tmpsumrcv[2];
	azcm = tmpsumrcv[3];

	axcm /= mtot;
	aycm /= mtot;
	azcm /= mtot;

	for(k=0;k<=nbodies;k++){
		ax[k] -= axcm;
		ay[k] -= aycm;
		az[k] -= azcm;
	}
}

void corrvel(double rcsign,int nbodies,double *tvel,double *tnow,double dtime,
		double *vx,double *vy,double *vz,double *ax,double *ay,
		double *az)
{
	int k;

	for(k=0;k<nbodies;k++){
		vx[k] += rcsign*ax[k]*0.5*dtime;
		vy[k] += rcsign*ay[k]*0.5*dtime;
		vz[k] += rcsign*az[k]*0.5*dtime;
	}

	*tvel += rcsign*0.5*dtime;
	*tnow = *tvel;
}

void endrun(double totime,int me)
{
	FILE *fp;

	if(me == 0){
		fp = fopen("SCFLOG","a");
		fprintf(fp,"     Total cpu time used (secs) = %15.7E\n",totime);

		fclose(fp);
	}
}

double gammln(double xx)
{
	int j;
	double cof[7],stp,half,one,fpf,x,tmp,ser,gam;

	cof[1]=76.18009173;
	cof[2]=-86.50532033;
	cof[3]=24.01409822;
	cof[4]=-1.231739516;
	cof[5]=0.120858003e-2;
	cof[6]=-0.536382e-5,
	stp=2.50662827465;
	half=0.5;
	one=1.0;
	fpf=5.5;

	x=xx-one;
	tmp=x+fpf;
	tmp=(x+half)*log(tmp)-tmp;
	ser=one;

	for(j=1;j<=6;j++){
		x=x+one;
		ser=ser+cof[j]/x;
	}

	gam=tmp+log(stp*ser);

	return gam;
}

void getcoef_CB(int nmax,int lmax,double pi,double *dblfact,
	double anltilde[][LLMAX+1],double coeflm[][LLMAX+1],double *twoalpha,
	double lmm[][LLMAX+1],double lplusm[][LLMAX+1],double c1[][LLMAX+1],
	double c2[][LLMAX+1],double *c3)
{
        int l,m,n;
        double deltam0,knl,fctrll,fctrln;

        dblfact[1]=1.0;

        if(lmax >= 2)
                for(l=2;l<=lmax;l++)
                        dblfact[l]=dblfact[l-1]*(2*l-1);

        for(n=0;n<=nmax;n++){
                fctrln=factrl(n);
                for(l=0;l<=lmax;l++){
                        knl=4*n*(n+2*l+2)+(2*l+1)*(2*l+3);
                        fctrll=factrl(l);
                        anltilde[n][l]=-pow(2,4*l+4)*((n+l+1)/(pi*knl))
                                *(fctrln*fctrll/factrl(n+2*l+1))*fctrll;
                }
        }

        for(l=0;l<=lmax;l++){
                twoalpha[l]=2*(l+1);
                for(m=0;m<=l;m++){
                        deltam0=2.0;
                        if(m == 0) deltam0=1.0;
                        coeflm[l][m]=(2*l+1)*deltam0*factrl(l-m)/factrl(l+m);
                        lmm[l][m]=l-m;
                        if(l != m) lmm[l][m]=1.0/lmm[l][m];
                        lplusm[l][m]=l+m-1.0;
                }
        }

        for(n=1;n<=nmax;n++){
                c3[n]=1.0/(n+1.0);
                for(l=0;l<=lmax;l++){
                        c1[n][l]=2*n+twoalpha[l];
                        c2[n][l]=n-1+twoalpha[l];
                }
        }
}

void getcoef_LH(int nmax,int lmax,double pi,double *dblfact,
	double anltilde[][LLMAX+1],double coeflm[][LLMAX+1],double *twoalpha,
	double lmm[][LLMAX+1],double lplusm[][LLMAX+1],double c1[][LLMAX+1],
	double c2[][LLMAX+1],double *c3)
{
        int l,m,n;
        double deltam0,knl,arggam,fctrln,fctrl2l;

        dblfact[1]=1.0;

        if(lmax >= 2)
                for(l=2;l<=lmax;l++)
                        dblfact[l]=dblfact[l-1]*(2*l-1);

        for(n=0;n<=nmax;n++){
                fctrln=factrl(n);
                for(l=0;l<=lmax;l++){
                        knl=0.5*n*(n+4*l+3)+(l+1)*(2*l+1);
			arggam=2*l+1.5;
                        fctrl2l=exp(gammln(arggam));
                        anltilde[n][l]=-pow(2,8*l+4)*((n+2*l+1.5)/(pi*knl))
                                *(fctrln*fctrl2l/factrl(n+4*l+2))*fctrl2l;
                }
        }

        for(l=0;l<=lmax;l++){
                twoalpha[l]=4*l+3;
                for(m=0;m<=l;m++){
                        deltam0=2.0;
                        if(m == 0) deltam0=1.0;
                        coeflm[l][m]=(2*l+1)*deltam0*factrl(l-m)
				/factrl(l+m);
                        lmm[l][m]=l-m;
                        if(l != m) lmm[l][m]=1.0/lmm[l][m];
                        lplusm[l][m]=l+m-1.0;
                }
        }

        for(n=1;n<=nmax;n++){
                c3[n]=1.0/(n+1);
                for(l=0;l<=lmax;l++){
                        c1[n][l]=2*n+twoalpha[l];
                        c2[n][l]=n-1+twoalpha[l];
                }
        }
}

void initpars(double *tpos,double *tvel,double *tnow)
{
	*tpos=*tnow;
	*tvel=*tnow;
}

void initsys(struct Params *param,int nbodies,
		double *tnow,double *tpos,double *tvel,double *dblfact,
		double *twoalpha,double c1[][LLMAX+1],
		double c2[][LLMAX+1],double *c3,double anltilde[][LLMAX+1],
		double lplusm[][LLMAX+1],double lmm[][LLMAX+1],
		double coeflm[][LLMAX+1],
		double *mass,double *x,double *y, double *z,
		double *vx,double *vy,double *vz,
		double *ax,double *ay,double *az,double *pot,
		double *etotal,double cputime,int me,int n_pes)
{
	int kstep=0;
	double pi;

	int lmax = param->lmax;
	int nmax = param->nmax;
	char *basisset = param->basisset;
	double dtime = param->dtime;
	int fixacc = param->fixacc;

	pi=4.0*atan(1.0);
	initpars(tpos,tvel,tnow);

	if(strcmp(basisset,"CB") == 0)
		getcoef_CB(nmax,lmax,pi,dblfact,anltilde,coeflm,twoalpha,
			lmm,lplusm,c1,c2,c3);
	else if(strcmp(basisset,"LH") ==0)
		getcoef_LH(nmax,lmax,pi,dblfact,anltilde,coeflm,twoalpha,
			lmm,lplusm,c1,c2,c3);

	accpot(param,nbodies,tnow,mass,x,y,z,ax,ay,az,pot,dblfact,twoalpha,
		c1,c2,c3,anltilde,lplusm,lmm,coeflm,me);

	if(fixacc == 1) corracc(nbodies,mass,ax,ay,az);

	outstate(param,kstep,nbodies,tvel,tnow,mass,x,y,z,vx,vy,vz,
			ax,ay,az,pot,etotal,cputime,me,n_pes);

	initvel(nbodies,vx,vy,vz,ax,ay,az,dtime,tnow,tvel);
}

void initvel(int nbodies,double *vx,double *vy,double *vz,double *ax,
		double *ay,double *az,double dtime,double *tnow,double *tvel)
{
	int k;

	for(k=0;k<nbodies;k++){
		vx[k] += 0.5*dtime*ax[k];
		vy[k] += 0.5*dtime*ay[k];
		vz[k] += 0.5*dtime*az[k];
	}

	*tvel += 0.5*dtime;
	*tnow = *tvel;
}

void inbods(int incont,int *nbodies,double *tnow,double *mass,
		double *x,double *y,double *z,double *vx,double *vy,
		double *vz,double *etotal,int me,int n_pes)
{
        int i,k,nb;
        double tt,etot,mm,xx,yy,zz,vxx,vyy,vzz;
	double temp01[NMAX],temp02[NMAX],temp03[NMAX],temp04[NMAX],
		temp05[NMAX],temp06[NMAX],temp07[NMAX];
        FILE *fp;

	if(me == 0){
		if((fp=fopen("SCFBI","r"))==NULL){
			printf("SCFBI not found.\n");
			exit(1);
		}

		if(incont == 0){
			fscanf(fp,"%d %lf",&nb,&tt);
			*nbodies=nb;
			*tnow=tt;
		}else{
			fscanf(fp,"%d %lf %lf",&nb,&tt,&etot);
			*nbodies=nb;
			*tnow=tt;
			*etotal=etot;
		}
	}
	shareint(nbodies,1);
	share(tnow,1);
	*nbodies=*nbodies/n_pes;
	if(me == 0){
        	for(k=0;k<*nbodies;k++){
                	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",
                        	&mm,&xx,&yy,&zz,&vxx,&vyy,&vzz);
                	mass[k]=mm;
                	x[k]=xx;
                	y[k]=yy;
                	z[k]=zz;
                	vx[k]=vxx;
                	vy[k]=vyy;
                	vz[k]=vzz;
		}
        }

	for(i=1;i<n_pes;i++){
		if(me == 0){
			for(k=0;k<*nbodies;k++){
				fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",
					&mm,&xx,&yy,&zz,&vxx,&vyy,&vzz);
				temp01[k]=mm;
				temp02[k]=xx;
				temp03[k]=yy;
				temp04[k]=zz;
				temp05[k]=vxx;
				temp06[k]=vyy;
				temp07[k]=vzz;
			}
		}
		mpiget(mass,temp01,*nbodies,i,0,me);
		mpiget(x   ,temp02,*nbodies,i,0,me);
		mpiget(y   ,temp03,*nbodies,i,0,me);
		mpiget(z   ,temp04,*nbodies,i,0,me);
		mpiget(vx  ,temp05,*nbodies,i,0,me);
		mpiget(vy  ,temp06,*nbodies,i,0,me);
		mpiget(vz  ,temp07,*nbodies,i,0,me);
	}

	if(me == 0)
		fclose(fp);
}

void inparams(struct Params *param,int me)
{
//************************* Input Parameters *****************************
//
//	incont    : option to choose new simulation (0) or continued
//	           simulation (1).
//	initsnap  : the starting number of snap files.
//	basisset  : option to choose Clutton-Brock's basis set (CB) or
//	           Hernquist-Ostriker's basis set (LH).
//	nmax      : number of radial eigenfunctions.
//	lmax      : nubmer of angular eigenfunctions.
//	nsteps	  : number of timesteps.
//	noutbod	  : output system state once every nsteps/noutbod steps.
//	noutlog	  : output logfile data once every nsteps/noutlog steps.
//	dtime	  : the timestep.
//	G	  : value of gravitational constant, in appropriate units.
//	abasis    : scale length of the basis functions.
//	outpcoef  : option to write-out expansion coefficients.
//	zeroodd   : option to zero all odd terms in the expansion.
//	zeroeven  : option to zero all even terms in the expansion.
//	fixacc    : option to force conservation of linear momentum
//	            by subtracting acceleration of c.o.m.
//	ubodlog   : option to write-out angular momentum and energy
//	            of individual stars.
//************************************************************************

        FILE *fp = fopen("SCFPAR","r");

	if(me == 0){
		if(!fp){
			printf("File not open\n");
			exit(1);
		}

		fscanf(fp,"%d", &(param->incont));
		fscanf(fp,"%d", &(param->initsnap));

		fscanf(fp,"%2s",param->basisset);

		fscanf(fp,"%d", &(param->nmax));
		fscanf(fp,"%d", &(param->lmax));
		fscanf(fp,"%d", &(param->nsteps));
		fscanf(fp,"%d", &(param->noutbod));
		fscanf(fp,"%d", &(param->noutlog));

		fscanf(fp,"%lf",&(param->dtime));
		fscanf(fp,"%lf",&(param->G));
		fscanf(fp,"%lf",&(param->abasis));

		fscanf(fp,"%d", &(param->outpcoef));
		fscanf(fp,"%d", &(param->zeroodd));
		fscanf(fp,"%d", &(param->zeroeven));
		fscanf(fp,"%d", &(param->fixacc));
		fscanf(fp,"%d", &(param->ubodlog));
		param->abasis2=(param->abasis)*(param->abasis);
		param->lskip=1;
		if(param->zeroodd == 1 || param->zeroeven == 1)
		 	param->lskip=2;

		param->lmin=0;
		if(param->zeroeven == 1)
			param->lmin=1;
	}
	MPI_Bcast(param,sizeof(struct Params),MPI_BYTE,0,MPI_COMM_WORLD);
}

void mpiget(double *rbuf,double *sbuf,int n,int idest,int isrc,int me)
{
	int itag;
	MPI_Status istatus;

	itag = isrc;

	if(me == isrc)
		MPI_Send(sbuf,n,MPI_DOUBLE,idest,itag,MPI_COMM_WORLD);
	if(me == idest)
		MPI_Recv(rbuf,n,MPI_DOUBLE,isrc,itag,MPI_COMM_WORLD,&istatus);
}

void writecoef(double *tnow,int nmax,int lmax,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],int me)
{
	static int firstc = 0;
	int l,m,n;
	FILE *fp;

	if(me == 0){
		if(firstc == 0){
			if((fp=fopen("SCFOCOEF","wx")) == NULL){
				printf("SCFOCOEF exists.\n");
				exit(1);
			}
		}else
			fp = fopen("SCFOCOEF","a");

		if(firstc == 0)
			firstc = 1;
		fprintf(fp,"%14.6E\n",*tnow);

		for(n=0;n<=nmax;n++)
			for(l=0;l<=lmax;l++)
				for(m=0;m<=l;m++){
					fprintf(fp,"%22.13E %22.13E\n",
					sinsum[n][l][m],cossum[n][l][m]);
				}
		fclose(fp);
	}
}

void outbods(int nbodies,int initsnap,double *tnow,double *mass,double *x,
		double *y,double *z,double *vx,double *vy,double *vz,
		double *pot,double *etotal,int ubodlog,int me,int n_pes)
{
	char a[50],b[50];
	static int firstc=0;
	static int isnap;
	int j,k,ntot;
	double lx,ly,lz,ek,ep,energy;
	double temp01[NMAX],temp02[NMAX],temp03[NMAX],temp04[NMAX],
		temp05[NMAX],temp06[NMAX],temp07[NMAX],temp08[NMAX];
	FILE *fp,*fpbod;

	if(firstc == 0){
		isnap=initsnap;
		firstc = 1;
	}

	if(me == 0){
		sprintf(a,"SNAP%03d",isnap);

		if((fp=fopen(a,"wx")) == NULL){
			printf("%s exists.\n",a);
			exit(1);
		}

		if(ubodlog == 1){
			sprintf(b,"SCFEL%03d",isnap);
			if((fpbod = fopen(b,"wx")) == NULL){
				printf("%s exists.\n",b);
				exit(1);
			}
		}
	}
	ntot = nbodies*n_pes;
	if(me == 0){
		fprintf(fp,"%9d %14.6E%14.6E\n",ntot,*tnow,*etotal);
		for(k=0;k<nbodies;k++)
			fprintf(fp,"%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",mass[k],x[k],y[k],z[k],vx[k],vy[k],vz[k],pot[k]);
		if(ubodlog == 1){
			fprintf(fpbod,"%9d %14.6E\n",ntot,*tnow);

			for(k=0;k<nbodies;k++){
				lx=mass[k]*(y[k]*vz[k]-z[k]*vy[k]);
				ly=mass[k]*(z[k]*vx[k]-z[k]*vz[k]);
				lz=mass[k]*(x[k]*vy[k]-z[k]*vx[k]);
				ek=0.5*mass[k]*(vx[k]*vx[k]+vy[k]*vy[k]
					+vz[k]*vz[k]);
				ep=mass[k]*pot[k];
				energy=ek+ep;
				fprintf(fpbod,"%14.6E%14.6E%14.6E%14.6E\n",
						lx,ly,lz,energy);
			}
		}
	}

	for(j=1;j<n_pes;j++){
		mpiget(temp01,mass,nbodies,0,j,me);
		mpiget(temp02,x   ,nbodies,0,j,me);
		mpiget(temp03,y   ,nbodies,0,j,me);
		mpiget(temp04,z   ,nbodies,0,j,me);
		mpiget(temp05,vx  ,nbodies,0,j,me);
		mpiget(temp06,vy  ,nbodies,0,j,me);
		mpiget(temp07,vz  ,nbodies,0,j,me);
		mpiget(temp08,pot ,nbodies,0,j,me);

                if(me == 0){
                        for(k=0;k<nbodies;k++){
				fprintf(fp,
			"%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",
					temp01[k],temp02[k],temp03[k],
					temp04[k],temp05[k],temp06[k],
					temp07[k],temp08[k]);
				if(ubodlog == 1){
					lx=temp01[k]*(temp03[k]*temp07[k]
						-temp04[k]*temp06[k]);
					ly=temp01[k]*(temp04[k]*temp05[k]
						-temp02[k]*temp07[k]);
					lz=temp01[k]*(temp02[k]*temp06[k]
						-temp03[k]*temp05[k]);
					ek=0.5*temp01[k]*(temp05[k]*temp05[k]
						+temp06[k]*temp06[k]
						+temp07[k]*temp07[k]);
					ep=0.5*temp01[k]*temp08[k];
					energy=ek+ep;
					fprintf(fpbod,"%14.6E%14.6E%14.6E%14.6E\n",
					lx,ly,lz,energy);
                                }
                        }
		}
	}

	if(me == 0){
		fclose(fp);
		if(ubodlog == 1)
			fclose(fpbod);
	}

	isnap++;
}

void outlog(struct Params *param,int nbodies,double *tnow,
		double *mass,double *x,double *y,double *z,double *vx,
		double *vy,double *vz,double *ax,double *ay,double *az,
		double *pot,double *etotal,double cputime0,int me,int n_pes)
{
	double cpux;
	static double cputime;
	static int firstc = 0;
	int k,nb=0,nbsum,nbsumall,nall;
	double mtot=0.0,xcm=0.0,ycm=0.0,zcm=0.0,vxcm=0.0,vycm=0.0,vzcm=0.0,
	       lxtot=0.0,lytot=0.0,lztot=0.0,etot=0.0,ektot=0.0,eptot=0.0,
	       epselfg=0.0,clausius=0.0,bclaus=0.0,vkin=0.0,vpot=0.0,
	       v2,potk,kink,ek,detot,m2tw,m2twb,m2claus,m2bclaus;
	double tmpsumsnd[17],tmpsumrcv[17];
	FILE *fp,*fpv;

	int incont = param->incont;
	int nmax = param->nmax;
	int lmax = param->lmax;

	cpux = MPI_Wtime();
	if(firstc == 0)
		cputime=cputime0;

	for(k=0;k<nbodies;k++){
		mtot += mass[k];
		xcm += mass[k]*x[k];
		ycm += mass[k]*y[k];
		zcm += mass[k]*z[k];
		vxcm += mass[k]*vx[k];
		vycm += mass[k]*vy[k];
		vzcm += mass[k]*vz[k];
		lxtot += mass[k]*(y[k]*vz[k]-z[k]*vy[k]);
		lytot += mass[k]*(z[k]*vx[k]-x[k]*vz[k]);
		lztot += mass[k]*(x[k]*vy[k]-y[k]*vx[k]);
		v2=vx[k]*vx[k]+vy[k]*vy[k]+vz[k]*vz[k];
		potk=0.5*mass[k]*pot[k];
		kink=0.5*mass[k]*v2;
		ek=0.5*v2+pot[k];
		if(ek < 0.0){
			vkin += kink;
			vpot += potk;
			bclaus += mass[k]*(x[k]*ax[k]+y[k]*ay[k]+z[k]*az[k]);
			nb += 1;
		}
		eptot += potk;
		epselfg += 0.5*mass[k]*pot[k];
		ektot += kink;
		clausius += mass[k]*(x[k]*ax[k]+y[k]*ay[k]+z[k]*az[k]);
	}

	tmpsumsnd[0]=mtot;
	tmpsumsnd[1]=xcm;
	tmpsumsnd[2]=ycm;
	tmpsumsnd[3]=zcm;
	tmpsumsnd[4]=vxcm;
	tmpsumsnd[5]=vycm;
	tmpsumsnd[6]=vzcm;
	tmpsumsnd[7]=lxtot;
	tmpsumsnd[8]=lytot;
	tmpsumsnd[9]=lztot;
	tmpsumsnd[10]=eptot;
	tmpsumsnd[11]=epselfg;
	tmpsumsnd[12]=ektot;
	tmpsumsnd[13]=vkin;
	tmpsumsnd[14]=vpot;
	tmpsumsnd[15]=clausius;
	tmpsumsnd[16]=bclaus;
	nbsum=nb;

	MPI_Allreduce(tmpsumsnd,tmpsumrcv,17,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&nbsum,&nbsumall,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	mtot    = tmpsumrcv[0];
	xcm     = tmpsumrcv[1];
	ycm     = tmpsumrcv[2];
	zcm     = tmpsumrcv[3];
	vxcm    = tmpsumrcv[4];
	vycm    = tmpsumrcv[5];
	vzcm    = tmpsumrcv[6];
	lxtot   = tmpsumrcv[7];
	lytot   = tmpsumrcv[8];
	lztot   = tmpsumrcv[9];
	eptot   = tmpsumrcv[10];
	epselfg = tmpsumrcv[11];
	ektot   = tmpsumrcv[12];
	vkin    = tmpsumrcv[13];
	vpot    = tmpsumrcv[14];
	clausius= tmpsumrcv[15];
	bclaus  = tmpsumrcv[16];

	xcm /= mtot;
	ycm /= mtot;
	zcm /= mtot;
	vxcm /= mtot;
	vycm /= mtot;
	vzcm /= mtot;

	etot=ektot+eptot;
	if(me == 0){
		if(firstc == 0 && incont == 0)
			*etotal=etot;
		detot=etot/(*etotal)-1.0;
		nbsum = nbsumall;
	}

	m2tw=-2.0*ektot/eptot;
	m2twb=-2.0*vkin/vpot;
	m2claus=-2.0*ektot/clausius;
	m2bclaus=-2.0*vkin/bclaus;
	nall = nbodies*n_pes;

	if(me == 0){
		if(firstc == 0){
			if((fp=fopen("SCFLOG","wx")) == NULL){
				printf("SCFLOG exists.\n");
				exit(1);
			}
		if((fpv=fopen("VIRIAL","w")) == NULL){
			printf("VIRIAL exists.\n");
			exit(1);
		}
		}else{
			fp=fopen("SCFLOG","a");
			fpv=fopen("VIRIAL","a");
		}
		if(*tnow == 0.0 || (incont != 0 && firstc == 0)){
			fprintf(fp,"%10d %10d\n",nmax,lmax);
			firstc = 1;
		}
		if(nbsum != nall)
			printf("nbound = %d\n",nbsum);

		fprintf(fp,"%18.10E\n",*tnow);
		fprintf(fp,"%18.10E%10d%10d\n",mtot,nall,nbsum);
		fprintf(fp,"%18.10E%18.10E%18.10E\n",xcm,ycm,zcm);
		fprintf(fp,"%18.10E%18.10E%18.10E\n",vxcm,vycm,vzcm);
		fprintf(fp,"%18.10E%18.10E%18.10E\n",lxtot,lytot,lztot);
		fprintf(fp,"%18.10E%18.10E%18.10E\n",ektot,eptot,epselfg);
		fprintf(fp,"%18.10E%18.10E%18.10E%18.10E\n",clausius,
				bclaus,etot,detot);
		fprintf(fp,"%18.10E%18.10E%18.10E%18.10E\n",m2tw,m2twb,
				m2claus,m2bclaus);
		fprintf(fp,"%18.10E\n",cpux-cputime);
//		fprintf(fp,"%18.10E\n",cpux-totime);
		fprintf(fp,"\n");
		fprintf(fpv,"%18.10E%18.10E%18.10E%18.10E%18.10E\n",
				*tnow,m2tw,m2twb,m2claus,m2bclaus);
		fclose(fp);
		fclose(fpv);
	}

	cputime=cpux;
}

void outstate(struct Params *param,int kstep,int nbodies,double *tvel,
		double *tnow,double *mass,double *x,double *y,double *z,
		double *vx,double *vy,double *vz,double *ax,double *ay,
		double *az,double *pot,double *etotal,double cputime,
		int me,int n_pes)
{
	int nlogmod,nbodmod;
	double rcsign;
	char message[18]=" step completed: ";

	int initsnap = param->initsnap;
	int noutbod = param->noutbod;
	int noutlog = param->noutlog;
	double dtime = param->dtime;
	int ubodlog = param->ubodlog;

	outterm(message,kstep,me);

	nlogmod=kstep-kstep/noutlog*noutlog;
	nbodmod=kstep-kstep/noutbod*noutbod;
	if(kstep == 0){
		outlog(param,nbodies,tnow,mass,x,y,z,vx,vy,vz,
				ax,ay,az,pot,etotal,cputime,me,n_pes);
		outbods(nbodies,initsnap,tnow,mass,x,y,z,vx,vy,vz,pot,
				etotal,ubodlog,me,n_pes);
	}else{
		if(nlogmod == 0 || nbodmod == 0){
			rcsign = -1.0;
			corrvel(rcsign,nbodies,tvel,tnow,dtime,vx,vy,vz,
				ax,ay,az);
			if(nlogmod == 0)
				outlog(param,nbodies,tnow,mass,x,y,z,
						vx,vy,vz,ax,ay,az,pot,
						etotal,cputime,me,n_pes);
			if(nbodmod == 0) outbods(nbodies,initsnap,tnow,
							mass,x,y,z,vx,vy,vz,
							pot,etotal,ubodlog,
							me,n_pes);
			rcsign = 1.0;
			corrvel(rcsign,nbodies,tvel,tnow,dtime,vx,vy,vz,
				ax,ay,az);
		}
	}
}

void outterm(char *message,int n,int me)
{
	FILE *fp;

	if(me == 0){
		fp = fopen("SCFOUT","a");
		if(n > 0){
			fprintf(fp,"%s %d\n",message,n);
			printf("%s %d\n",message,n);
		}

		fclose(fp);
	}
}

void scfcoef_CB(struct Params *param,int nbodies,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double c1[][LLMAX+1],double c2[][LLMAX+1],
		double *c3,double anltilde[][LLMAX+1],double *dblfact,
		double *twoalpha,double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1])
{
	int k,l,m,n,icount;
	double rxy2d,rxy2dinv,r2,r,costh,sinth,xi,cosmphi[LLMAX+1],
		sinmphi[LLMAX+1],ultrasp[NNMAX+1][LLMAX+1],
		ultraspt[NNMAX+1][LLMAX+1],un,unm1,ratio,a2r2,
		plm[LLMAX+1][LLMAX+1],plm1m,plm2m,temp3,temp4,temp5,ttemp5;

	int nmax = param->nmax;
	int lmax = param->lmax;
	int lmin = param->lmin;
	int lskip = param->lskip;
	double abasis = param->abasis;
	double abasis2 = param->abasis2;
	double tmpsumsnd[(NNMAX+1)*(LLMAX+1)*(LLMAX+2)];
	double tmpsumrcv[(NNMAX+1)*(LLMAX+1)*(LLMAX+2)];

	for(k=0;k<nbodies;k++){
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

			for(n=1;n<nmax;n++){
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
				tmpsumsnd[icount-2]=sinsum[n][l][m];
				tmpsumsnd[icount-1]=cossum[n][l][m];
			}
	MPI_Allreduce(tmpsumsnd,tmpsumrcv,icount,MPI_DOUBLE,MPI_SUM,
			MPI_COMM_WORLD);
	icount = 0;
	for(l=lmin;l<=lmax;l += lskip)
		for(m=0;m<=l;m++)
			for(n=0;n<=nmax;n++){
				icount += 2;
				sinsum[n][l][m]=tmpsumrcv[icount-2];
				cossum[n][l][m]=tmpsumrcv[icount-1];
			}
}

void scfcoef_LH(struct Params *param,int nbodies,
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
	double tmpsumsnd[(NNMAX+2)*(LLMAX+1)*(LLMAX+2)];
	double tmpsumrcv[(NNMAX+1)*(LLMAX+1)*(LLMAX+2)];

	for(k=0;k<nbodies;k++){
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

			for(n=1;n<nmax;n++){
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
}

void share(double *x,int n)
{
        MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


void shareint(int *x,int n)
{
        MPI_Bcast(x,n,MPI_INT,0,MPI_COMM_WORLD);
}

void startout(int me)
{
	FILE *fp;

	if(me == 0){
		if((fp = fopen("SCFOUT","wx"))==NULL){
			printf("SCFOUT exists.\n");
			exit(1);
		}

		fprintf(fp,"  Start of output\n");
		fclose(fp);
	}
}

void steppos(int nbodies,double *x,double *y,double *z,double *vx,
		double *vy,double *vz,double dtime,double *tnow,
		double *tpos)
{
	int k;

	for(k=0;k<nbodies;k++){
		x[k] += vx[k]*dtime;
		y[k] += vy[k]*dtime;
		z[k] += vz[k]*dtime;
	}

	*tpos += dtime;
	*tnow=*tpos;
}

void stepsys(struct Params *param,int kstep,int nbodies,
		double *tnow,double *tpos,double *tvel,
		double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double *vx,double *vy,double *vz,
		double *ax,double *ay,double *az,double *pot,
		double *etotal,double cputime,int me,int n_pes)
{
	double dtime = param->dtime;
	int fixacc = param->fixacc;

	steppos(nbodies,x,y,z,vx,vy,vz,dtime,tnow,tpos);

	accpot(param,nbodies,tnow,mass,x,y,z,ax,ay,az,pot,dblfact,twoalpha,
		c1,c2,c3,anltilde,lplusm,lmm,coeflm,me);

	if(fixacc == 1)
		corracc(nbodies,mass,ax,ay,az);

	stepvel(nbodies,vx,vy,vz,ax,ay,az,dtime,tnow,tvel);

	outstate(param,kstep,nbodies,tvel,tnow,mass,x,y,z,vx,vy,vz,
			ax,ay,az,pot,etotal,cputime,me,n_pes);
}

void stepvel(int nbodies,double *vx,double *vy,double *vz,double *ax,
		double *ay,double *az,double dtime,double *tnow,
		double *tvel)
{
	int k;

	for(k=0;k<nbodies;k++){
		vx[k] += ax[k]*dtime;
		vy[k] += ay[k]*dtime;
		vz[k] += az[k]*dtime;
	}

	*tvel += dtime;
	*tnow=*tvel;
}
