// #include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include "scfdhmpi.h"

#include "timer.hpp"

#define scfcoef_CB scfcoef_CB_avx
#define accp_CB    accp_CB_avx
#define scfcoef_LH scfcoef_LH_avx
#define accp_LH    accp_LH_avx

void accpot(struct Params *param,int ndisk,int nhalo,int nbodies,
		double *tnow,double *mass,double *x,double *y,double *z,
		double *ax,double *ay,double *az,double *pot,
		double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1])
{
	int k,l,m,n;
	static int nzero=0;

	const char *basisset = param->basisset;
	int nmax = param->nmax;
	int lmax = param->lmax;
	double sinsumhd[NNMAX+1][LLMAX+1][LLMAX+1];
	double cossumhd[NNMAX+1][LLMAX+1][LLMAX+1];
	double sinsumhh[NNMAX+1][LLMAX+1][LLMAX+1];
	double cossumhh[NNMAX+1][LLMAX+1][LLMAX+1];

	if(nhalo > 0){
		Timer::beg(Timer::SCF_COEF);
		for(l=0;l<=lmax;l++)
			for(m=0;m<=l;m++)
				for(n=0;n<=nmax;n++){
					sinsumhd[n][l][m]=0.0;
					cossumhd[n][l][m]=0.0;
					sinsumhh[n][l][m]=0.0;
					cossumhh[n][l][m]=0.0;
				}

		if(strcmp(basisset,"CB") == 0){
			scfcoef_CB(param,nzero,nhalo,sinsumhd,cossumhd,mass,x,y,z,
					c1,c2,c3,anltilde,dblfact,twoalpha,lplusm,lmm,
					coeflm);
			scfcoef_CB(param,nhalo,nbodies,sinsumhh,cossumhh,mass,x,y,z,
					c1,c2,c3,anltilde,dblfact,twoalpha,lplusm,lmm,
					coeflm);
			for(l=0;l<=lmax;l++)
				for(m=0;m<=l;m++)
					for(n=0;n<=nmax;n++){
						sinsumhh[n][l][m] += sinsumhd[n][l][m];
						cossumhh[n][l][m] += cossumhd[n][l][m];
					}
		}else if(strcmp(basisset,"LH")==0){
			scfcoef_LH(param,nzero,nhalo,sinsumhd,cossumhd,mass,x,y,z,
					c1,c2,c3,anltilde,dblfact,twoalpha,lplusm,lmm,
					coeflm);
			scfcoef_LH(param,nhalo,nbodies,sinsumhh,cossumhh,mass,x,y,z,
					c1,c2,c3,anltilde,dblfact,twoalpha,lplusm,lmm,
					coeflm);
			for(l=0;l<=lmax;l++)
				for(m=0;m<=l;m++)
					for(n=0;n<=nmax;n++){
						sinsumhh[n][l][m] += sinsumhd[n][l][m];
						cossumhh[n][l][m] += cossumhd[n][l][m];
					}
		}else if(strcmp(basisset,"EX")==0){
			// nothing to do
		}else if(strcmp(basisset,"NA")==0){
			// nothing to do
		}else{
			assert(0);
		}
		Timer::end(Timer::SCF_COEF);
	}

	Timer::beg(Timer::SCF_ACCP);
	for(k=0;k<nhalo;k++){
		ax[k]=0.0;
		ay[k]=0.0;
		az[k]=0.0;
		pot[k]=0.0;
	}

	if(param->outpcoef != 0){
		outputcoef(nmax,lmax,tnow,sinsumhh,cossumhh);
	}

        //	accp_tree(param,ndisk,mass+nhalo,x+nhalo,y+nhalo,z+nhalo,
        //			ax+nhalo,ay+nhalo,az+nhalo,pot+nhalo);

        //        std::cerr<<"ax[nhalo]= "<<ax[nhalo]<<" pot[nhalo]= "<<pot[nhalo]<<std::endl;
        
	if(strcmp(basisset,"CB") == 0){
		accp_CB(param,nhalo,nbodies,mass,x,y,z,ax,ay,az,pot,dblfact,
			twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm,
			sinsumhd,cossumhd);
		accp_CB(param,nzero,nhalo,mass,x,y,z,ax,ay,az,pot,dblfact,
			twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm,
			sinsumhh,cossumhh);
	}else if(strcmp(basisset,"LH") == 0){
		accp_LH(param,nhalo,nbodies,mass,x,y,z,ax,ay,az,pot,dblfact,
			twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm,
			sinsumhd,cossumhd);
		accp_LH(param,nzero,nhalo,mass,x,y,z,ax,ay,az,pot,dblfact,
			twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm,
			sinsumhh,cossumhh);
	}else if(strcmp(basisset,"EX")==0){
		assert(0 == nhalo);
		assert(ndisk == nbodies);

		const double rs     = param->rmax / param-> cnfw;
		const double cplus1 = 1.0 + param->cnfw;
		const double cpot   = 1.0/cplus1;
		const double cf     = param->cnfw / (log(cplus1) - param->cnfw*cpot);
		const double Gpot   = cf*(param->G)*(param->mh) / param->rmax;
                const double Gforce = Gpot*(param->cnfw) / param->rmax;

		double e_ext = 0.0;
		for(k=0;k<ndisk;k++){
			double r2     = x[k]*x[k]+y[k]*y[k]+z[k]*z[k];
			double r      = sqrt(r2);
			double rnorm  = r/rs;
			double rplus1 = rnorm+1.0;
			double rpot   = log(rplus1)/rnorm;
			double rforce = Gforce*(rpot-1.0/rplus1)/rnorm;
			ax[k]  -= rforce/r * x[k];
			ay[k]  -= rforce/r * y[k];
			az[k]  -= rforce/r * z[k];
			pot[k] -= Gpot*rpot;
			e_ext -= mass[k] * Gpot * rpot;
		}
		Common::ext_pot_energy = e_ext;
	}else if(strcmp(basisset,"NA")==0){
			// nothing to do
	}else{
		assert(0);
	}
	Timer::end(Timer::SCF_ACCP);
}

#if 0
void accp_tree(struct Params *param,int ndisk,
	       double *mass,double *x,double *y,double *z,
		double *ax,double *ay,double *az,double *pot)
{
	static int firstc=0;
	int i,j;
	double sx,sy,sz,spot,xij,yij,zij,rij2,rij,cacc,cpot;
	double G = param->G;
	double eps = param->eps;
	static double eps2;
	double tmpsumsnd[4],tmpsumrcv[4];

	if(firstc == 0){
		eps2=eps*eps;
		firstc = 1;
	}

	for(i=0;i<ndisk;i++){
		sx = 0.0;
		sy = 0.0;
		sz = 0.0;
		spot = 0.0;
		for(j=0;j<ndisk;j++){
			if(i != j){
				xij=x[j]-x[i];
				yij=y[j]-y[i];
				zij=z[j]-z[i];
				rij2 = xij*xij+yij*yij+zij*zij+eps2;
				rij = sqrt(rij2);
				cpot=mass[j]/rij;
				cacc=cpot/rij2;
				sx += cacc*xij;
				sy += cacc*yij;
				sz += cacc*zij;
				spot += cpot;
			}
		}
		tmpsumsnd[0]=sx;
		tmpsumsnd[1]=sy;
		tmpsumsnd[2]=sz;
		tmpsumsnd[3]=cpot;
		MPI_Allreduce(tmpsumsnd,tmpsumrcv,4,MPI_DOUBLE,
				MPI_SUM,MPI_COMM_WORLD);
		sx = tmpsumrcv[0];
		sy = tmpsumrcv[1];
		sz = tmpsumrcv[2];
		cpot = tmpsumrcv[3];

		ax[i] = G*sx;
		ay[i] = G*sy;
		az[i] = G*sz;
		pot[i] = -G*spot;
	}

}
#endif

void corracc(int nbodies,double *mass,double *ax,double *ay,double *az)
{
	int k;
	double mtot=0.0,axcm=0.0,aycm=0.0,azcm=0.0;
	// double tmpsumsnd[4],tmpsumrcv[4];
	double tmpsum[4];

	for(k=0;k<nbodies;k++){
		mtot += mass[k];
		axcm += mass[k]*ax[k];
		aycm += mass[k]*ay[k];
		azcm += mass[k]*az[k];
	}

	tmpsum[0]=mtot;
	tmpsum[1]=axcm;
	tmpsum[2]=aycm;
	tmpsum[3]=azcm;

	// MPI_Allreduce(tmpsumsnd,tmpsumrcv,4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	scfmpi_sum_double(tmpsum, 4);

	mtot = tmpsum[0];
	axcm = tmpsum[1];
	aycm = tmpsum[2];
	azcm = tmpsum[3];

	axcm /= mtot;
	aycm /= mtot;
	azcm /= mtot;

	for(k=0;k<nbodies;k++){
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

void initpars(double *tpos,double *tvel,double *tnow)
{
	*tpos=*tnow;
	*tvel=*tnow;
}


void initvel(
		int nbodies,
		double * __restrict vx, double * __restrict vy, double * __restrict vz,
		double * __restrict ax, double * __restrict ay, double * __restrict az,
		double dtime,double *tnow,double *tvel)
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

void steppos(
		int nbodies,
		double * __restrict x,  double * __restrict y,  double * __restrict z,
		double * __restrict vx, double * __restrict vy, double * __restrict vz,
		double dtime, double *tnow, double *tpos)
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

void stepsys(struct Params *param,int kstep,int ndisk,int nhalo,int nbodies,
                double *tnow,double *tpos,double *tvel,
                double *dblfact,double *twoalpha,
                double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
                double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
                double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
                double *mass,double *x,double *y,double *z,
                double *vx,double *vy,double *vz,
                double *ax,double *ay,double *az,double *pot, int *id,
                double *etotal,double cputime,int me,int n_pes)
{
        double dtime = param->dtime;
        int fixacc = param->fixacc;

        //steppos(nbodies,x,y,z,vx,vy,vz,dtime,tnow,tpos);

	// correct disk acc befor evaluating external potential
        if(nhalo==0 && fixacc == 1)
                corracc(nbodies,mass,ax,ay,az);

        accpot(param,ndisk,nhalo,nbodies,tnow,mass,x,y,z,ax,ay,az,pot,
                dblfact,twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm);

        if(nhalo>0 && fixacc == 1)
                corracc(nbodies,mass,ax,ay,az);

	Timer::beg(Timer::STEP_VEL);
        stepvel(nbodies,vx,vy,vz,ax,ay,az,dtime,tnow,tvel);
	Timer::end(Timer::STEP_VEL);
#if 1
	Timer::beg(Timer::IO_WRITE);
        outstate(param,kstep,ndisk,nhalo,nbodies,tvel,tnow,mass,x,y,z,
                        vx,vy,vz,ax,ay,az,pot,id,cputime,etotal,me,n_pes);
	scfmpi_barr();
	Timer::end(Timer::IO_WRITE);
#endif
}

void stepvel(
		int nbodies,
		double * __restrict vx, double * __restrict vy, double * __restrict vz,
		double * __restrict ax, double * __restrict ay, double * __restrict az,
		double dtime,double *tnow,
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
