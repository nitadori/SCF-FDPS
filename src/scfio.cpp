#include <mpi.h>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>
#include "scfdhmpi.h"
#include "timer.hpp"

void endrun(double t01mpi, double t12mpi, int me)
{
	FILE *fp;

	if(me == 0){
		fp = fopen("SCFLOG","a");
		fprintf(fp, "     Initsys time        (secs) = %15.7E\n", t01mpi);
		fprintf(fp, "     Total cpu time used (secs) = %15.7E\n", t12mpi);

		fclose(fp);
	}
}

void initsys(struct Params *param,int ndisk,int nhalo,int nbodies,
		double *tnow,double *tpos,double *tvel,double *dblfact,
		double *twoalpha,double c1[][LLMAX+1],
		double c2[][LLMAX+1],double *c3,double anltilde[][LLMAX+1],
		double lplusm[][LLMAX+1],double lmm[][LLMAX+1],
		double coeflm[][LLMAX+1],
		double *mass,double *x,double *y, double *z,
		double *vx,double *vy,double *vz,
		double *ax,double *ay,double *az,double *pot, int *id,
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

        if(nhalo==0 && fixacc == 1)
                corracc(nbodies,mass,ax,ay,az);

	accpot(param,ndisk,nhalo,nbodies,tnow,mass,x,y,z,ax,ay,az,pot,
		dblfact,twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm);

        if(nhalo>0 && fixacc == 1)
                corracc(nbodies,mass,ax,ay,az);

	Timer::beg(Timer::IO_WRITE);
	outstate(param,kstep,ndisk,nhalo,nbodies,tvel,tnow,mass,x,y,z,
			vx,vy,vz,ax,ay,az,pot,id,cputime,etotal,me,n_pes);
	scfmpi_barr();
	Timer::end(Timer::IO_WRITE);

	initvel(nbodies,vx,vy,vz,ax,ay,az,dtime,tnow,tvel);
}

void inbods(int incont,int *ndisk,int *nhalo,double *tnow,double *mass,
		double *x,double *y,double *z,double *vx,double *vy,
		double *vz,int *id,double *etotal,int me,int n_pes,
		int skip_halo)
{
        int i,k,nb,nbodies;
        double tt,etot,mm,xx,yy,zz,vxx,vyy,vzz;
	//static double temp01[NMAX],temp02[NMAX],temp03[NMAX],temp04[NMAX],
	//	temp05[NMAX],temp06[NMAX],temp07[NMAX];
	double *base_buf=0, *temp01=0, *temp02=0, *temp03=0, *temp04=0, *temp05=0, *temp06=0, *temp07=0;
        FILE *fpd=0,*fph=0;

	if(me == 0){
		if((fph=fopen("SCFBI","r"))==NULL){
			printf("SCFBI not found.\n");
			exit(1);
		}

		if((fpd=fopen("TREEBI","r"))==NULL){
			printf("TREEBI not found.\n");
			exit(1);
		}

		if(incont == 0){
			fscanf(fph,"%d %lf",&nb,&tt);
			*nhalo=nb;
			*tnow=tt;
			printf("in inbods: nhalo = %d\n",*nhalo);

			if(skip_halo){
				*nhalo = 0;
				printf("in inbods: set nhalo = 0\n");
			}

			fscanf(fpd,"%d %lf",&nb,&tt);
			*ndisk=nb;
			printf("in inbods: ndisk = %d\n",*ndisk);
		}else{
			fscanf(fph,"%d %lf %lf",&nb,&tt,&etot);
			*nhalo=nb;
			// *tnow=tt;
			// *etotal=etot;

			fscanf(fpd,"%d %lf %lf",&nb,&tt,&etot);
			*ndisk=nb;
			*tnow=tt;
			*etotal=etot;
		}
		assert(base_buf = (double *)malloc(sizeof(double) * NMAX * 8));
	}
	nbodies = *ndisk + *nhalo;
	shareint(nhalo,1);
	shareint(ndisk,1);
	shareint(&nbodies,1);
	share(tnow,1);
//	printf("1 me = %d: nhalo = %d, ndisk = %d, tnow = %lf\n",me,*nhalo,
//		*ndisk,*tnow);

	assert(0 == *ndisk%n_pes);
	assert(0 == *nhalo%n_pes);
	*ndisk = *ndisk/n_pes;
	*nhalo = *nhalo/n_pes;
	nbodies = nbodies/n_pes;

//	printf("2 ndisk = %d  nhalo = %d  nbodies = %d\n",
//		*ndisk,*nhalo,nbodies);
	if(me == 0){
		for(k=0;k<*nhalo;k++){
			assert(k < NMAX);
			fscanf(fph,"%lf %lf %lf %lf %lf %lf %lf",
				&mm,&xx,&yy,&zz,&vxx,&vyy,&vzz);
			mass[k]=mm;
			x[k]=xx;
			y[k]=yy;
			z[k]=zz;
			vx[k]=vxx;
			vy[k]=vyy;
			vz[k]=vzz;
		}
		for(k=*nhalo;k<nbodies;k++){
			assert(k < NMAX);
			if(incont){
				int idk;
				fscanf(fpd, "%d", &idk);
				id[k] = idk;
			}
			fscanf(fpd,"%lf %lf %lf %lf %lf %lf %lf",
				&mm,&xx,&yy,&zz,&vxx,&vyy,&vzz);
			if(incont){
				double potk;
				fscanf(fpd, "%lf", &potk);
			}
			mass[k]=mm;
			x[k]=xx;
			y[k]=yy;
			z[k]=zz;
			vx[k]=vxx;
			vy[k]=vyy;
			vz[k]=vzz;
		}
	}
//		for(k=0;k<nbodies;k++)
//			printf("%12.3E%12.3E%12.3E%12.3E%12.3E%12.3E%12.3E\n",
//				mass[k],x[k],y[k],z[k],vx[k],vy[k],vz[k]);
	for(i=1;i<n_pes;i++){
		int *itemp01=0;
		if(me == 0){
			temp01 = base_buf;
			temp02 = temp01 + NMAX;
			temp03 = temp02 + NMAX;
			temp04 = temp03 + NMAX;
			temp05 = temp04 + NMAX;
			temp06 = temp05 + NMAX;
			temp07 = temp06 + NMAX;
			itemp01 = (int *)(temp07 + NMAX);

			for(k=0;k<*nhalo;k++){
				assert(k < NMAX);
				fscanf(fph,"%lf %lf %lf %lf %lf %lf %lf",
					&mm,&xx,&yy,&zz,&vxx,&vyy,&vzz);
				temp01[k]=mm;
				temp02[k]=xx;
				temp03[k]=yy;
				temp04[k]=zz;
				temp05[k]=vxx;
				temp06[k]=vyy;
				temp07[k]=vzz;
			}
//			for(k=0;k<*nhalo;k++)
//				printf("%d%12.3E%12.3E%12.3E%12.3E%12.3E%12.3E%12.3E\n",k,temp01[k],temp02[k],temp03[k],temp04[k],temp05[k],temp06[k],
//				temp07[k]);
			for(k=*nhalo;k<nbodies;k++){
				assert(k < NMAX);
				if(incont){
					int idk;
					fscanf(fpd, "%d", &idk);
					itemp01[k] = idk;
				}
				fscanf(fpd,"%lf %lf %lf %lf %lf %lf %lf",
					&mm,&xx,&yy,&zz,&vxx,&vyy,&vzz);
				if(incont){
					double potk;
					fscanf(fpd, "%lf", &potk);
				}
				temp01[k]=mm;
				temp02[k]=xx;
				temp03[k]=yy;
				temp04[k]=zz;
				temp05[k]=vxx;
				temp06[k]=vyy;
				temp07[k]=vzz;
			}
//			for(k=*nhalo;k<*nhalo+*ndisk;k++)
//				printf("%d%12.3E%12.3E%12.3E%12.3E%12.3E%12.3E%12.3E\n",k,temp01[k],temp02[k],temp03[k],temp04[k],temp05[k],temp06[k],
//				temp07[k]);
		}
	
		mpiget(mass,temp01,nbodies,i,0,me);
		mpiget(x   ,temp02,nbodies,i,0,me);
		mpiget(y   ,temp03,nbodies,i,0,me);
		mpiget(z   ,temp04,nbodies,i,0,me);
		mpiget(vx  ,temp05,nbodies,i,0,me);
		mpiget(vy  ,temp06,nbodies,i,0,me);
		mpiget(vz  ,temp07,nbodies,i,0,me);
		if(incont){
			mpiiget(id,itemp01,nbodies,i,0,me);
		}
	}

	if(me == 0){
		fclose(fph);
		fclose(fpd);
		free(base_buf);
	}
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
//	zeroodd   : option to zero all odd terms in the expansion.
//	zeroeven  : option to zero all even terms in the expansion.
//	fixacc    : option to force conservation of linear momentum
//	            by subtracting acceleration of c.o.m.
//	ubodlog   : option to write-out angular momentum and energy
//	            of individual stars.
//      theta     : opening angle for tree
//      n_leaf_limit:  maximum # of particles in a leaf cell
//      n_group_limit: maximum # of particles in i-particle group
//      : maximum # of particles in i-particle group
//      cnfw      : concentration parameter of the halo
//      rmax      : maximum cut-off radius of the halo
//      mh        : halo mass
//************************************************************************


	if(me == 0){
		FILE *fp = fopen("SCFPAR","r");
        	if(!fp){
			printf("File not open\n");
			scfmpi_abort(1);
        	}

// from https://marycore.jp/prog/c-lang/scanf-string-safely/
#define TAIL "%*[^\n]%*c"
		fscanf(fp,"%d"  TAIL, &(param->incont));
		fscanf(fp,"%d"  TAIL, &(param->initsnap));

		fscanf(fp,"%2s" TAIL, param->basisset);

		fscanf(fp,"%d"  TAIL, &(param->nmax));
		fscanf(fp,"%d"  TAIL, &(param->lmax));
		fscanf(fp,"%d"  TAIL, &(param->nsteps));
		fscanf(fp,"%d"  TAIL, &(param->noutbod));
		fscanf(fp,"%d"  TAIL, &(param->noutlog));

		fscanf(fp,"%lf" TAIL, &(param->dtime));
		fscanf(fp,"%lf" TAIL, &(param->G));
		fscanf(fp,"%lf" TAIL, &(param->eps));
		fscanf(fp,"%lf" TAIL, &(param->abasis));

		fscanf(fp,"%d"  TAIL, &(param->zeroodd));
		fscanf(fp,"%d"  TAIL, &(param->zeroeven));
		fscanf(fp,"%d"  TAIL, &(param->fixacc));
		fscanf(fp,"%d"  TAIL, &(param->ubodlog));



		int ret = fscanf(fp,"%d" TAIL, &(param->outpcoef));
		assert(1 == ret);
		if(param->outpcoef){
			fprintf(stderr, "Experimental: outpcoef=%d", param->outpcoef);
		}

		param->abasis2=(param->abasis)*(param->abasis);
		param->lskip=1;

                fscanf(fp,"%3s" TAIL, param->softening_type);
                fscanf(fp,"%d"  TAIL, &(param->mm_order));
                fscanf(fp,"%lf" TAIL, &(param->theta));
                fscanf(fp,"%d"  TAIL, &(param->n_leaf_limit));
                fscanf(fp,"%d"  TAIL, &(param->n_group_limit));
                assert( (0 == strcmp(param->softening_type,"PLM"))
                     || (0 == strcmp(param->softening_type,"SPL")) );
                assert( param->mm_order <= 2 );

                // std::cerr<<"param->n_group_limit= "<<param->n_group_limit<<std::endl;
		fprintf(stderr, "param->n_group_limit= %d\n", param->n_group_limit);
		if(0 == strcmp(param->basisset,"EX")){
			int nret = 0;
			nret += fscanf(fp,"%lf" TAIL, &(param->cnfw));
			nret += fscanf(fp,"%lf" TAIL, &(param->rmax));
			nret += fscanf(fp,"%lf" TAIL, &(param->mh));
			assert(3 == nret);
		}else if(0 == strcmp(param->basisset,"CB")){
			// nothing to do
		}else if(0 == strcmp(param->basisset,"LH")){
			// nothing to do
		}else if(0 == strcmp(param->basisset,"NA")){
			// nothing to do
		}else{
			fprintf(stderr, "Invalid basisset name!\n");
			scfmpi_abort(-1);
		}

		if(param->zeroodd == 1 || param->zeroeven == 1)
			param->lskip=2;

		param->lmin=0;
		if(param->zeroeven == 1)
			param->lmin=1;

		fclose(fp);
	}

	MPI_Bcast(param,sizeof(struct Params),MPI_BYTE,0,MPI_COMM_WORLD);
#undef TAIL
}

static_assert(8 == sizeof(double));
static_assert(8 == sizeof(int64_t));

void outbods(int ndisk,int nhalo,int nbodies,int initsnap,double *tnow,
		double *mass,double *x,double *y,double *z,double *vx,
		double *vy,double *vz,double *pot, int *id, double *etotal,int ubodlog,
		int me,int n_pes)
{
	static char a[50],b[50],c[50],d[50];
	static int firstc=0;
	static int isnap;
	int j,k,ndtot,nhtot;
	double lx,ly,lz,ek,ep,energy;
	//static double temp01[NMAX],temp02[NMAX],temp03[NMAX],temp04[NMAX],
	//	temp05[NMAX],temp06[NMAX],temp07[NMAX],temp08[NMAX];
	double *base_buf=0, 
	       *temp01=0, *temp02=0, *temp03=0, *temp04=0, 
	       *temp05=0, *temp06=0, *temp07=0, *temp08=0;;
	int *itemp01=0;
	FILE *fpd=0, *fph=0, *fpdbod=0, *fphbod=0;

	if(firstc == 0){
		isnap=initsnap;
		firstc = 1;
	}

	if(me == 0){
		sprintf(a,"DSNAP%03d",isnap);
		sprintf(b,"HSNAP%03d",isnap);

#if !defined(NO_EXCLUSIVE_IO)
		if(((fpd=fopen(a,"wx")) == NULL) ||
			((nhalo>0 ) && ((fph=fopen(b,"wx")) == NULL)))
		{
				printf("SNAP file exists.\n");
				exit(1);
		}
#else
                fpd=fopen(a,"w");
                fph=fopen(b,"w");
#endif
		if(ubodlog == 1){
			sprintf(c,"DSCFEL%03d",isnap);
			sprintf(d,"HSCFEL%03d",isnap);

#if !defined(NO_EXCLUSIVE_IO)
			if(((fpdbod=fopen(c,"wx")) == NULL) ||
				((fphbod=fopen(d,"wx")) == NULL)){
				printf("SCFEL file exists.\n");
				exit(1);
			}
#else
                        fpdbod=fopen(c,"w");
                        fphbod=fopen(d,"w");
#endif
		}
		assert(base_buf = (double *)malloc(sizeof(double) * NMAX * 9));
	}

	//ndtot = ndisk*n_pes;
	MPI_Allreduce(&ndisk, &ndtot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	//nhtot = nhalo*n_pes;
	MPI_Allreduce(&nhalo, &nhtot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        int ndisks[n_pes];
        int nhalos[n_pes];
        int nbodiess[n_pes];
        MPI_Allgather(&ndisk, 1, MPI_INT, ndisks, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&nhalo, 1, MPI_INT, nhalos, 1, MPI_INT, MPI_COMM_WORLD);
        for(k=0; k<n_pes; k++){
            nbodiess[k] = ndisks[k] + nhalos[k];
        }
        
	if(me == 0){
#ifndef OUTBODS_BINARY
		(nhalo>0) && fprintf(fph,"%9d %14.6E %14.6E\n",nhtot,*tnow,*etotal);
		fprintf(fpd,"%9d %14.6E %14.6E\n",ndtot,*tnow,*etotal);
#else
		if(nhalo > 0){
			int64_t nn = nhtot;
			fwrite(&nn,    sizeof(int64_t), 1, fph);
			fwrite(tnow,   sizeof(double),  1, fph);
			fwrite(etotal, sizeof(double),  1, fph);
		}
		{
			int64_t nn = ndtot;
			fwrite(&nn,    sizeof(int64_t), 1, fpd);
			fwrite(tnow,   sizeof(double),  1, fpd);
			fwrite(etotal, sizeof(double),  1, fpd);
		}
#endif
	}
	
	if(me == 0){
		for(k=0;k<nhalo;k++)
#ifndef OUTBODS_BINARY
			(nhalo>0) && fprintf(fph,"%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",mass[k],x[k],y[k],z[k],vx[k],vy[k],vz[k],pot[k]);
#else
		if(nhalo > 0){
			double buf[8] = {mass[k], x[k], y[k], z[k], vx[k], vy[k], vz[k], pot[k]};
			fwrite(buf, sizeof(double), 8, fph);
		}
#endif

		if(ubodlog == 1){
			fprintf(fphbod,"%9d %14.6E\n",nhalo,*tnow);

			for(k=0;k<nhalo;k++){
				lx=mass[k]*(y[k]*vz[k]-z[k]*vy[k]);
				ly=mass[k]*(z[k]*vx[k]-z[k]*vz[k]);
				lz=mass[k]*(x[k]*vy[k]-z[k]*vx[k]);
				ek=0.5*mass[k]*(vx[k]*vx[k]+vy[k]*vy[k]
					+vz[k]*vz[k]);
				ep=0.5*mass[k]*pot[k];
				energy=ek+ep;
				fprintf(fphbod,"%14.6E%14.6E%14.6E%14.6E\n",
					lx,ly,lz,energy);
			}
		}
		for(k=nhalo;k<nbodies;k++){
#ifndef OUTBODS_BINARY
			fprintf(fpd,"%8d%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",id[k],mass[k],x[k],y[k],z[k],vx[k],vy[k],vz[k],pot[k]);
#else
			double buf[9] = {0.0, mass[k], x[k], y[k], z[k], vx[k], vy[k], vz[k], pot[k]};
			// *(int64_t *)buf = id[k];
			int64_t idk = id[k];
			memcpy(buf, &idk, sizeof(int64_t));
			fwrite(buf, sizeof(double), 9, fpd);
#endif
		}

		if(ubodlog == 1){
			fprintf(fpdbod,"%9d %14.6E\n",ndisk,*tnow);

			for(k=nhalo;k<nbodies;k++){
				lx=mass[k]*(y[k]*vz[k]-z[k]*vy[k]);
				ly=mass[k]*(z[k]*vx[k]-z[k]*vz[k]);
				lz=mass[k]*(x[k]*vy[k]-z[k]*vx[k]);
				ek=0.5*mass[k]*(vx[k]*vx[k]+vy[k]*vy[k]
					+vz[k]*vz[k]);
				ep=0.5*mass[k]*pot[k];
				energy=ek+ep;
				fprintf(fpdbod,"%14.6E%14.6E%14.6E%14.6E\n",
					lx,ly,lz,energy);
			}
		}
	}

	if(me == 0){
		temp01 = base_buf;
		temp02 = temp01 + NMAX;
		temp03 = temp02 + NMAX;
		temp04 = temp03 + NMAX;
		temp05 = temp04 + NMAX;
		temp06 = temp05 + NMAX;
		temp07 = temp06 + NMAX;
		temp08 = temp07 + NMAX;
		itemp01 = (int *)(temp08 + NMAX);
	}

	for(j=1;j<n_pes;j++){
            int nb_tmp = nbodiess[j];
            int nh_tmp = nhalos[j];
            /*
		mpiget(temp01,mass,nbodies,0,j,me);
		mpiget(temp02,x   ,nbodies,0,j,me);
		mpiget(temp03,y   ,nbodies,0,j,me);
		mpiget(temp04,z   ,nbodies,0,j,me);
		mpiget(temp05,vx  ,nbodies,0,j,me);
		mpiget(temp06,vy  ,nbodies,0,j,me);
		mpiget(temp07,vz  ,nbodies,0,j,me);
		mpiget(temp08,pot ,nbodies,0,j,me);
            */
		mpiget(temp01,mass,nb_tmp,0,j,me);
		mpiget(temp02,x   ,nb_tmp,0,j,me);
		mpiget(temp03,y   ,nb_tmp,0,j,me);
		mpiget(temp04,z   ,nb_tmp,0,j,me);
		mpiget(temp05,vx  ,nb_tmp,0,j,me);
		mpiget(temp06,vy  ,nb_tmp,0,j,me);
		mpiget(temp07,vz  ,nb_tmp,0,j,me);
		mpiget(temp08,pot ,nb_tmp,0,j,me);
		mpiiget(itemp01, id ,nb_tmp,0, j, me);
		if(me == 0){
                    //for(k=0;k<nhalo;k++){
                    for(k=0;k<nh_tmp;k++){
#ifndef OUTBODS_BINARY
				(nhalo>0) && fprintf(fph,
			"%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",
					temp01[k],temp02[k],temp03[k],
					temp04[k],temp05[k],temp06[k],
					temp07[k],temp08[k]);
#else
				if(nhalo > 0){
					double buf[8] = {temp01[k],temp02[k],temp03[k],
					                 temp04[k],temp05[k],temp06[k],
					                 temp07[k],temp08[k]};
					fwrite(buf, sizeof(double), 8, fph);
				}
#endif
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
					fprintf(fphbod,"%14.6E%14.6E%14.6E%14.6E\n",
						lx,ly,lz,energy);
				}
			}

                    //for(k=nhalo;k<nbodies;k++){
                    for(k=nh_tmp;k<nb_tmp;k++){
#ifndef OUTBODS_BINARY
				fprintf(fpd,
			"%8d%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",
					itemp01[k],
					temp01[k],temp02[k],temp03[k],
					temp04[k],temp05[k],temp06[k],
					temp07[k],temp08[k]);
#else
				double buf[9] = {0.0,       temp01[k],temp02[k],
				                 temp03[k], temp04[k],temp05[k],
						 temp06[k], temp07[k],temp08[k]};
				int64_t idk = itemp01[k];
				memcpy(buf, &idk, sizeof(int64_t));
				fwrite(buf, sizeof(double), 9, fpd);
#endif
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
					fprintf(fpdbod,"%14.6E%14.6E%14.6E%14.6E\n",
						lx,ly,lz,energy);
				}
			}
		}
	}

	if(me == 0){
		(nhalo>0) && fclose(fph);
		fclose(fpd);
		if(ubodlog == 1){
			fclose(fphbod);
			fclose(fpdbod);
		}
		free(base_buf);
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
	eptot   += 0.5 * Common::ext_pot_energy;
	epselfg -= 0.5 * Common::ext_pot_energy;

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
	//nall = nbodies*n_pes;
	MPI_Allreduce(&nbodies, &nall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if(me == 0){
            if(firstc == 0){
#if !defined(NO_EXCLUSIVE_IO)
                if((fp=fopen("SCFLOG","wx")) == NULL){
                    printf("SCFLOG exists.\n");
                    exit(1);
                }
#else
                fp=fopen("SCFLOG","w");
#endif
#if !defined(NO_EXCLUSIVE_IO)
		if((fpv=fopen("VIRIAL","wx")) == NULL){
                    printf("VIRIAL exists.\n");
                    exit(1);
		}
#else
                fpv=fopen("VIRIAL","w");
#endif                
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
            fprintf(fp,"%18.10E%18.10E%18.10E\n",ektot, eptot, epselfg);
            fprintf(fp,"%18.10E%18.10E%18.10E%18.10E\n",clausius, bclaus, etot, detot);
            fprintf(fp,"%18.10E%18.10E%18.10E%18.10E\n",m2tw, m2twb, m2claus, m2bclaus);
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

void outstate(struct Params *param,int kstep,int ndisk,int nhalo,
		int nbodies,double *tvel,double *tnow,double *mass,
		double *x,double *y,double *z,
		double *vx,double *vy,double *vz,double *ax,double *ay,
		double *az,double *pot, int *id, double cputime,double *etotal,
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
#if 1                
		outbods(ndisk,nhalo,nbodies,initsnap,tnow,mass,x,y,z,
				vx,vy,vz,pot,id,etotal,ubodlog,me,n_pes);
#endif
	}else{
		if(nlogmod == 0 || nbodmod == 0){
			rcsign = -1.0;
			corrvel(rcsign,nbodies,tvel,tnow,dtime,vx,vy,vz,
				ax,ay,az);
			if(nlogmod == 0)
				outlog(param,nbodies,tnow,mass,x,y,z,
						vx,vy,vz,ax,ay,az,pot,
						etotal,cputime,me,n_pes);
#if 1
			if(nbodmod == 0) outbods(ndisk,nhalo,nbodies,initsnap,
							tnow,mass,x,y,z,vx,vy,
							vz,pot,id,etotal,ubodlog,
							me,n_pes);
#endif
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

void startout(int me)
{
	FILE *fp;

	if(me == 0){
#if !defined(NO_EXCLUSIVE_IO)
		if((fp = fopen("SCFOUT","wx"))==NULL){
			printf("SCFOUT exists.\n");
			exit(1);
		}
#else
                fp = fopen("SCFOUT","w");
#endif
		fprintf(fp,"  Start of output\n");
		fclose(fp);
	}
}

void outputcoef(int nmax,int lmax,double *tnow,
                double sinsumhh[][LLMAX+1][LLMAX+1],
                double cossumhh[][LLMAX+1][LLMAX+1])
{
        static int firstc = 0;
        int l,m,n;
        FILE *fp;

	int me = scfmpi_rank();
	if(me) return;

        if(firstc == 0){
                firstc = 1;
                if((fp=fopen("SCFOCOEF","wx")) == NULL){
                        printf("SCFOCOEF exists.\n");
                        exit(1);
                }
        }else{
                fp=fopen("SCFOCOEF","a");
        }

        fprintf(fp,"%14.6E\n",*tnow);
        for(n=0;n<=nmax;n++)
                for(l=0;l<=lmax;l++)
                        for(m=0;m<=l;m++)
                                fprintf(fp,"%22.13E %22.13E\n",
                                                sinsumhh[n][l][m],
                                                cossumhh[n][l][m]);
        fclose(fp);
}
