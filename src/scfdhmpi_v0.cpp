#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include "user_defined.hpp"
#include "specializations.hpp" // SIMD kernels for FDPS
#include <particle_simulator.hpp>
//#include <omp.h>
//omp_get_wtime();

#include "scfdhmpi.h"

#include "timer.hpp"


void TreeProf::show_local(FILE *fp=stdout){
	double eplen = niepjsum / nisum;
	double splen = nispjsum / nisum;
	fprintf(fp, "local:  <#EPJ>=%f, <#SPJ>=%f\n", eplen, splen);
	// nisum = niepjsum = nispjsum = 0.0;
}

void TreeProf::show_global(FILE *fp=stdout){
	double buf[3] = {nisum, niepjsum, nispjsum};
	scfmpi_sum_double(buf, 3);
	double eplen = buf[1] / buf[0];
	double splen = buf[2] / buf[0];
	if(!scfmpi_rank()){
		fprintf(fp, "global: <#EPJ>=%f, <#SPJ>=%f\n", eplen, splen);
	}
	// nisum = niepjsum = nispjsum = 0.0;
}
 
void print_average_etc(double val){
	int n = scfmpi_size();

	double avr = val;
	scfmpi_sum_double(&avr, 1);
	avr /= n;

	double var = val - avr;
	var *= var;
	scfmpi_sum_double(&var, 1);
	var = sqrt(var / n);

	double maxmin[] = {val, -val};
	MPI_Allreduce(MPI_IN_PLACE, maxmin, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	double max = +maxmin[0];
	double min = -maxmin[1];

	if(0 == scfmpi_rank()) {
		printf("avr: %e, var: %e, max %e, min %e\n", avr, var, max, min);
	}
}

template<typename Tsys>
void dump_initial_accp(const Tsys &sys){
	int rank = scfmpi_rank();
	int stat = 0;
	FILE *fp = 0;
	if(0 == rank){
		fp = fopen("init_accp.dat", "wx");
		if(!fp){
			fprintf(stderr, "File exists: init_accp.dat\n");
		       	stat = 1;
		}
	}
	shareint(&stat, 1);
	if(stat) return;
	
	int nloc = sys.getNumberOfParticleLocal();
	int ntot = nloc;
	scfmpi_sum_int(&ntot, 1);

	struct ForceID : Force{
		long id;
	};

	ForceID *fo = new ForceID[rank ? nloc : ntot];

	if(rank){
		for(int i=0; i<nloc; i++){
			fo[i].acc = sys[i].acc;
			fo[i].pot = sys[i].pot;
			fo[i].id  = sys[i].id;
		}
		MPI_Send(&nloc, 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
		MPI_Send(fo, sizeof(ForceID)*nloc, MPI_BYTE, 0, 102, MPI_COMM_WORLD);
	}else{ // root
		for(int i=0; i<nloc; i++){
			fo[i].acc = sys[i].acc;
			fo[i].pot = sys[i].pot;
			fo[i].id  = sys[i].id;
		}
		int icurr = nloc;
		int nprocs = scfmpi_size();
		for(int src=1; src<nprocs; src++){
			MPI_Recv(&nloc, 1, MPI_INT, src, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(fo+icurr, sizeof(ForceID)*nloc, MPI_BYTE, src, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			icurr += nloc;
		}

		std::sort(fo, fo+ntot, [](auto &a, auto &b){ return a.id < b.id; } );

		fprintf(fp, "%d\n", ntot);
		for(int i=0; i<ntot; i++){
			assert(i == fo[i].id);
			fprintf(fp, "%ld\t%A\t%A\t%A\t%A\n", 
					(long)fo[i].id, fo[i].acc.x, fo[i].acc.y, fo[i].acc.z, fo[i].pot);
		}
		// puts("TEST PASSED");
		fclose(fp);
	}

	delete [] fo;

	scfmpi_barr();
}

template<typename Tsys>
void set_particle_id(Tsys & sys, const int n_disk){
	PS::S64 base = n_disk * scfmpi_rank();
	for(int i=0; i<n_disk; i++){
		sys[i].id = base + i;
	}
}

template<typename Tsys>
void copy_ptcls_scf2fdps(Tsys & sys,
                         const int & n_halo,
                         const int & n_disk,
                         const double m[],
                         const double x[],
                         const double y[],
                         const double z[],
                         const double vx[],
                         const double vy[],
                         const double vz[],
			 const int id[]){
    sys.setNumberOfParticleLocal(n_disk);
    assert(n_halo + n_disk <= NMAX);
    for(auto i=0, nd=n_disk; i<nd; i++){
        sys[i].mass  = m[n_halo+i];
        sys[i].pos.x = x[n_halo+i];
        sys[i].pos.y = y[n_halo+i];
        sys[i].pos.z = z[n_halo+i];
        sys[i].vel.x = vx[n_halo+i];
        sys[i].vel.y = vy[n_halo+i];
        sys[i].vel.z = vz[n_halo+i];
        sys[i].id    = id[n_halo+i];
    }
}

template<typename Tsys>
void copy_ptcls_fdps2scf(const Tsys & sys,
                         const int & n_halo,
                         int & n_disk,
                         double m[],
                         double x[],
                         double y[],
                         double z[],
                         double vx[],
                         double vy[],
                         double vz[],
			 int    id[]){
    n_disk = sys.getNumberOfParticleLocal();

    if(n_halo + n_disk > NMAX){
	    int rank = scfmpi_rank();
	    printf("Overflow, rank=%d, NMAX=%d, n_halo+n_disk=%d, n_disk=%d\n",
			    rank, NMAX, n_halo+n_disk, n_disk);
    }

    assert(n_halo + n_disk <= NMAX);
    for(auto i=0, nd=n_disk; i<nd; i++){
        m[n_halo+i]  = sys[i].mass;
        x[n_halo+i]  = sys[i].pos.x;
        y[n_halo+i]  = sys[i].pos.y;
        z[n_halo+i]  = sys[i].pos.z;
        vx[n_halo+i] = sys[i].vel.x;
        vy[n_halo+i] = sys[i].vel.y;
        vz[n_halo+i] = sys[i].vel.z;
        id[n_halo+i] = sys[i].id;
    }
}


template<typename Tsys>
void copy_forces_fdps2scf(const Tsys & sys,
                          const int & n_halo,
                          int & n_disk,
                          double ax[],
                          double ay[],
                          double az[],
                          double pot[]){
    n_disk = sys.getNumberOfParticleLocal();
    assert(n_halo + n_disk <= NMAX);
    for(auto i=0, nd=n_disk; i<nd; i++){
        ax[n_halo+i]  = sys[i].acc.x;
        ay[n_halo+i]  = sys[i].acc.y;
        az[n_halo+i]  = sys[i].acc.z;
        pot[n_halo+i] = sys[i].pot;
    }
}

using TreeType = PS::TreeForForceLong<Force, EPI, EPJ>::Quadrupole;
using SpjType  = PS::SPJQuadrupole;

int main(int argc, char **argv)
{
	struct Params param;
	double t0mpi,t1mpi,t2mpi,t01mpi, t12mpi;
	int n,ndisk,nhalo,nbodies;
	int me,n_pes;
	double tnow,tpos,tvel,etotal;
	static double mass[NMAX],x[NMAX],y[NMAX],z[NMAX],vx[NMAX],
		vy[NMAX],vz[NMAX],ax[NMAX],ay[NMAX],az[NMAX],
		pot[NMAX];
	static int id[NMAX];
	double dblfact[LLMAX+1],twoalpha[LLMAX+1],c1[LLMAX+1][LLMAX+1],
		c2[LLMAX+1][LLMAX+1],c3[NNMAX+1],anltilde[NNMAX+1][LLMAX+1],
		lplusm[LLMAX+1][LLMAX+1],lmm[LLMAX+1][LLMAX+1],
		coeflm[LLMAX+1][LLMAX+1];

	//MPI_Init(&argc,&argv);
        PS::Initialize(argc, argv);

	Timer::initialize();
	Timer::beg(Timer::TOTAL);

	MPI_Comm_size(MPI_COMM_WORLD,&n_pes);

	MPI_Comm_rank(MPI_COMM_WORLD,&me);

	t0mpi = MPI_Wtime();

	startout(me);

	inparams(&param,me);

        //void (*CalcForceSp)() = CalcForceSpMono.*operator();
        //if(param.mm_order == 2){
        //    CalcForceSp = (void *)CalcForceSpQuad();
        //}
        
	int incont = param.incont;
	int nsteps = param.nsteps;
	EPI::eps = param.eps;

	if(0 == strcmp(param.softening_type,"PLM")){
		TreeParams::set_plummer(param.eps);
	}
	if(0 == strcmp(param.softening_type,"SPL")){
		TreeParams::set_spline(param.eps);
	}

	Timer::beg(Timer::IO_READ);
	int skip_halo = 0 == strcmp(param.basisset, "EX");
	inbods(incont,&ndisk,&nhalo,&tnow,mass,x,y,z,vx,vy,vz,id,&etotal,me,n_pes, skip_halo);
	MPI_Barrier(MPI_COMM_WORLD);
	Timer::end(Timer::IO_READ);

        /*
        auto CalcForceEp = &CalcGravity<EPI, EPJ, Force>::operator();
        auto CalcForceSp = &CalcGravity<EPI, SpjType, Force, 0>::operator();
        if(param.mm_order == 2){
            CalcForceSp = &CalcGravity<EPI, SpjType, Force, 2>::operator();
        }
        */
        PS::ParticleSystem<FP> system_disk;
        system_disk.initialize();
        system_disk.setNumberOfParticleLocal(ndisk);
        copy_ptcls_scf2fdps(system_disk, nhalo, ndisk, mass, x, y, z, vx, vy, vz, id);
	if(param.incont == 0){
		set_particle_id(system_disk, ndisk);
	}
        const auto coef_ema = 0.3;
        PS::DomainInfo dinfo;
        dinfo.initialize(coef_ema);
	{
		int size = scfmpi_size();
		int ny = (int)sqrt(0.1 + size);
		for(; ny>=1; ny--){
			if(0 == size%ny) break;
		}
		int nx = size / ny;
		fprintf(stderr, "Domain: nx=%d, ny=%d\n", nx, ny);
		dinfo.setNumberOfDomainMultiDimension(nx, ny);
	}
        dinfo.decomposeDomainAll(system_disk);
        system_disk.exchangeParticle(dinfo);
        copy_ptcls_fdps2scf(system_disk, nhalo, ndisk, mass, x, y, z, vx, vy, vz, id);

        auto CalcGravityEp = &CalcGravityP0<EPI, EPJ, Force>;
        auto CalcGravitySp = &CalcGravityP0<EPI, SpjType, Force>;
        if(strcmp(param.softening_type,"SPL") == 0){
            CalcGravityEp = &CalcGravitySpl<EPI, EPJ, Force>;
            if(param.mm_order == 2){
                CalcGravitySp = &CalcGravityP2Eps0<EPI, SpjType, Force>;
            } else {
                CalcGravitySp = &CalcGravityP0Eps0<EPI, SpjType, Force>;
            }
            
        } else {
            if(param.mm_order == 2){
                CalcGravitySp = &CalcGravityP2<EPI, SpjType, Force>;
            }
        }
        TreeType tree_disk;
        tree_disk.initialize(ndisk, param.theta, param.n_leaf_limit, param.n_group_limit);
        tree_disk.calcForceAllAndWriteBack(CalcGravityEp,
                                           CalcGravitySp,
                                           system_disk,
                                           dinfo);
	double tree_weight = TreeProf::costsum;
	TreeProf::costsum = 0;
        std::cerr<<"system_disk[0].acc= "<<system_disk[0].acc<<" pot= "<<system_disk[0].pot<<std::endl;
        copy_forces_fdps2scf(system_disk, nhalo, ndisk, ax, ay, az, pot);

	Timer::beg(Timer::IO_WRITE);
	dump_initial_accp(system_disk);
	Timer::end(Timer::IO_WRITE);
        
//	printf("passed inbods\n");
//	printf("me = %d; ndisk = %d  nhalo = %d\n",me,ndisk,nhalo);
//
	nbodies=ndisk+nhalo;

	initsys(&param,ndisk,nhalo,nbodies,&tnow,&tpos,&tvel,dblfact,
		twoalpha,c1,c2,c3,anltilde,lplusm,lmm,coeflm,
		mass,x,y,z,vx,vy,vz,ax,ay,az,pot,id,&etotal,t0mpi,me,n_pes);
        
	t1mpi = MPI_Wtime();
	t01mpi = t1mpi - t0mpi;
	// printf("initsys time = %lf on proc %d\n",t01mpi,me);

	for(n=1;n<=nsteps;n++){
            double dtime = param.dtime;
	    Timer::beg(Timer::STEP_POS);
            steppos(nbodies,x,y,z,vx,vy,vz,dtime,&tnow,&tpos);
	    Timer::end(Timer::STEP_POS);

            copy_ptcls_scf2fdps(system_disk, nhalo, ndisk, mass, x, y, z, vx, vy, vz, id);

	    Timer::beg(Timer::TREE_PREP);
            // dinfo.decomposeDomainAll(system_disk);
            dinfo.decomposeDomainAll(system_disk, tree_weight);
            system_disk.exchangeParticle(dinfo);
	    Timer::end(Timer::TREE_PREP);


            copy_ptcls_fdps2scf(system_disk, nhalo, ndisk, mass, x, y, z, vx, vy, vz, id);

	    Timer::beg(Timer::TREE_WALK);
            tree_disk.calcForceAllAndWriteBack(CalcGravityEp,
                                               CalcGravitySp,
                                               system_disk,
                                               dinfo);
	    tree_weight = TreeProf::costsum;
	    TreeProf::costsum = 0;
	    Timer::end(Timer::TREE_WALK);
#if 0
	    if(0 == strcmp(TreeParams::softening_type, "Spline")){
		    int max_len = TreeProf::max_list_length;
		    scfmpi_max_int(&max_len, 1);
		    if(0 == me){
			    fprintf(stderr, "spline: max_list_length: %d\n", max_len);
		    }
		    TreeProf::max_list_length = 0;
	    }
#endif

	    Timer::beg(Timer::TREE_BARR, true);
	    // print_average_etc(tree_weight);
	    MPI_Barrier(MPI_COMM_WORLD);
	    Timer::end(Timer::TREE_BARR);

            copy_forces_fdps2scf(system_disk, nhalo, ndisk, ax, ay, az, pot);
            nbodies=ndisk+nhalo;
            stepsys(&param,n,ndisk,nhalo,nbodies,&tnow,&tpos,&tvel,
                    dblfact,twoalpha,c1,c2,c3,anltilde,lplusm,lmm,
                    coeflm,mass,x,y,z,vx,vy,vz,ax,ay,az,pot,id,&etotal,
                    t0mpi,me,n_pes);
        }
        
	t2mpi = MPI_Wtime();
	t12mpi = t2mpi-t1mpi;
//	printf("total runtime = %d on proc %lf\n",me,t12mpi);

	endrun(t01mpi,t12mpi,me);

	Timer::end(Timer::TOTAL);
	if(0==me){
		Timer::show(stdout);
	}
	{
		char fname[256];
		sprintf(fname, "prof.%03d.%03d", scfmpi_rank(), scfmpi_size());
		FILE *fp = fopen(fname, "w");
		if(fp){
			Timer::show(fp);
			TreeProf::show_local(fp);
		}
		fclose(fp);
	}
	TreeProf::show_global();

	//MPI_Finalize();
        PS::Finalize();

	return 0;
}

