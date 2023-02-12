#include <mpi.h>
#include "scfdhmpi.h"

void scfmpi_sum_double(double *ptr, int nwords)
{
	MPI_Allreduce(MPI_IN_PLACE, ptr, nwords, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void scfmpi_sum_int(int *ptr, int nwords)
{
	MPI_Allreduce(MPI_IN_PLACE, ptr, nwords, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void scfmpi_max_int(int *ptr, int nwords)
{
	MPI_Allreduce(MPI_IN_PLACE, ptr, nwords, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
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
void mpiiget(int *rbuf,int *sbuf,int n,int idest,int isrc,int me)
{
	int itag;
	MPI_Status istatus;

	itag = isrc;

	if(me == isrc)
		MPI_Send(sbuf,n,MPI_INT,idest,itag,MPI_COMM_WORLD);
	if(me == idest)
		MPI_Recv(rbuf,n,MPI_INT,isrc,itag,MPI_COMM_WORLD,&istatus);
}


void share(double *x,int n)
{
	MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
}


void shareint(int *x,int n)
{
	MPI_Bcast(x,n,MPI_INT,0,MPI_COMM_WORLD);
}

void scfmpi_barr(){
	MPI_Barrier(MPI_COMM_WORLD);
}

int scfmpi_rank(void){
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
}

int scfmpi_size(void){
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	return size;
}

void scfmpi_abort(int e){
	MPI_Abort(MPI_COMM_WORLD, e);
}
