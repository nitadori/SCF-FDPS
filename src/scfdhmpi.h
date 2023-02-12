//#define NMAX 60000
// #define NMAX 100000
#define NMAX 1'100'000
#define NNMAX 64
#define LLMAX 32
//#define NO_EXCLUSIVE_IO
#define OUTBODS_BINARY

struct Params{
        enum {TEXT_LEN = 8};
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
        double eps;
        double abasis;
        double abasis2;
        int zeroodd;
        int zeroeven;
        int fixacc;
        int ubodlog;
	int outpcoef;
	// for NFW external halo
	double cnfw;
	double rmax;
	double mh;

        // for tree
        char softening_type[TEXT_LEN];
        double theta;
        int n_leaf_limit;
        int n_group_limit;
        int mm_order;
};

struct Common{
	static inline double ext_pot_energy = 0.0;
};

void accp_CB(struct Params *param,int ni,int ne,double *mass,
		double *x,double *y,double *z,double *ax,double *ay,
		double *az,double *pot,double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1]);

void accp_LH(struct Params *param,int ni,int ne,double *mass,
		double *x,double *y,double *z,double *ax,double *ay,
		double *az,double *pot,double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1]);

void accpot(struct Params *param,int ndisk,int nhalo,int nbodies,
		double *tnow,double *mass,double *x,double *y,double *z,
		double *ax,double *ay,double *az,double *pot,
		double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1]);

void accp_tree(struct Params *param,int ndisk,
		double *mass,double *x,double *y,double *z,
		double *ax,double *ay,double *az,double *pot);

void corracc(int nbodies,double *mass,double *ax,double *ay,double *az);

void corrvel(double rcsign,int nbodies,double *tvel,double *tnow,double dtime,
		double *vx,double *vy,double *vz,double *ax,double *ay,
		double *az);

void endrun(double t01mpi, double t12mpi, int me);

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

void inbods(int incont,int *ndisk,int *nhalo,double *tnow,double *mass,
		double *x,double *y,double *z,double *vx,double *vy,
		double *vz,int *id,double *etotal,int me,int n_pes,
		int skip_halo);

void initpars(double *tpos,double *tvel,double *tnow);

void initsys(struct Params *param,int ndisk,int nhalo,int nbodies,
		double *tnow,double *tpos,double *tvel,double *dblfact,
		double *twoalpha,double c1[][LLMAX+1],
		double c2[][LLMAX+1],double *c3,double anltilde[][LLMAX+1],
		double lplusm[][LLMAX+1],double lmm[][LLMAX+1],
		double coeflm[][LLMAX+1],
		double *mass,double *x,double *y, double *z,
		double *vx,double *vy,double *vz,
		double *ax,double *ay,double *az,double *pot, int *id,
		double *etotal,double cputime,int me,int n_pes);

void initvel(int nbodies,double *vx,double *vy,double *vz,double *ax,
		double *ay,double *az,double dtime,double *tnow,double *tvel);

void inparams(struct Params *param,int me);

void mpiget(double *rbuf,double *sbuf,int n,int idest,int isrc,int me);
void mpiiget(int *rbuf,int *sbuf,int n,int idest,int isrc,int me);

void outbods(int ndisk,int nhalo,int nbodies,int initsnap,double *tnow,
		double *mass,double *x,double *y,double *z,double *vx,
		double *vy,double *vz,double *pot, int *id, double *etotal,int ubodlog,
		int me,int n_pes);

void outlog(struct Params *param,int nbodies,double *tnow,
		double *mass,double *x,double *y,double *z,double *vx,
		double *vy,double *vz,double *ax,double *ay,double *az,
		double *pot,double *etotal,double cputime,int me,int n_pes);

void outstate(struct Params *param,int kstep,int ndisk,int nhalo,
		int nbodies,double *tvel,double *tnow,double *mass,
		double *x,double *y,double *z,
		double *vx,double *vy,double *vz,double *ax,double *ay,
		double *az,double *pot, int *id, double cputime,double *etotal,
		int me,int n_pes);

void outterm(char *message,int n,int me);

void outputcoef(int nmax,int lmax,double *tnow,
                double sinsumhh[][LLMAX+1][LLMAX+1],
                double cossumhh[][LLMAX+1][LLMAX+1]);

void scfcoef_CB(struct Params *param,int ni,int ne,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double c1[][LLMAX+1],double c2[][LLMAX+1],
		double *c3,double anltilde[][LLMAX+1],double *dblfact,
		double *twoalpha,double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1]);

void scfcoef_LH(struct Params *param,int ni,int ne,
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

void stepsys(struct Params *param,int kstep,int ndisk,int nhalo,int nbodies,
		double *tnow,double *tpos,double *tvel,
		double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double *vx,double *vy,double *vz,
		double *ax,double *ay,double *az,double *pot, int *id,
		double *etotal,double cputime,int me,int n_pes);

void stepvel(int nbodies,double *vx,double *vy,double *vz,double *ax,
		double *ay,double *az,double dtime,double *tnow,
		double *tvel);

typedef void fun_scfcoef(struct Params *param,int ni,int ne,
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1],
		double *mass,double *x,double *y,double *z,
		double c1[][LLMAX+1],double c2[][LLMAX+1],
		double *c3,double anltilde[][LLMAX+1],double *dblfact,
		double *twoalpha,double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1]);

fun_scfcoef scfcoef_CB_mod;
fun_scfcoef scfcoef_LH_mod;
fun_scfcoef scfcoef_CB_avx;
fun_scfcoef scfcoef_LH_avx;

typedef void fun_accp(struct Params *param,int ni,int ne,double *mass,
		double *x,double *y,double *z,double *ax,double *ay,
		double *az,double *pot,double *dblfact,double *twoalpha,
		double c1[][LLMAX+1],double c2[][LLMAX+1],double *c3,
		double anltilde[][LLMAX+1],double lplusm[][LLMAX+1],
		double lmm[][LLMAX+1],double coeflm[][LLMAX+1],
		double sinsum[][LLMAX+1][LLMAX+1],
		double cossum[][LLMAX+1][LLMAX+1]);
fun_accp accp_CB_mod;
fun_accp accp_LH_mod;
fun_accp accp_CB_avx;
fun_accp accp_LH_avx;

void scfmpi_sum_int(int *ptr, int nwords);
void scfmpi_max_int(int *ptr, int nwords);
void scfmpi_sum_double(double *ptr, int nwords);
int  scfmpi_rank(void);
int  scfmpi_size(void);
void scfmpi_abort(int);
void scfmpi_barr();
void print_average_etc(double);

extern "C" {
	int ext_isnan(double x);
	int ext_isfinite(double x);
}
