#include <cstdio>
#include <ps_types.hpp>
#include <tree.hpp>
namespace PS = ParticleSimulator;

struct Force{
	PS::F64vec acc;
	PS::F64    pot;    

	void clear() {
		acc = 0.0;
		pot = 0.0;
	}

};

struct FP{
	PS::S64    id;
	PS::F64    mass;
	PS::F64vec pos;
	PS::F64vec vel;
	PS::F64vec acc;
	PS::F64    pot;    

	void copyFromForce(const Force &force) {
		acc = force.acc;
		pot = force.pot;
	}
	PS::F64vec getPos() const {
		return pos;
	}

	void writeAscii(FILE* fp) const {
		fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
				this->id, this->mass,
				this->pos.x, this->pos.y, this->pos.z,
				this->vel.x, this->vel.y, this->vel.z);
	}

	void readAscii(FILE* fp) {
		fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
				&this->id, &this->mass,
				&this->pos.x, &this->pos.y, &this->pos.z,
				&this->vel.x, &this->vel.y, &this->vel.z);
	}

};

struct EPI{
	PS::F64vec pos;
	PS::F64    mass;

	static inline PS::F64 eps = 0.0;

	void copyFromFP(const FP &fp){ 
		mass = fp.mass;
		pos  = fp.pos;
	}

	PS::F64vec getPos() const {
		return pos;
	}

	PS::F64 getCharge() const {
		return mass;
	}


};

using EPJ = EPI;
using SPJ = ParticleSimulator::SPJQuadrupole;

struct TreeProf{
	static inline double nisum    = 0.0;
	static inline double niepjsum = 0.0;
	static inline double nispjsum = 0.0;

	static void show_local(FILE *);
	static void show_global(FILE *);

	static inline double cost_ratio = 2.4; //ratio for spj/epj
	static inline double costsum = 0.0;

	static inline int max_list_length = 0;

	static void count_epj(int ni, int nj){
		nisum += ni;
		niepjsum += ni*nj;

		costsum += ni*nj;
	}
	static void count_spj(int ni, int nj){
		nisum += ni;
		nispjsum += ni*nj;

		costsum += cost_ratio * ni*nj;
	}
};

template<typename Tepj>
void nbody_m256d(const double eps2, const double epsinv, const EPI epi[], const int ni, const Tepj epj[], const int nj, Force force[]);
// void nbody_m256d_nj1(const double eps2, const EPI epi[], const int ni, const SPJ epj[], const int nj, Force force[]);
void nbody_m256d_quad(const double eps2, const EPI epi[], const int ni, const SPJ epj[], const int nj, Force force[]);
void nbody_m256d_spline(const double eps2, const double epsinv, const EPI epi[], const int ni, const EPJ epj[], const int nj, Force force[]);

struct TreeParams{
	static inline decltype(&nbody_m256d<EPJ>) kernel_short;
	static inline decltype(&nbody_m256d_quad) kernel_long;
	
	static inline double eps2_short;
	static inline double epsinv;
	static inline double eps2_long;
	static inline const char *softening_type = 0;

	static void set_plummer(const double eps){
		eps2_short = eps2_long = eps*eps;
		epsinv = 1.0 / eps;

		kernel_short = nbody_m256d;
		kernel_long  = nbody_m256d_quad;

		TreeProf::cost_ratio = 2.4;

		softening_type = "Plummer";
	}

	static void set_spline(const double eps){
		eps2_short = eps*eps;
		eps2_long  = 0.0;
		epsinv     = 1.0 / eps;

		kernel_short = nbody_m256d_spline;
		kernel_long  = nbody_m256d_quad;

		TreeProf::cost_ratio = 2.2;

		softening_type = "Spline";
	}
};

#if 0
template <typename Tpi, typename Tpj, typename Tforce, int p=0>
struct CalcGravity{
	void operator() (const Tpi * ep_i,
                         const PS::S32 n_ip,
                         const Tpj * ep_j,
                         const PS::S32 n_jp,
                         Tforce * force) {
		const auto eps2 = Tpi::eps * Tpi::eps;
		(void)eps2;
		if constexpr ( p == 0 || p == 1){
#if 0
			for(PS::S32 i = 0; i < n_ip; i++){
				PS::F64vec xi = ep_i[i].getPos();
				PS::F64vec ai = 0.0;
				PS::F64 poti = 0.0;
				for(PS::S32 j = 0; j < n_jp; j++){
					PS::F64vec rij    = xi - ep_j[j].getPos();
					if(rij.x == 0.0 && rij.y == 0.0  && rij.z == 0.0) continue;
					PS::F64    r3_inv = rij * rij + eps2;
					PS::F64    r_inv  = 1.0/sqrt(r3_inv);
					r3_inv  = r_inv * r_inv;
					r_inv  *= ep_j[j].getCharge();
					r3_inv *= r_inv;
					ai     -= r3_inv * rij;
					poti   -= r_inv;
				}
				force[i].acc += ai;
				force[i].pot += poti;
			}
#else
			// const auto epsinv = 1.0 / Tpi::eps;
			// nbody_m256d(eps2, epsinv, ep_i, n_ip, ep_j, n_jp, force);
			TreeParams::kernel_short(TreeParams::eps2_short, TreeParams::epsinv, 
					ep_i, n_ip, ep_j, n_jp, force);
#endif
			TreeProf::count_epj(n_ip, n_jp);
		}
		else if constexpr ( p == 2){
#if 0
			for(PS::S32 ip=0; ip<n_ip; ip++){
				PS::F64vec xi = ep_i[ip].pos;
				PS::F64vec ai = 0.0;
				PS::F64 poti = 0.0;
				for(PS::S32 jp=0; jp<n_jp; jp++){
					PS::F64 mj = ep_j[jp].mass;
					PS::F64vec xj= ep_j[jp].pos;
					PS::F64vec rij= xi - xj;
					PS::F64 r2 = rij * rij + eps2;
					PS::F64mat qj = ep_j[jp].quad;
					PS::F64 tr = qj.getTrace();
					PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
							(qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
							(qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
					PS::F64 qrr = qr * rij;
					PS::F64 r_inv = 1.0f/sqrt(r2);
					PS::F64 r2_inv = r_inv * r_inv;
					PS::F64 r3_inv = r2_inv * r_inv;
					PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
					PS::F64 qrr_r5 = r5_inv * qrr;
					PS::F64 qrr_r7 = r2_inv * qrr_r5;
					PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
					PS::F64 B = -2.0*r5_inv;
#if 1
					ai -= A*rij + B*qr;
					poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
#else
					ai -= (mj*r3_inv)*rij;
					poti -= mj*r_inv;
#endif
				}
				force[ip].acc += ai;
				force[ip].pot += poti;
			}
#else
			// nbody_m256d_nj1(eps2, ep_i, n_ip, ep_j, n_jp, force);
			// nbody_m256d_quad(eps2, ep_i, n_ip, ep_j, n_jp, force);
			TreeParams::kernel_long(TreeParams::eps2_short, ep_i, n_ip, ep_j, n_jp, force);
#endif
			TreeProf::count_spj(n_ip, n_jp);
		}
		else{
			static_assert(p==0 || p==2);
			/*
			if(PS::Comm::getRank()==0){
				std::cerr<<"p must be 0, 1, or 2"<<std::endl;
			}
			PS::Comm::barrier();
			PS::Abort();
			*/
		}
	}
};
#endif

template <typename Tpi, typename Tpj, typename Tforce>
void CalcGravityImplP0(const PS::F64 eps2,
                       const Tpi * ep_i,
                       const PS::S32 n_ip,
                       const Tpj * ep_j,
                       const PS::S32 n_jp,
                       Tforce * force){
    // assert(0); // check for specialization
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
            PS::F64vec rij    = xi - ep_j[j].getPos();
            if(rij.x == 0.0 && rij.y == 0.0  && rij.z == 0.0) continue;
            PS::F64    r3_inv = rij * rij + eps2;
            PS::F64    r_inv  = 1.0/sqrt(r3_inv);
            r3_inv  = r_inv * r_inv;
            r_inv  *= ep_j[j].getCharge();
            r3_inv *= r_inv;
            ai     -= r3_inv * rij;
            poti   -= r_inv;
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}

template <typename Tpi, typename Tpj, typename Tforce>
void CalcGravityImplSpl(const PS::F64 rcut,
                        const Tpi * ep_i,
                        const PS::S32 n_ip,
                        const Tpj * ep_j,
                        const PS::S32 n_jp,
                        Tforce * force){
    // assert(0); // check for specialization
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
            PS::F64vec rij    = xi - ep_j[j].getPos();
            if(rij.x == 0.0 && rij.y == 0.0  && rij.z == 0.0) continue;
            PS::F64 r_sq = rij*rij;
            PS::F64 r_inv = r_sq > 0.0 ? 1.0 / sqrt(r_sq) : 0.0;
            PS::F64 r_cut_half = rcut * 0.5;
            //PS::F64 r_cut_half = Tpi::eps;
            PS::F64 r_cut_half_inv = r_sq > 0.0 ? 1.0 / r_cut_half : 0.0;
            PS::F64 r = r_sq * r_inv;
            PS::F64 u = r * r_cut_half_inv;
            PS::F64 r3_inv = r_inv * r_inv * r_inv;
            PS::F64 u2 = u * u;
            PS::F64 u3 = u * u2;
            PS::F64 u4 = u2 * u2;
            PS::F64 u5 = u2 * u3;
            PS::F64 u6 = u3 * u3;
            PS::F64 fr = r_inv;
            PS::F64 gr = r3_inv;
            if(2.0 > u && 1.0 <= u){
                gr = r3_inv * (-1.0/15.0 + (8.0/3.0)*u3 - 3.0*u4 + (6.0/5.0)*u5 - (1.0/6.0)*u6);
                fr = -1.0/15.0*r_inv
                    - r_cut_half_inv*((4.0/3.0)*u2 - u3 + (3.0/10.0)*u4 - (1.0/30.0)*u5)
                    + 8.0/5.0*r_cut_half_inv;
            }
            else if(1.0 > u){
                PS::F64 r_cut_half_inv = 1.0 / r_cut_half;
                PS::F64 r_cut_half3_inv = r_cut_half_inv * r_cut_half_inv * r_cut_half_inv;
                gr = r_cut_half3_inv*(4.0/3.0 - (6.0/5.0)*u2 + (1.0/2.0)*u3);
                fr = -2.0*r_cut_half_inv*((1.0/3.0)*u2 - (3.0/20.0)*u4 + (1.0/20.0)*u5) + 7.0/5.0*r_cut_half_inv;
            }
            ai   -= ep_j[j].mass * gr * rij;
            poti -= ep_j[j].mass * fr;
        }
        force[i].acc += ai;
        force[i].pot += poti;
    }
}

template <typename Tpi, typename Tpj, typename Tforce>
void CalcGravityImplP2(const PS::F64 eps2,
                       const Tpi * ep_i,
                       const PS::S32 n_ip,
                       const Tpj * ep_j,
                       const PS::S32 n_jp,
                       Tforce * force){
    // assert(0); // check for specialization
    for(PS::S32 ip=0; ip<n_ip; ip++){
        PS::F64vec xi = ep_i[ip].pos;
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 jp=0; jp<n_jp; jp++){
            PS::F64 mj = ep_j[jp].mass;
            PS::F64vec xj= ep_j[jp].pos;
            PS::F64vec rij= xi - xj;
            PS::F64 r2 = rij * rij + eps2;
            PS::F64mat qj = ep_j[jp].quad;
            PS::F64 tr = qj.getTrace();
            PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                           (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                           (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
            PS::F64 qrr = qr * rij;
            PS::F64 r_inv = 1.0f/sqrt(r2);
            PS::F64 r2_inv = r_inv * r_inv;
            PS::F64 r3_inv = r2_inv * r_inv;
            PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
            PS::F64 qrr_r5 = r5_inv * qrr;
            PS::F64 qrr_r7 = r2_inv * qrr_r5;
            PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
            PS::F64 B = -2.0*r5_inv;
#if 1
            ai -= A*rij + B*qr;
            poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
#else
            ai -= (mj*r3_inv)*rij;
            poti -= mj*r_inv;
#endif
        }
        force[ip].acc += ai;
        force[ip].pot += poti;
    }
}



template <typename Tpi, typename Tpj, typename Tforce>
void CalcGravityP0(const Tpi * ep_i,
                   const PS::S32 n_ip,
                   const Tpj * ep_j,
                   const PS::S32 n_jp,
                   Tforce * force){
    const auto eps2 = Tpi::eps*Tpi::eps;
    CalcGravityImplP0(eps2, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_epj(n_ip, n_jp);
}

template <typename Tpi, typename Tpj, typename Tforce>
void CalcGravityP2(const Tpi * ep_i,
                    const PS::S32 n_ip,
                    const Tpj * ep_j,
                    const PS::S32 n_jp,
                    Tforce * force){
    const auto eps2 = Tpi::eps*Tpi::eps;
    CalcGravityImplP2(eps2, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_spj(n_ip, n_jp);
}

template <typename Tpi, typename Tpj, typename Tforce>
void CalcGravitySpl(const Tpi * ep_i,
                    const PS::S32 n_ip,
                    const Tpj * ep_j,
                    const PS::S32 n_jp,
                    Tforce * force){
    const auto rcut = 2.0 * Tpi::eps;
    CalcGravityImplSpl(rcut, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_epj(n_ip, n_jp);
}

template <typename Tpi, typename Tpj, typename Tforce>
void CalcGravityP0Eps0(const Tpi * ep_i,
                       const PS::S32 n_ip,
                       const Tpj * ep_j,
                       const PS::S32 n_jp,
                       Tforce * force){
    const auto eps2 = 0.0;
    CalcGravityImplP0(eps2, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_epj(n_ip, n_jp);
}

template <typename Tpi, typename Tpj, typename Tforce>
void CalcGravityP2Eps0(const Tpi * ep_i,
                       const PS::S32 n_ip,
                       const Tpj * ep_j,
                       const PS::S32 n_jp,
                       Tforce * force){
    const auto eps2 = 0.0;
    CalcGravityImplP2(eps2, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::count_spj(n_ip, n_jp);
}

// Specializations -> <specializations.hpp>
/*
template <>
void CalcGravityP0(const EPI * ep_i,
                   const PS::S32 n_ip,
                   const EPJ * ep_j,
                   const PS::S32 n_jp,
                   Force * force){
    const auto eps2 = EPI::eps*EPI::eps;
    const auto epsinv = 1.0 / EPI::eps;
    nbody_m256d(eps2, epsinv, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::nisum += n_ip;
    TreeProf::niepjsum += n_ip * n_jp;
}

template <>
void CalcGravityP2(const EPI * ep_i,
                    const PS::S32 n_ip,
                    const SPJ * ep_j,
                    const PS::S32 n_jp,
                    Force * force){
    const auto eps2 = EPI::eps*EPI::eps;
    nbody_m256d_quad(eps2, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::nisum += n_ip;
    TreeProf::niepjsum += n_ip * n_jp;
}

template <>
void CalcGravitySpl(const EPI * ep_i,
                    const PS::S32 n_ip,
                    const EPJ * ep_j,
                    const PS::S32 n_jp,
                    Force * force){
    const auto eps2 = EPI::eps*EPI::eps;
    const auto epsinv = 1.0 / EPI::eps;
    nbody_m256d_spline(eps2, epsinv, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::nisum += n_ip;
    TreeProf::niepjsum += n_ip * n_jp;
}

template <>
void CalcGravityP0Eps0(const EPI * ep_i,
                       const PS::S32 n_ip,
                       const EPJ * ep_j,
                       const PS::S32 n_jp,
                       Force * force){
    const auto eps2 = 0.0;
    const auto epsinv = 0.0;

    nbody_m256d(eps2, epsinv, ep_i, n_ip, ep_j, n_jp, force);
    TreeProf::nisum += n_ip;
    TreeProf::niepjsum += n_ip * n_jp;
}

template <>
void CalcGravityP2Eps0(const EPI * ep_i,
                       const PS::S32 n_ip,
                       const SPJ * ep_j,
                       const PS::S32 n_jp,
                       Force * force){
    const auto eps2 = 0.0;
    nbody_m256d_quad(eps2, ep_i, n_ip, ep_j, n_jp, force);

    TreeProf::nisum += n_ip;
    TreeProf::niepjsum += n_ip * n_jp;
}
*/



    

