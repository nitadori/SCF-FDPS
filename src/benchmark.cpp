#include <cstdio>
#include <cstdlib>
#include "user_defined.hpp"
#include "timer.hpp"

using SPJ = ParticleSimulator::SPJQuadrupole;

int main(){
	enum{
		NI = 2048,
		NJ = 2048,
	};

	const double eps = 1./20.;
	const double eps2 = eps*eps;
	const double epsinv = 1.0/eps;

	static EPI epi[NI];
	static Force force[NI];
	static EPJ epj[NJ];
	static SPJ spj[NJ];

	srand48(20220526);
	for(int i=0; i<NI; i++){
		epi[i].pos.x = drand48() - 0.5;
		epi[i].pos.y = drand48() - 0.5;
		epi[i].pos.z = drand48() - 0.5;
	}

	for(int j=0; j<NJ; j++){
		epj[j].pos.x = spj[j].pos.x = drand48() - 0.5;
		epj[j].pos.y = spj[j].pos.y = drand48() - 0.5;
		epj[j].pos.z = spj[j].pos.z = drand48() - 0.5;
		epj[j].mass  = spj[j].mass  = (1./NJ) * (drand48() + 0.5);

		spj[j].quad.xx = drand48();
		spj[j].quad.yy = drand48();
		spj[j].quad.zz = drand48();
		spj[j].quad.xy = drand48();
		spj[j].quad.yz = drand48();
		spj[j].quad.xz = drand48();
	}

	// warm-up
	nbody_m256d(eps2, epsinv, epi, NI, epj, NJ, force);
	nbody_m256d_quad(eps2, epi, NI, spj, NJ, force);
	nbody_m256d_spline(eps2, epsinv, epi, NI, epj, NJ, force);

	for(int i=0; i<10; i++){
		auto tick0 = get_utime();
		nbody_m256d(eps2, epsinv, epi, NI, epj, NJ, force);
		auto tick1 = get_utime();
		nbody_m256d_quad(eps2, epi, NI, spj, NJ, force);
		auto tick2 = get_utime();
		nbody_m256d_spline(eps2, epsinv, epi, NI, epj, NJ, force);
		auto tick3 = get_utime();

		double dt_epj = tick2second(tick1 - tick0);
		double dt_spj = tick2second(tick2 - tick1);
		double dt_spl = tick2second(tick3 - tick2);

		double ns1 = dt_epj / NI / NJ * 1.e9;
		double ns2 = dt_spj / NI / NJ * 1.e9;
		double ns3 = dt_spl / NI / NJ * 1.e9;

		printf("epj: %.3f nsec, spj: %.3f nsec, ratio = %.3f\n", ns1, ns2, ns2/ns1);
		printf("spl: %.3f nsec, spj: %.3f nsec, ratio = %.3f\n", ns3, ns2, ns2/ns3);
		printf("spline: max_length = %d\n", TreeProf::max_list_length);
	}
}
