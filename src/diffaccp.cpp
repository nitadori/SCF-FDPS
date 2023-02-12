#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

struct Force{
	long id;
	double ax, ay, az;
	double pot;

	int read(FILE *fp){
		return fscanf(fp, "%ld %lf %lf %lf %lf", 
				&id, &ax, &ay, &az, &pot); 
	}
};

struct FDiff{
	double aerr_abs;
	double aerr_rel;
	double perr_abs;
	double perr_rel;

	FDiff(const Force &val, const Force &ref){
		auto norm = [](auto x, auto y, auto z) { return std::sqrt(x*x +  y*y + z*z); };
		auto dx = val.ax - ref.ax;
		auto dy = val.ay - ref.ay;
		auto dz = val.az - ref.az;
		auto dp = val.pot - ref.pot;

		aerr_abs = norm(dx, dy, dz);
		aerr_rel = aerr_abs / norm (ref.ax, ref.ay, ref.az);
		perr_abs = std::fabs(dp);
		perr_rel = perr_abs / std::fabs(ref.pot);
	}
};

std::vector<Force> readfile(const char *filename){
	FILE *fp = fopen(filename, "r");
	assert(fp);
	int n;
	fscanf(fp, "%d", &n);
	fprintf(stderr, "n = %d in %s\n", n, filename);
	
	std::vector<Force> fv(n);
	for(int i=0; i<n; i++){
		fv[i].read(fp);
	}
	fclose (fp);

	return fv;
}

int main(int argc, char **argv){
	assert(argc >= 3);

	auto v1 = readfile(argv[1]);
	auto v2 = readfile(argv[2]);

	assert(v1.size() == v2.size());

	const int n = v1.size();
	std::vector<double> aa(n), ar(n), pa(n), pr(n);
	for(int i=0; i<n; i++){
		FDiff d(v1[i], v2[i]);
		aa[i] = d.aerr_abs;
		ar[i] = d.aerr_rel;
		pa[i] = d.perr_abs;
		pr[i] = d.perr_rel;
	}
	auto sortvec = [](auto &v){ std::sort(v.begin(), v.end()); };
	sortvec(aa);
	sortvec(ar);
	sortvec(pa);
	sortvec(pr);

	for(double x=0.0; x<1.0; x+=0.01){
		int i = n*x;
		if(i >= n) break;
		printf("%2.2f  %.3e  %.3e  %.3e  %.3e\n", 
				x, aa[i], ar[i], pa[i], pr[i]);
	}

	return 0;
}
