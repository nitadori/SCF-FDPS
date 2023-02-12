#pragma once
#include <cstdio>
#include <cstdint>

#ifdef TIMER_SILENT
#  define TIMER_NUM 0
#else
#  define TIMER_NUM 1
#endif

#ifdef __linux__
#  include <time.h>
static int64_t get_utime(){
	timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);

	return ts.tv_nsec + ts.tv_sec*1000000000ll;
}
static double tick2second(uint64_t tick){
	return 1.e-9 * (double)tick;
}
#else
#  error
#endif

struct Timer_template_base {
	enum Items{
		SCF_COEF = 0,
		SCF_ACCP,
		STEP_POS,
		STEP_VEL,
		TREE_PREP,
		TREE_WALK,
		TREE_BARR,
		IO_READ,
		IO_WRITE,
		// don't touch below
		TOTAL,
		MISC,
		NUM_ITEMS,
	};
};
template <int>
struct Timer_template: Timer_template_base {
	static void flush(){}
	static void initialize(){}
	static void beg(const Items /*elem*/, bool const /*reuse*/ = false){}
	static void end(const Items /*elem*/, bool const /*reuse*/ = false, int64_t /*acc*/ = 0ll){}
	static void show(FILE *, char const * /*fmt*/ = ""){}
};

template<>
struct Timer_template<1> : Timer_template_base{
	static char const *name(int const i){
		static const char *strs[NUM_ITEMS] = {
			"SCF_COEF",
			"SCF_ACCP",
			"STEP_POS",
			"STEP_VEL",
			"TREE_PREP",
			"TREE_WALK",
			"TREE_BARR",
			"IO_READ",
			"IO_WRITE",
			// don't touch below
			"TOTAL",
			"MISC",
		};
		return strs[i];
	}

	static void flush(){
		for(int i=0; i<NUM_ITEMS; i++){
			time (i) = 0ll;
		}
	}

	static void initialize(){
		flush();
		tprev(1) = get_utime();
	}

	static void beg(const Items elem, bool const reuse = false){
		// puts(name(elem)); // nice for debugging
		if(reuse) time(elem) -= tprev();
		else time(elem) -= (tprev() = get_utime());
	}

	static void end(const Items elem, bool const reuse = false){
		if(reuse) time(elem) += tprev();
		else time(elem) += (tprev() = get_utime());
	}
	static void show(
			FILE *fp = stderr,
			const char *fmt = " %-12s : %e sec : %6.2f %%\n")

	{
		fflush(fp);

		time(MISC) = time(TOTAL);

		for(int i=0; i<NUM_ITEMS-2; i++){
			time(MISC) -= time(i);
		}

		for(int i=0; i<NUM_ITEMS; i++){
			double ratio = 100.0 * time(i)/ time(TOTAL);
			fprintf(fp, fmt, name(i), rtime(i), ratio);
		}

		fflush(fp);
	}

	private:
	static int64_t &time(int const i){
		static int64_t buf[NUM_ITEMS];
		return buf[i];
	}
	static int64_t &tprev(int const ch=0){
		static int64_t t[2]; /* 0 : previous time */
		                     /* 1 : initial time  */
		return t[ch];
	}
	static double rtime(int const i){
		return tick2second(time(i));
	}
};

using Timer = Timer_template<TIMER_NUM>;

