#include <stdio.h>
#include <stdint.h>
#include <string.h>

int main(int argc, char **argv){
	if(argc < 3){
		fprintf(stderr, "usage: ./b2ad <infile(binary)> <outfile(ascii)>\n");
		return 1;
	}
	FILE *fpin  = fopen(argv[1], "r");
	if(!fpin){
		fprintf(stderr, "open failed: %s\n", argv[1]);
		return 2;
	}
	FILE *fpout = fopen(argv[2], "wx");
	if(!fpout){
		fprintf(stderr, "open failed: %s\n", argv[2]);
		return 3;
	}

	int64_t ndtot;
	double tnow, etotal;

	fread(&ndtot,  8, 1, fpin);
	fread(&tnow,   8, 1, fpin);
	fread(&etotal, 8, 1, fpin);

	fprintf(fpout,"%9d %14.6E %14.6E\n", (int)ndtot, tnow, etotal);

	for(int k=0; k<ndtot; k++){
		double buf[9];
		fread(buf, 8, 9, fpin);
		int64_t id;
		memcpy(&id, buf, 8);
		fprintf(fpout, "%8d%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",
				(int)id, buf[1], buf[2], buf[3], buf[4], buf[5], buf[6], buf[7], buf[8]);
	}
}
