double unl[nmax+1][lmax+1];
double plm[lmax+1][lmax+1];
double cosmphi[lmax+1], sinmphi[lmax+1];

for(int l=0; l<=lmax; l++){
	for(int m=0; m<=l; m++){
		double ctemp = plm[l][m] * cosmphi[m];
		double stemp = plm[l][m] * sosmphi[m];
		for(int n=0; n<=nmax; n++){
			cossum[m][l][n] += unl[n][l] * ctmp;
			sinsum[m][l][n] += unl[n][l] * stmp;
		}
	}
}
