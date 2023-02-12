#include <cstdio>
#include <cmath>
#include "scfdhmpi.h"

double factrl(int n) {
	int ntop,j;
	double a[34],arggam,ftrl;

	ntop=0;
	a[1]=1.0;

	if(n < 0){
		printf("negative factorial\n");
		exit(1);
	}else if(n <= ntop)
		ftrl=a[n+1];
	else if(n <= 32){
		for(j=ntop+1;j<=n;j++)
			a[j+1]=j*a[j];
		ntop=n;
		ftrl=a[n+1];
	}else{
		arggam=n+1.0;
		ftrl=exp(gammln(arggam));
	}

	return ftrl;
}	

double gammln(double xx)
{
	int j;
	double cof[7],stp,half,one,fpf,x,tmp,ser,gam;

	cof[1]=76.18009173;
	cof[2]=-86.50532033;
	cof[3]=24.01409822;
	cof[4]=-1.231739516;
	cof[5]=0.120858003e-2;
	cof[6]=-0.536382e-5,
	stp=2.50662827465;
	half=0.5;
	one=1.0;
	fpf=5.5;

	x=xx-one;
	tmp=x+fpf;
	tmp=(x+half)*log(tmp)-tmp;
	ser=one;

	for(j=1;j<=6;j++){
		x=x+one;
		ser=ser+cof[j]/x;
	}

	gam=tmp+log(stp*ser);

	return gam;
}

void getcoef_CB(int nmax,int lmax,double pi,double *dblfact,
	double anltilde[][LLMAX+1],double coeflm[][LLMAX+1],double *twoalpha,
	double lmm[][LLMAX+1],double lplusm[][LLMAX+1],double c1[][LLMAX+1],
	double c2[][LLMAX+1],double *c3)
{
        int l,m,n;
        double deltam0,knl,fctrll,fctrln;

        dblfact[1]=1.0;

        if(lmax >= 2)
                for(l=2;l<=lmax;l++)
                        dblfact[l]=dblfact[l-1]*(2*l-1);

        for(n=0;n<=nmax;n++){
                fctrln=factrl(n);
                for(l=0;l<=lmax;l++){
                        knl=4*n*(n+2*l+2)+(2*l+1)*(2*l+3);
                        fctrll=factrl(l);
                        anltilde[n][l]=-pow(2,4*l+4)*((n+l+1)/(pi*knl))
                                *(fctrln*fctrll/factrl(n+2*l+1))*fctrll;
                }
        }

        for(l=0;l<=lmax;l++){
                twoalpha[l]=2*(l+1);
                for(m=0;m<=l;m++){
                        deltam0=2.0;
                        if(m == 0) deltam0=1.0;
                        coeflm[l][m]=(2*l+1)*deltam0*factrl(l-m)/factrl(l+m);
                        lmm[l][m]=l-m;
                        if(l != m) lmm[l][m]=1.0/lmm[l][m];
                        lplusm[l][m]=l+m-1.0;
                }
        }

        for(n=1;n<=nmax;n++){
                c3[n]=1.0/(n+1.0);
                for(l=0;l<=lmax;l++){
                        c1[n][l]=2*n+twoalpha[l];
                        c2[n][l]=n-1+twoalpha[l];
                }
        }
}

void getcoef_LH(int nmax,int lmax,double pi,double *dblfact,
	double anltilde[][LLMAX+1],double coeflm[][LLMAX+1],double *twoalpha,
	double lmm[][LLMAX+1],double lplusm[][LLMAX+1],double c1[][LLMAX+1],
	double c2[][LLMAX+1],double *c3)
{
        int l,m,n;
        double deltam0,knl,arggam,fctrln,fctrl2l;

        dblfact[1]=1.0;

        if(lmax >= 2)
                for(l=2;l<=lmax;l++)
                        dblfact[l]=dblfact[l-1]*(2*l-1);

        for(n=0;n<=nmax;n++){
                fctrln=factrl(n);
                for(l=0;l<=lmax;l++){
                        knl=0.5*n*(n+4*l+3)+(l+1)*(2*l+1);
			arggam=2*l+1.5;
                        fctrl2l=exp(gammln(arggam));
                        anltilde[n][l]=-pow(2,8*l+4)*((n+2*l+1.5)/(pi*knl))
                                *(fctrln*fctrl2l/factrl(n+4*l+2))*fctrl2l;
                }
        }

        for(l=0;l<=lmax;l++){
                twoalpha[l]=4*l+3;
                for(m=0;m<=l;m++){
                        deltam0=2.0;
                        if(m == 0) deltam0=1.0;
                        coeflm[l][m]=(2*l+1)*deltam0*factrl(l-m)
				/factrl(l+m);
                        lmm[l][m]=l-m;
                        if(l != m) lmm[l][m]=1.0/lmm[l][m];
                        lplusm[l][m]=l+m-1.0;
                }
        }

        for(n=1;n<=nmax;n++){
                c3[n]=1.0/(n+1);
                for(l=0;l<=lmax;l++){
                        c1[n][l]=2*n+twoalpha[l];
                        c2[n][l]=n-1+twoalpha[l];
                }
        }
}
