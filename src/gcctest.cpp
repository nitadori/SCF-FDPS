void mulv(const int n, double  * __restrict v, const double s){
	for(int i=0; i<n; i++){
		v[i] *= s;
	}
}

void scale (
		const int n, 
		double  * __restrict dst, 
		const double  * __restrict src, 
		const double s)
{
	for(int i=0; i<n; i++){
		dst[i] = s * src[i];
	}
}
