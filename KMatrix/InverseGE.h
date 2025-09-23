#ifndef Headder_KMatrix_GEInverse
#define Headder_KMatrix_GEInverse


template<int N> inline
KMatrix<double,N,N> invGE(const KMatrix<double,N,N>& aA)
{ // https://thira.plavox.info/blog/2008/06/_c.html
	//掃き出し法, Gaussian Elimination
	const double EPS = 1.0e-9;
	KMatrix<double,N,N> A = aA ;
	KMatrix<double,N,N> B      ; //ここに逆行列が入る
	iniM_IS(&B,1.);
	for (int i = 0; i<N; i++) {
		double buf = 1./(std::max)(EPS,A(i,i));
		for (int j = 0; j<N; j++) {
			A(i,j) *= buf;
			B(i,j) *= buf;
		}
		for (int j = 0; j<N; j++) {
			if (i!=j) {
				buf = A(j,i);
				for (int k = 0; k<N; k++) {
					A(j,k) -= A(i,k)*buf;
					B(j,k) -= B(i,k)*buf;
				}
			}
		}
	}
	return B;
}


#endif


