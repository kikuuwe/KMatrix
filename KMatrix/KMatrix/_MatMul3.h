
//#########################################################################################
template<int K,int L,int M,int N> inline
void iniM_mulMTMOMO(KMatrix<double,K,L>* const pOut
					, const KMatrix<double,M,K>& A
					, const KMatrix<double,M,N>& B
					, const KMatrix<double,N,L>& C)
{
	KMatrix<double,K,N> tAB;
	iniM_mulMTMO(&tAB,A,B);
	iniM_mulMOMO(pOut,tAB,C);
}
//#########################################################################################
template<int K,int L,int M,int N> inline
void iniM_mulMOMOMT(KMatrix<double,K,L>* const pOut
					, const KMatrix<double,K,M>& A
					, const KMatrix<double,M,N>& B
					, const KMatrix<double,L,N>& C)
{
	KMatrix<double,K,N> tAB;
	iniM_mulMOMO(&tAB,A,B);
	iniM_mulMOMT(pOut,tAB,C);
}
template<int K,int L,int M,int N> inline
void addM_mulMOMOMT(KMatrix<double,K,L>* const pOut
					, const KMatrix<double,K,M>& A
					, const KMatrix<double,M,N>& B
					, const KMatrix<double,L,N>& C)
{
	KMatrix<double,K,N> tAB;
	iniM_mulMOMO(&tAB,A,B);
	addM_mulMOMT(pOut,tAB,C);
}
//#########################################################################################

template<class T,int M,int N>inline
void iniM_mulUOMOUO(KMatrix<T,M,N>* const pO,  const KUSymMat<T,M>& A, const KMatrix<T,M,N>& B, const KUSymMat<T,N>& C)
{
	KMatrix<T,M,N> tBC  ;
	init_mulOO(&tBC	,		B	 , C  );
	init_mulOO(pO	  ,  A , tBC  );
}
