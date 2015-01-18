/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

#define ___xxxV_mulUOVO___( _eqsini_, _eqs_ ) {                              \
	T* pO = ((T*)pOut);                                                      \
	T* pA = ((T*)&A) ;                                                       \
	T* pB = ((T*)&B) ;                                                       \
	{                                                                        \
		                                    (*pO) _eqsini_ (*pA) * (*pB) ;   \
		for(int k=0+1;k< N;k++){pA++ ;pB++; (*pO) _eqs_    (*pA) * (*pB) ;}  \
	}                                                                        \
	for(int i=1;i<N;i++)                                                     \
	{                                                                        \
		pO++;                                                                \
		pA = ((T*)&A) + i  ;                                                 \
		pB = ((T*)&B) ;                                                      \
		                                    (*pO) _eqsini_ (*pA) * (*pB) ;   \
		for(int k=1  ;k<=i;k++){pA+=N;pB++; (*pO) _eqs_    (*pA) * (*pB) ;}  \
		for(int k=i+1;k< N;k++){pA++ ;pB++; (*pO) _eqs_    (*pA) * (*pB) ;}  \
	}                                                                        \
}                                                                            \

template<class T,int N> inline void
iniV_mulUOVO(KVector<T,N>* const pOut, const KUSymMat<T,N>& A,const KVector<T,N>& B){___xxxV_mulUOVO___( =, += )}
template<class T,int N> inline void
addV_mulUOVO(KVector<T,N>* const pOut, const KUSymMat<T,N>& A,const KVector<T,N>& B){___xxxV_mulUOVO___(+=, += )}
template<class T,int N> inline void
minV_mulUOVO(KVector<T,N>* const pOut, const KUSymMat<T,N>& A,const KVector<T,N>& B){___xxxV_mulUOVO___(-=, -= )}

#undef ___xxxV_mulUOVO___


/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
template<class T,int N> inline
void addU_squVOVT(KUSymMat<T,N>* const pOut, const KVector<T,N>& A)
{
	for(int i=0;i<N;i++)
	{
		double& Ai= A(i);
		for(int j=i;j<N;j++) (*pOut)(i,j) += (Ai) * A(j) ;
	}
}
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
template<class T,int N> inline
void iniU_squVOSVT(KUSymMat<T,N>* pOut, const KVector<T, N>& A, const T& a)
{
	for(int i=0;i<N;i++)
	{
		double t = (A(i))* a;
		for(int j=i;j<N;j++) (*pOut)(i,j) = t * A(j) ;
	}
}

template<class T,int N> inline
void addU_squVOSVT(KUSymMat<T,N>* pOut, const KVector<T, N>& A, const T& a)
{
	for(int i=0;i<N;i++)
	{
		double t = (A(i))* a;
		for(int j=i;j<N;j++) (*pOut)(i,j) += t * A(j) ;
	}
}

/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
#define ___xxxS_mulVTUOVO___( _eqsini_, _eqs_ ) {                            \
	const T* pA = (const T*)&A;                                              \
	const T* pB = (const T*)&B;                                              \
	const T* pC = (const T*)&C;                                              \
	{                                                                        \
		                                   (*pOut)_eqsini_(*pA)*(*pB)*(*pC); \
		for(int j=1  ;j< N;j++){pC++;pB++; (*pOut)_eqs_   (*pA)*(*pB)*(*pC);}\
	}                                                                        \
	for(int i=1;i<N;i++)                                                     \
	{                                                                        \
		pA++;pC = (const T*)&C;pB=((const T*)&B)+i;                          \
		                                   (*pOut)_eqs_   (*pA)*(*pB)*(*pC); \
		for(int j=1  ;j<=i;j++){pC++;pB+=N;(*pOut)_eqs_   (*pA)*(*pB)*(*pC);}\
		for(int j=i+1;j< N;j++){pC++;pB++ ;(*pOut)_eqs_   (*pA)*(*pB)*(*pC);}\
	}                                                                        \
}                                                                            \

template<class T,int N> inline void
iniS_mulVTUOVO(T* const pOut,const KVector<T,N>& A,const KUSymMat<T,N>& B,const KVector<T,N>& C)
{
	___xxxS_mulVTUOVO___( = , += );
}

template<class T,int N> inline void
addS_mulVTUOVO(T* const pOut,const KVector<T,N>& A,const KUSymMat<T,N>& B,const KVector<T,N>& C)
{
	___xxxS_mulVTUOVO___( += , += );
}

#undef ___xxxS_mulVTUOVO___
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

template<class T,int N> inline
void iniS_quadVTUOVO(T* const pOut,  const KUSymMat<T,N>& B, const KVector<T, N>& A)
{
	const T* pA = (const T*)&A;
	const T* pB = (const T*)&B;
	const T* pC = (const T*)&A;
	{
		                                 (*pOut) =    (*pA)*(*pB)*(*pC);
		for(int j=1  ;j<N;j++){pC++;pB++;(*pOut)+= 2.*(*pA)*(*pB)*(*pC);}
	}
	for(int i=1;i<N;i++)
	{
		pA++;pC =pA;pB+=(i+1);
		                                 (*pOut) +=    (*pA)*(*pB)*(*pC) ;
		for(int j=i+1;j<N;j++){pC++;pB++;(*pOut) += 2.*(*pA)*(*pB)*(*pC) ;}
	}

}




/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
