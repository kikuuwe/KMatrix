 
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

#define ___xxxM_mulVOVT___( _eqs_ ) {                             \
	T* pO = (T*)pOut;                                             \
	T* pA = (T*)&A;                                               \
	T* pB = (T*)&B;                                               \
	for(int i=0;i<M-1;i++)                                        \
	{                                                             \
		for(int j=0;j<N-1;j++){(*pO) _eqs_ (*pA)*(*pB);pO++;pB++;}\
		                      {(*pO) _eqs_ (*pA)*(*pB);}          \
    	pA++;pO++;pB=(T*)&B;                                      \
	}                                                             \
	{                                                             \
		for(int j=0;j<N-1;j++){(*pO) _eqs_ (*pA)*(*pB);pO++;pB++;}\
		                      {(*pO) _eqs_ (*pA)*(*pB);}          \
	}                                                             \
}                                                                 \

template<class T, int M,int N> inline void
iniM_mulVOVT(KMatrix<T,M,N>* const pOut, const KVector<T,M>& A,const KVector<T,N>& B){___xxxM_mulVOVT___(= )}
template<class T, int M,int N> inline void
addM_mulVOVT(KMatrix<T,M,N>* const pOut, const KVector<T,M>& A,const KVector<T,N>& B){___xxxM_mulVOVT___(+=)}
template<class T, int M,int N> inline void
minM_mulVOVT(KMatrix<T,M,N>* const pOut, const KVector<T,M>& A,const KVector<T,N>& B){___xxxM_mulVOVT___(-=)}


template<class T,int N> inline
void addU_mulVOVT(KUSymMat<T,N>* pOut, const KVector<T, N>& A, const KVector<T, N>& B)
{
	for(int i=0;i<N;i++)
	{
		double& Ai= A(i);
		for(int j=i;j<N;j++){(*pOut)(i,j) += Ai * B(j) ;}
	}
}

#undef ___xxxM_mulVOVT___
 
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

#define ___xxxM_mulVOSVT___( _eqs_ ) {                            \
	T* pO = (T*)pOut;                                             \
	T* pA = (T*)&A;                                               \
	T* pB = (T*)&B;                                               \
	double pAk;                                                   \
	{                                                             \
		pAk = (*pA)*k;                                            \
		                    {          (*pO) _eqs_ (pAk)*(*pB);}  \
		for(int j=1;j<N;j++){pO++;pB++;(*pO) _eqs_ (pAk)*(*pB);}  \
	}                                                             \
	for(int i=1;i<M;i++)                                          \
	{                                                             \
    	pA++;pO++;pB=(T*)&B;                                      \
		pAk = (*pA)*k;                                            \
		                    {          (*pO) _eqs_ (pAk)*(*pB);}  \
		for(int j=1;j<N;j++){pO++;pB++;(*pO) _eqs_ (pAk)*(*pB);}  \
	}                                                             \
}                                                                 \


template<class T,int M,int N> inline void
iniM_mulVOSVT(KMatrix<T,M,N>* const pOut,const KVector<T,M>& A,const T& k,const KVector<T,N>& B){___xxxM_mulVOSVT___(=)}


#undef ___xxxM_mulVOSVT___
 
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

#define ___xxxV_mulMOVO___( _eqsini_, _eqs_ ) {                    \
	T* pO = ((T*)pOut);                                            \
	T* pA = ((T*)&A );                                             \
	T* pB = ((T*)&B );                                             \
	{                                                              \
		                    {            *pO _eqsini_ (*pA)*(*pB);}\
		for(int k=1;k<N;k++){pB++;pA++ ; *pO _eqs_    (*pA)*(*pB);}\
	}                                                              \
	for(int i=1;i<M;i++)                                           \
	{                                                              \
		pO++;                                                      \
		pA++;                                                      \
		pB  = ((T*)&B );                                           \
		                    {            *pO _eqsini_ (*pA)*(*pB);}\
		for(int k=1;k<N;k++){pB++;pA++ ; *pO _eqs_    (*pA)*(*pB);}\
	}                                                              \
}                                                                  \



template<class T,int M,int N> inline void
iniV_mulMOVO(KVector<T,M>* const pOut,const KMatrix<T,M,N>& A,const KVector<T,N>& B){___xxxV_mulMOVO___(  =, += );}
template<class T,int M,int N> inline void
addV_mulMOVO(KVector<T,M>* const pOut,const KMatrix<T,M,N>& A,const KVector<T,N>& B){___xxxV_mulMOVO___( +=, += );}
template<class T,int M,int N> inline void
minV_mulMOVO(KVector<T,M>* const pOut,const KMatrix<T,M,N>& A,const KVector<T,N>& B){___xxxV_mulMOVO___( -=, -= );}



template<class T,int M,int N> inline KVector<T,M>
operator * (const KMatrix<T,M,N>& A,const KVector<T,N>& B)
{KVector<T,M> out; iniV_mulMOVO(&out,A,B); return out;}


#undef ___xxxV_mulMOVO___


/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

#define ___xxxV_mulMTVO___( _eqsini_, _eqs_ ) {                    \
	T* pO = ((T*)pOut);                                            \
	T* pA = ((T*)&A );                                             \
	T* pB = ((T*)&B );                                             \
	{                                                              \
		                    {            *pO _eqsini_ (*pA)*(*pB);}\
		for(int k=1;k<N;k++){pB++;pA+=M; *pO _eqs_    (*pA)*(*pB);}\
	}                                                              \
	for(int i=1;i<M;i++)                                           \
	{                                                              \
		pO++;                                                      \
		pA -= (M*(N-1)-1)   ;                                      \
		pB  = ((T*)&B );                                           \
		                    {            *pO _eqsini_ (*pA)*(*pB);}\
		for(int k=1;k<N;k++){pB++;pA+=M; *pO _eqs_    (*pA)*(*pB);}\
	}                                                              \
}                                                                  \

template<class T,int M,int N> inline void
iniV_mulMTVO(KVector<T,M>* const pOut,const KMatrix<T,N,M>& A, const KVector<T,N>& B){___xxxV_mulMTVO___(  = , += )}
template<class T,int M,int N> inline void
addV_mulMTVO(KVector<T,M>* const pOut,const KMatrix<T,N,M>& A, const KVector<T,N>& B){___xxxV_mulMTVO___( += , += )}
template<class T,int M,int N> inline void
minV_mulMTVO(KVector<T,M>* const pOut,const KMatrix<T,N,M>& A, const KVector<T,N>& B){___xxxV_mulMTVO___( -= , -= )}



#undef ___xxxV_mulMTVO___


/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/
