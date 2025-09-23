//#########################################################################################
template<class T,int M,int N,int K> inline void
iniM_mulMOMO(KMatrix<T,M,N>* const pOut, const KMatrix<T,M,K>& A,const KMatrix<T,K,N>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)
  {
                         (*pOut)(i,j)  = A(i,0) * B(0,j) ;
    for(int k=1;k<K;k++) (*pOut)(i,j) += A(i,k) * B(k,j) ;
  }
}
//#########################################################################################
template<class T,int M,int N,int K> inline void
addM_mulMOMO(KMatrix<T,M,N>* const pOut, const KMatrix<T,M,K>& A,const KMatrix<T,K,N>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)
  {
                         (*pOut)(i,j) += A(i,0) * B(0,j) ;
    for(int k=1;k<K;k++) (*pOut)(i,j) += A(i,k) * B(k,j) ;
  }
}

//#########################################################################################
template<class T,int M,int N,int K> inline void
iniM_mulMTMO(KMatrix<T,M,N>* const pOut, const KMatrix<T,K,M>& A,const KMatrix<T,K,N>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)
  {
                         (*pOut)(i,j)  = A(0,i) * B(0,j) ;
    for(int k=1;k<K;k++) (*pOut)(i,j) += A(k,i) * B(k,j) ;
  }
}
template<class T,int M,int N,int K> inline void
addM_mulMTMO(KMatrix<T,M,N>* const pOut, const KMatrix<T,K,M>& A,const KMatrix<T,K,N>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)for(int k=0;k<K;k++) (*pOut)(i,j) += A(k,i) * B(k,j) ;
}
template<class T,int M,int N,int K> inline void
minM_mulMTMO(KMatrix<T,M,N>* const pOut, const KMatrix<T,K,M>& A,const KMatrix<T,K,N>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)for(int k=0;k<K;k++) (*pOut)(i,j) -= A(k,i) * B(k,j) ;
}
//#########################################################################################

template<class T,int M,int N,int K> inline void
iniM_mulMOMT(KMatrix<T,M,N>* const pOut, const KMatrix<T,M,K>& A,const KMatrix<T,N,K>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)
  {
                         (*pOut)(i,j)  = A(i,0) * B(j,0) ;
    for(int k=1;k<K;k++) (*pOut)(i,j) += A(i,k) * B(j,k) ;
  }
}

template<class T,int M,int N,int K> inline void
addM_mulMOMT(KMatrix<T,M,N>* const pOut, const KMatrix<T,M,K>& A,const KMatrix<T,N,K>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)for(int k=0;k<K;k++) (*pOut)(i,j) += A(i,k) * B(j,k) ;
}

template<class T,int M,int N,int K> inline void
minM_mulMOMT(KMatrix<T,M,N>* const pOut, const KMatrix<T,M,K>& A,const KMatrix<T,N,K>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)for(int k=0;k<K;k++) (*pOut)(i,j) -= A(i,k) * B(j,k) ;
}

//#########################################################################################

#define ___xxxU_mulMOMT___( _eqsini_, _eqs_ ) {             \
  for(int i=0;i<M;i++)for(int j=i;j<M;j++)                  \
  {                                                         \
      T* pO =&((*pOut)(i,j))  ;                             \
                           (*pO) _eqsini_ A(i,0) * B(j,0) ; \
      for(int k=1;k<K;k++) (*pO) _eqs_    A(i,k) * B(j,k) ; \
  }                                                         \
}                                                           \


template<class T,int M,int K> inline void
iniU_mulMOMT(KUSymMat<T,M>* const pOut, const KMatrix<T,M,K>& A,const KMatrix<T,M,K>& B){___xxxU_mulMOMT___(=,+=)}

template<class T,int M,int K> inline void
addU_mulMOMT(KUSymMat<T,M>* const pOut, const KMatrix<T,M,K>& A,const KMatrix<T,M,K>& B){___xxxU_mulMOMT___(+=,+=)}

template<class T,int M,int K> inline void
minU_mulMOMT(KUSymMat<T,M>* const pOut, const KMatrix<T,M,K>& A,const KMatrix<T,M,K>& B){___xxxU_mulMOMT___(-=,-=)}

#undef ___xxxU_mulMOMT___
//#########################################################################################
template<class T,int M,int K> inline
void iniU_mulMTMO(KUSymMat<T,M>* const pOut, const KMatrix<T,K,M>& A,const KMatrix<T,K,M>& B)
{
  for(int i=0;i<M;i++)for(int j=i;j<M;j++)
  {
                         (*pOut)(i,j)  = A(0,i) * B(0,j) ;
    for(int k=1;k<K;k++) (*pOut)(i,j) += A(k,i) * B(k,j) ;
  }
}

template<class T,int M,int K> inline
void addU_mulMTMO(KUSymMat<T,M>* const pOut, const KMatrix<T,K,M>& A,const KMatrix<T,K,M>& B)
{
  for(int i=0;i<M;i++)for(int j=i;j<M;j++)for(int k=0;k<K;k++) (*pOut)(i,j) += A(k,i) * B(k,j) ;
}

template<class T,int M,int K> inline
void minU_mulMTMO(KUSymMat<T,M>* const pOut, const KMatrix<T,K,M>& A,const KMatrix<T,K,M>& B)
{
  for(int i=0;i<M;i++)for(int j=i;j<M;j++)for(int k=0;k<K;k++) (*pOut)(i,j) -= A(k,i) * B(k,j) ;
}

//#########################################################################################

template<class T,int M,int N> inline
void iniM_mulUOMO(KMatrix<T,M,N>* const pOut, const KUSymMat<T,M>& A,const KMatrix<T,M,N>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)
  {
    (*pOut)(i,j) = 0 ;
    for(int ka=0;ka<i;ka++) (*pOut)(i,j) += A(ka,i ) * B(ka,j) ;
    for(int kb=i;kb<M;kb++) (*pOut)(i,j) += A(i ,kb) * B(kb,j) ;
  }
}

template<class T,int M,int N> inline
void iniM_mulMOUO(KMatrix<T,M,N>* const pOut, const KMatrix<T,M,N>& A, const KUSymMat<T,N>& B)
{
  for(int i=0;i<M;i++)for(int j=0;j<N;j++)
  {
    (*pOut)(i,j) = 0 ;
    for(int ka=0;ka<j;ka++) (*pOut)(i,j) += A(i,ka) * B(ka,j) ;
    for(int kb=j;kb<N;kb++) (*pOut)(i,j) += A(i,kb) * B(j,kb) ;
  }
}

template<class T,int M> inline
void iniM_mulUOUO(KMatrix<T,M,M>* const pOut, const KUSymMat<T,M>& A,const KUSymMat<T,M>& B)
{
  for(int i=0;i<M;i++)
  {
    for(int j=0;j<i;j++)
    {
      (*pOut)(i,j) = 0.0 ;
      for(int k=0;k<j;k++)  (*pOut)(i,j) += A(k,i) * B(k,j) ;
      for(int k=j;k<i;k++)  (*pOut)(i,j) += A(k,i) * B(j,k) ;
      for(int k=i;k<M;k++)  (*pOut)(i,j) += A(i,k) * B(j,k) ;
    }
    for(int j=i;j<M;j++)
    {
      (*pOut)(i,j) = 0.0 ;
      for(int k=0;k<i;k++)  (*pOut)(i,j) += A(k,i) * B(k,j) ;
      for(int k=i;k<j;k++)  (*pOut)(i,j) += A(i,k) * B(k,j) ;
      for(int k=j;k<M;k++)  (*pOut)(i,j) += A(i,k) * B(j,k) ;
    }
  }
}

//#########################################################################################

template<class T,int K,int N,int M> inline
KMatrix<T,M,K> operator*(const KMatrix<T,M,N>& a, const KMatrix<T,N,K>& b)
{
  KMatrix<T,M,K> out;
  iniM_mulMOMO(&out , a,b);
  return out;
}

//#########################################################################################
