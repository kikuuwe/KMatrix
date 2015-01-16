
//#########################################################################################
template<class T,int M,int N>struct KMatrix
{
  T elem[M*N];
  void  zero   ()                 {memset(this,0,sizeof(*this));}
  T& operator()(int i, int j)const{return *(((T*)this)+ (i * N + j));}
  #define _TYPE_ KMatrix<T,M,N>
  #define _NUM_  (M*N)
  #include "KMatrix/_elementwise_member.h"
  #undef _TYPE_
  #undef _NUM_
};

#define _NUM_          (M*N)
#define _DECTEMPLATE_  template<class T,int M,int N> 
#define _TYPE_         KMatrix<T,M,N>
#include "KMatrix/_elementwise_global.h"
#undef _DECTEMPLATE_
#undef _TYPE_
#undef _NUM_
//#########################################################################################
template<class T,int M> inline
void iniS_trace(T *pOut,const KMatrix<T,M,M>& A)
{
  *pOut = A(0,0); for(int i=1;i<M;i++) (*pOut)+=A(i,i);
}

template<class T, int M> inline
T trace(const KMatrix<T,M,M>& A)
{
  T out ;
  iniS_trace(&out,A);
  return out;
}

//#########################################################################################
//#########################################################################################
template<class T, int M> inline void
addM_IS(KMatrix<T,M,M>* const pOut, const T& val)
{
  T* pO = (T*)pOut;
                      {              (*pO) += val;}
  for(int i=1;i<M;i++){ pO += (M+1); (*pO) += val;}
}
template<class T, int M> inline void
minM_IS(KMatrix<T,M,M>* const pOut, const T& val)
{
  T* pO = (T*)pOut;
                      {              (*pO) -= val;}
  for(int i=1;i<M;i++){ pO += (M+1); (*pO) -= val;}
}
//#########################################################################################
//#########################################################################################

#define __xxxM_MT__( _eqs_ ) {                            \
  T* pO = (T*)pOut;                                       \
  T* pA = (T*)&A;                                         \
  for(int i=0;i<M-1;i++)                                  \
  {                                                       \
    for(int j=0;j<N-1;j++){(*pO) _eqs_ (*pA);pO++;pA+=M;} \
                          {(*pO) _eqs_ (*pA);pO++;}       \
    pA-=(M*N-M-1);                                        \
  }                                                       \
  {                                                       \
    for(int j=0;j<N-1;j++){(*pO) _eqs_ (*pA);pO++;pA+=M;} \
                          {(*pO) _eqs_ (*pA);}            \
  }                                                       \
}                                                         \

//#########################################################################################

template<class T,int M,int N> inline
void addM_MT(KMatrix<T,M,N>* const pOut, const KMatrix<T,N,M>& A)
{
  __xxxM_MT__(+=)
}
template<class T,int M,int N> inline
void iniM_MT(KMatrix<T,M,N>* const pOut, const KMatrix<T,N,M>& A)
{
  __xxxM_MT__(=)
}
template<class T,int M,int N> inline
void negM_MT(KMatrix<T,M,N>* const pOut, const KMatrix<T,N,M>& A)
{
  __xxxM_MT__(=-)
}

template<class T,int M,int N> inline
void minM_MT(KMatrix<T,M,N>* const pOut, const KMatrix<T,N,M>& A)
{
  __xxxM_MT__(-=)
}

#undef __xxxM_MT__

//#########################################################################################

#undef macro_check

