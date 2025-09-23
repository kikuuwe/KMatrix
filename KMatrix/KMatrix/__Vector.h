//#########################################################################################
template<class T, int N>struct KVector
{
  T elem[N];
  inline void zero    ()            {memset(this,0,sizeof(*this));}
  inline T&   operator()(int i)const{return *(((T*)this)+(i));}
  #define _TYPE_ KVector<T,N>
  #define _NUM_  (N)
  #include "_elementwise_member.h"
  #undef _TYPE_
  #undef _NUM_
};

#define _TYPE_         KVector<T,N>
#define _NUM_          (N)
#define _DECTEMPLATE_  template<class T,int N> 
#include "_elementwise_global.h"
#undef _DECTEMPLATE_
#undef _TYPE_
#undef _NUM_

//#########################################################################################

#define __xxxS_mulVTVO__(_ini_eqs_, _eqs_ ) {                        \
  T* pA = (T*)&A;                                                    \
  T* pB = (T*)&B;                                                    \
                      {           *(pOut) _ini_eqs_ (*pA) * (*pB) ;} \
  for(int i=1;i<N;i++){pA++;pB++; *(pOut) _eqs_     (*pA) * (*pB) ;} \
}                                                                    \

template<class T, int N> inline void
iniS_mulVTVO(T* const pOut, const KVector<T,N>& A, const KVector<T,N>& B){ __xxxS_mulVTVO__( = , += );}
template<class T, int N> inline void
addS_mulVTVO(T* const pOut, const KVector<T,N>& A, const KVector<T,N>& B){ __xxxS_mulVTVO__( += , += );}
template<class T, int N> inline void
minS_mulVTVO(T* const pOut, const KVector<T,N>& A, const KVector<T,N>& B){ __xxxS_mulVTVO__( -= , -= );}

#undef __xxxS_mulVTVO__

//#########################################################################################


#define __xxxS_squVTVO__(_ini_eqs_, _eqs_ ) {                   \
  T* pA = (T*)&A;                                               \
                      {       *(pOut) _ini_eqs_ (*pA) * (*pA) ;}\
  for(int i=1;i<N;i++){pA++;  *(pOut) _eqs_     (*pA) * (*pA) ;}\
}                                                               \


template<class T, int N> inline void 
iniS_squVTVO(T* const pOut, const KVector<T,N>& A){__xxxS_squVTVO__( =,+=)}
template<class T, int N> inline void 
addS_squVTVO(T* const pOut, const KVector<T,N>& A){__xxxS_squVTVO__(+=,+=)}
template<class T, int N> inline void 
minS_squVTVO(T* const pOut, const KVector<T,N>& A){__xxxS_squVTVO__(-=,-=)}

#undef __xxxS_squVTVO__


//#########################################################################################

template<class T, int N> inline
T inner(const KVector<T,N>& A,const KVector<T,N>& B){T a;iniS_mulVTVO(&a,A,B);return a;}

template<class T, int N> inline
T vabs (const KVector<T,N>& A){T a;iniS_squVTVO(&a,A);return (T)sqrt(a);}


//#########################################################################################
