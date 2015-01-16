//#########################################################################################

template<class T, int N>struct KUSymMat
{
	T elem[N*N];

	void  zero    () {memset(this,0,sizeof(*this));}
	T& operator()(int i, int j)const {return *(((T*)this)+ (i * N + j));}

	KMatrix<T,N,N>&  asMatrix()const {return *((KMatrix<T,N,N>*)this)  ;}
	void  restore(){for(int i=0;i<N;i++)for(int j=i+1;j<N;j++){(*this)(j,i)  = (*this)(i,j) ;}};

	#define __USYM_MODE__
	#define _TYPE_ KUSymMat<T,N>
	#define _UNUM_  (N)
	#include "KMatrix/_elementwise_member.h"
	#undef _TYPE_
	#undef _UNUM_
	#undef __USYM_MODE__
};


#define __USYM_MODE__
#define _UNUM_          (N)
#define _DECTEMPLATE_  template<class T,int N> 
#define _TYPE_         KUSymMat<T,N>
#include "KMatrix/_elementwise_global.h"
#undef _DECTEMPLATE_
#undef _TYPE_
#undef _UNUM_
#undef __USYM_MODE__

//#########################################################################################

