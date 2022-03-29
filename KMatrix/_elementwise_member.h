//##########################################
#ifndef __USYM_MODE__
//##########################################
#define macro_ope_zero(__enzan__){        \
  T* pO   = (T*)this  ;                   \
  for(int i=0;i<_NUM_-1;i++){             \
    __enzan__ ;                           \
    pO++;                                 \
  }                                       \
  {                                       \
    __enzan__ ;                           \
  }                                       \
}
#define macro_ope_one(__enzan__){         \
  T* pO   = (T*)this  ;                   \
  const T* pA   = (T*)&A    ;             \
  for(int i=0;i<_NUM_-1;i++){             \
    __enzan__ ;                           \
    pO++;pA++;                            \
  }                                       \
  {                                       \
    __enzan__ ;                           \
  }                                       \
}
#define macro_ope_two(__enzan__){         \
        T* pO   = (T*)this  ;             \
  const T* pA   = (T*)&A    ;             \
  const T* pB   = (T*)&B    ;             \
  for(int i=0;i<_NUM_-1;i++){             \
    __enzan__ ;                           \
    pO++;pA++;pB++;                       \
  }                                       \
  {                                       \
    __enzan__ ;                           \
  }                                       \
}
//##########################################
#else
//##########################################
#define macro_ope_zero(__enzan__){  \
  T* pO   = (T*)this  ;             \
  for(int i =0;i<_UNUM_-1;i++){     \
    for(int j=i;j<_UNUM_-1;j++){    \
      __enzan__ ;                   \
      pO++;                         \
    }                               \
    __enzan__ ;                     \
    pO+=i+2;                        \
  }                                 \
  {                                 \
    __enzan__ ;                     \
  }                                 \
}
#define macro_ope_one(__enzan__){   \
        T* pO   = (T*)this  ;       \
  const T* pA   = (T*)&A    ;       \
  for(int i =0;i<_UNUM_-1;i++){     \
    for(int j=i;j<_UNUM_-1;j++){    \
      __enzan__ ;                   \
      pO++;pA++;                    \
    }                               \
    __enzan__ ;                     \
    pO+=i+2;pA+=i+2;                \
  }                                 \
  {                                 \
    __enzan__ ;                     \
  }                                 \
}
#define macro_ope_two(__enzan__){   \
        T* pO   = (T*)this  ;       \
  const T* pA   = (T*)&A    ;       \
  const T* pB   = (T*)&B    ;       \
  for(int i =0;i<_UNUM_-1;i++){     \
    for(int j=i;j<_UNUM_-1;j++){    \
      __enzan__ ;                   \
      pO++;pA++;pB++;               \
    }                               \
    __enzan__ ;                     \
    pO+=i+2;pA+=i+2;pB+=i+2;        \
  }                                 \
  {                                 \
    __enzan__ ;                     \
  }                                 \
}
//##########################################
#endif
//##########################################

_TYPE_ operator =(const _TYPE_& A)
{
	macro_ope_one({
		 *pO = (*pA) ;
	});
	return *this;
}
_TYPE_ operator +=(const _TYPE_& A)
{
	macro_ope_one({
		 *pO += (*pA) ;
	});
	return *this;
}
_TYPE_ operator -=(const _TYPE_& A)
{
	macro_ope_one({
		 *pO -= (*pA) ;
	});
	return *this;
}
_TYPE_ operator *=(const T& ko)
{
	macro_ope_zero({
		 *pO *= ko ;
	});
	return *this;
}
_TYPE_ operator /=(const T& ko)
{
	(*this) *= (1./ko);
	return *this;
}

#ifndef __USYM_MODE__
_TYPE_ operator -() const
{
	_TYPE_ A ;
  const T* pO   = (T*)this  ;  
        T* pA   = (T*)&A    ;  
  for(int i=0;i<(_NUM_)-1;i++){
    (*pA) = -(*pO) ;           
    pO++;pA++;                 
  }                            
  {                            
    (*pA) = -(*pO) ;           
  }                            
	return A;
}
#endif
//############################################################################################



#undef macro_ope_zero
#undef macro_ope_one
#undef macro_ope_two
