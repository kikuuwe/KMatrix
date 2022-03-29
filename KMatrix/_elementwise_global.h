//##########################################
#ifndef __USYM_MODE__
//##########################################
#define macro_ope_zero(__enzan__){        \
  T* pO   = (T*)pOut  ;                   \
  for(int i=0;i<_NUM_-1;i++){             \
    __enzan__ ;                           \
    pO++;                                 \
  }                                       \
  {                                       \
    __enzan__ ;                           \
  }                                       \
}
#define macro_ope_one(__enzan__){         \
        T* pO   = (      T*)pOut  ;       \
  const T* pA   = (const T*)&A    ;       \
  for(int i=0;i<_NUM_-1;i++){             \
    __enzan__ ;                           \
    pO++;pA++;                            \
  }                                       \
  {                                       \
    __enzan__ ;                           \
  }                                       \
}
#define macro_ope_two(__enzan__){         \
        T* pO   = (      T*)pOut  ;       \
  const T* pA   = (const T*)&A    ;       \
  const T* pB   = (const T*)&B    ;       \
  for(int i=0;i<_NUM_-1;i++){             \
    __enzan__ ;                           \
    pO++;pA++;pB++;                       \
  }                                       \
  {                                       \
    __enzan__ ;                           \
  }                                       \
}
#define macro_ope_three(__enzan__){       \
        T* pO   = (      T*)pOut  ;       \
  const T* pA   = (const T*)&A    ;       \
  const T* pB   = (const T*)&B    ;       \
  const T* pC   = (const T*)&C    ;       \
  for(int i=0;i<_NUM_-1;i++){             \
    __enzan__ ;                           \
    pO++;pA++;pB++;pC++;                  \
  }                                       \
  {                                       \
    __enzan__ ;                           \
  }                                       \
}
//##########################################
#else
//##########################################
#define macro_ope_zero(__enzan__){  \
  T* pO   = (T*)pOut  ;             \
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
  T* pO   = (T*)pOut  ;             \
  T* pA   = (T*)&A    ;             \
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
  T* pO   = (T*)pOut  ;             \
  T* pA   = (T*)&A    ;             \
  T* pB   = (T*)&B    ;             \
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
#define macro_ope_three(__enzan__){ \
        T* pO   = (      T*)pOut  ; \
  const T* pA   = (const T*)&A    ; \
  const T* pB   = (const T*)&B    ; \
  const T* pC   = (const T*)&C    ; \
  for(int i =0;i<_UNUM_-1;i++){     \
    for(int j=i;j<_UNUM_-1;j++){    \
      __enzan__ ;                   \
      pO++;pA++;pB++;pC++;          \
    }                               \
    __enzan__ ;                     \
    pO+=i+2;pA+=i+2;pB+=i+2;pC+=i+2;\
  }                                 \
  {                                 \
    __enzan__ ;                     \
  }                                 \
}
//##########################################
#endif
//##########################################
_DECTEMPLATE_ inline
void iniX_addXOXO(_TYPE_* const pOut, const _TYPE_& A, const _TYPE_& B)
{
	macro_ope_two({
		 *pO = (*pA) + (*pB);
	});
}
_DECTEMPLATE_ inline
void iniX_minXOXO(_TYPE_* const pOut, const _TYPE_& A, const _TYPE_& B)
{
	macro_ope_two({
		 *pO = (*pA) - (*pB);
	});
}
_DECTEMPLATE_ inline
void iniX_addXOXOXO(_TYPE_* const pOut, const _TYPE_& A, const _TYPE_& B, const _TYPE_& C)
{
	macro_ope_three({
		 *pO = (*pA) + (*pB) + (*pC);
	});
}


_DECTEMPLATE_ inline
void iniX_mulXOS(_TYPE_* const pOut, const _TYPE_& A, const T& ka )
{
	macro_ope_one({
		 *pO = (*pA)*ka ;
	});
}
_DECTEMPLATE_ inline
void iniX_divXOS(_TYPE_* const pOut, const _TYPE_& A, const T& ka )        
{
	iniX_mulXOS(pOut,A, 1./ka);
}

_DECTEMPLATE_ inline
void iniX_mulXX_cmpntwise(_TYPE_* const pOut, const _TYPE_& A, const _TYPE_& B )
{
	macro_ope_two({
		 *pO = (*pA)*(*pB) ;
	});
}
//##########################################
_DECTEMPLATE_ inline
void addX_XO(_TYPE_* const pOut, const _TYPE_& A)
{
	macro_ope_one({
		 *pO += (*pA);
	});
}
_DECTEMPLATE_ inline
void minX_XO(_TYPE_* const pOut, const _TYPE_& A)
{
	macro_ope_one({
		 *pO -= (*pA);
	});
}
//##########################################
_DECTEMPLATE_ inline
void iniX_XO(_TYPE_* const pOut, const _TYPE_& A)
{
	macro_ope_one({
		 *pO = (*pA);
	});
}
_DECTEMPLATE_ inline
void negX_XO(_TYPE_* const pOut, const _TYPE_& A)
{
	macro_ope_one({
		 *pO =- (*pA);
	});
}
//##########################################
_DECTEMPLATE_ inline
void mulX_S(_TYPE_* const pOut, const T& ka)              
{
	macro_ope_zero({
		 *pO *=ka ;
	});
}
// mulX_S_addX_XO
_DECTEMPLATE_ inline
void pile_times_plus(_TYPE_* const pOut, T ko, const _TYPE_& A)      // Mul_AddProdOf
{
	macro_ope_one({
		 *pO *=ko    ;
		 *pO +=(*pA) ;
	});
}

_DECTEMPLATE_ inline
void addX_mulXOS(_TYPE_* const pOut, const _TYPE_& A,  const T& ka)  
{
	macro_ope_one({
		 *pO +=(*pA)*ka ;
	});
}
_DECTEMPLATE_ inline
void minX_mulXOS(_TYPE_* const pOut, const _TYPE_& A,  const T& ka) 
{
	macro_ope_one({
		 *pO -=(*pA)*ka ;
	});
}
//##########################################
_DECTEMPLATE_ inline
void init_wplus_wplus(_TYPE_* const pOut, const _TYPE_& A, const T& ka, const _TYPE_& B, const T& kb)  
{
	macro_ope_two({
		 *pO = (*pA)*ka + (*pB)*kb;
	});
}

_DECTEMPLATE_ inline
void init_plus_wplus(_TYPE_* const pOut,  const _TYPE_& A, const _TYPE_& B, const T& kb)  // ___WeightedSumOf
{
	macro_ope_two({
		 *pO = (*pA) + (*pB)*kb ;
	});
}

// mulX_S_addX_mulXOS
_DECTEMPLATE_ inline
void pile_times_wplus(_TYPE_* const pOut, T ko, const _TYPE_& A,  const T& ka)      // Mul_AddProdOf
{
	macro_ope_one({
		 *pO *=ko       ;
		 *pO +=(*pA)*ka ;
	});
}
//////////////////////////////////////

//############################################################################################
_DECTEMPLATE_ inline
_TYPE_ operator+  (const _TYPE_& a, const _TYPE_& b) 
{
	_TYPE_ out;
	::iniX_addXOXO(&out,a, b);
	return out;
}
//-----------------------------------------
_DECTEMPLATE_ inline
_TYPE_ operator-  (const _TYPE_& a, const _TYPE_& b) 
{
	_TYPE_ out;
	::iniX_minXOXO(&out,a, b);
	return out;
}
//-----------------------------------------
_DECTEMPLATE_ inline
_TYPE_ operator*  (const _TYPE_& a, const T&  b) 
{
	_TYPE_ out;
	::iniX_mulXOS(&out,a, b);
	return out;
}
//-----------------------------------------
_DECTEMPLATE_ inline
_TYPE_ operator* (const T&  b, const _TYPE_& a)          
{
  _TYPE_ out;
  iniX_mulXOS(&out, a, b);
  return out;
}
//-----------------------------------------
_DECTEMPLATE_ inline
_TYPE_ operator/  (const _TYPE_& a, const T& b) 
{
	_TYPE_ out;
	::iniX_mulXOS(&out,a, 1./b);
	return out;
}
//############################################################################################

//_TYPE_ operator- () const    
//{
//  _TYPE_ out;
//  ::iniX_mulXOS(&out,*this, (T)(-1.));
//  return out;
//}


#undef macro_ope_zero
#undef macro_ope_one
#undef macro_ope_two
#undef macro_ope_three
