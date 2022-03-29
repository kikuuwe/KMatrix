template<class T,int M> inline
void addM_UO(KMatrix<T,M,M>* const pOut, const KUSymMat<T,M>& A)
{
  const T * pA  = (const T*) &A   ;
        T * pO1 = (      T*) pOut ;
        T * pO2 = (      T*) pOut ;
  {
                                             (*pO1)+=(*pA);
    for(int j=0+1;j<M;j++){pO1++;pO2+=M;pA++;(*pO1)+=(*pA);(*pO2)+=(*pA);}
  }                        
  for(int i=1;i<M;i++)     
  {                        
    pA +=(1+i) ;pO1+=(1+i);pO2 =((T*)pOut)+(M+1)*i;  
                                             (*pO1)+=(*pA);
    for(int j=i+1;j<M;j++){pO1++;pO2+=M;pA++;(*pO1)+=(*pA);(*pO2)+=(*pA);}
  }
}
/************************************************************************/

template<class T,int M> inline
void minU_MO(KUSymMat<T,M>* const pOut, const KMatrix<T,M,M>& A)
{
        T * pO = (      T*) pOut ;
  const T * pA = (const T*) &A   ;
  {
                                     (*pO) -= (*pA);
    for(int j=0+1;j<M;j++){pO++;pA++;(*pO) -= (*pA);}
  }
  for(int i=1;i<M;i++)
  {
    pO+=(1+i); pA +=(1+i) ;
                                     (*pO) -= (*pA);
    for(int j=i+1;j<M;j++){pO++;pA++;(*pO) -= (*pA);}
  }
}

template<class T,int M> inline
void minU_MT(KUSymMat<T,M>* const pOut, const KMatrix<T,M,M>& A)
{
        T * pO = (      T*) pOut ;
  const T * pA = (const T*) &A   ;
  {
                      (*pO) -= (*pA);
    for(int j=  1;j<M;j++){pO++;pA+=M;(*pO) -= (*pA);}
  }
  for(int i=1;i<M;i++)
  {
    pO+=(1+i); pA =((const double*) &A)+ ((i)*(M+1)) ;
                                        (*pO) -= (*pA);
    for(int j=i+1;j<M;j++){pO++;pA+=M;(*pO) -= (*pA);}
  }
}

/************************************************************************/

template<class T, int N, int M> inline
void iniU_squMOMT(KUSymMat<T,M> * const pO, const KMatrix<T,M,N>& A)
{
  for(int i=0;i<M;i++)for(int j=i;j<M;j++)
  {
    (*pO)(i,j) = 0 ;
    for(int k=0;k<N;k++) (*pO)(i,j) += A(i,k) * A(j,k) ;
  }
}

template<class T, int M, int N> inline
void iniU_squMTMO(KUSymMat<T,M> * const pOut, const KMatrix<T,N,M>& A)
{
  double *pAT = (double*)&A;
  double *pAO = (double*)&A;
  for(int i=0;i<M;i++)
  {
    double *pO  = (double*)&((*pOut)(i,i));
    for(int j=i;j<M;j++)
    {
      pAT = (double*)&(A(0,i));
      pAO = (double*)&(A(0,j));
                                          (*pO)  = (*pAT) * (*pAO) ;
      for(int k=1;k<N;k++){pAT+=M;pAO+=M; (*pO) += (*pAT) * (*pAO) ;}
      pO ++;
    }
  }
}
/************************************************************************/

template<int M,int N> inline
void iniU_quadMOUOMT(KUSymMat<double,M>* const pOut, const KUSymMat<double,N>& A, const KMatrix<double,M,N>& B)
{
  KMatrix<double,M,N> tBA  ;
  iniM_mulMOUO(&tBA  , B ,    A  );
  iniU_mulMOMT(pOut  ,  tBA , B  );
}

template<int M,int N> inline
void addU_quadMOUOMT(KUSymMat<double,M>* const pOut, const KUSymMat<double,N>& A, const KMatrix<double,M,N>& B)
{
  KMatrix<double,M,N> tBA  ;
  iniM_mulMOUO(&tBA  , B ,    A  );
  addU_mulMOMT(pOut  ,  tBA , B  );
}

template<int M,int N> inline
void iniU_quadMTUOMO(KUSymMat<double,M>* const pOut, const KUSymMat<double,N>& A, const KMatrix<double,N,M>& B)
{
  KMatrix<double,N,M> tAB  ;
  iniM_mulUOMO(&tAB  , A ,    B  );
  iniU_mulMTMO(pOut   , B ,  tAB );
}
