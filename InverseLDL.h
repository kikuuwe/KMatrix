#ifndef Headder_KMatrix_InverseLDL
#define Headder_KMatrix_InverseLDL

/*************************************************************/
namespace KNSforLDInverse
{
  template<int N> struct KDummySize {};
  template<int N,int M> inline
  double LD_Decompose_Here(KUSymMat<double,N>* pL, KDummySize<M>, int* idx=NULL) 
  {
    #define __swap(x,y) {{tmp = (x); (x) = (y); (y) = tmp; }}
    for(int k=0; k<M-1;k++)
    {
      double tmp  =  (*pL)( k , k);
      int imax = k;
      for(int i=k+1;i<M;i++)
      {
        if(tmp < (*pL)( i , i))
        {
          tmp = (*pL)( i , i) ;
          imax = i;
        }
      }
      if(idx!=NULL) idx[k] = imax;
      if(imax!=k)
      {
        double tmp;
        int i;
        for(i=0     ;i<k   ;i++)__swap((*pL)(   k,   i), (*pL)(imax, i  ));
        for(i=k+1   ;i<imax;i++)__swap((*pL)(   k,   i), (*pL)(   i,imax));
        for(i=imax+1;i<M   ;i++)__swap((*pL)(   k,   i), (*pL)(imax,   i));
                                __swap((*pL)(   k,   k), (*pL)(imax,imax));
      }
      if((*pL)(  k , k ) <=0.  ) return (*pL)(  k , k );
      tmp = 1./(*pL)(  k , k ) ;
      int j;
      for(j=k+1;j<M;j++) (*pL)(  j,  k) = (*pL)(  k,j) * tmp ;  // k<j
      for(j=k+1;j<M;j++)
      {
        tmp = (*pL)( j ,k )  * (*pL)( k, k) ;
        for(int i=j;i<M;i++) (*pL)(j ,i) -= (*pL)(i,k) * tmp ; // k<j<i
      }
    }
    if((*pL)( M-1 , M-1)<=0.){return (*pL)( M-1 , M-1);}
    return 1.0;
    #undef __swap
  }
  /*************************************************************/
  template<int N,int M> inline
  void LD_inverse_and_Composite_Here(KUSymMat<double,N>* pL, const int * idx, KDummySize<M>)
  {
    double tmp;
    int    idxk;
    #define __swap(x,y) {{tmp = (x); (x) = (y); (y) = tmp; }}
  // processing transpose-inverse of L at upper triagonal region.
  // Diagonal elements are inversed on-the-spot.
  // U = L\inv\T
  // A = LDL\T
  // A\inv = L\inv\T D\inv L\inv = U D\inv U\T
    for(int i=0;i<M;i++)
    {
      for(int j=i+1;j<M;j++)  // i < k < j
      {
                                   (*pL)(i,j) = -(*pL)(j,i) ;
        for(int k= i +1;k < j;k++) (*pL)(i,j) -= (*pL)(i,k) * (*pL)(j,k) ;
      }
      (*pL)(i,i) = 1.0/(*pL)(i,i);
    }
    // processing (U D\inv)\T in lower triagonal region  i < j
    for(int i=0;i<M;i++)for(int j=i+1;j<M;j++) (*pL)(j, i)= (*pL)(i, j) * (*pL)(j, j);
    // processing U * (U D\inv)\T in upper triagonal region
    for(int j=0;j<M;j++)  //  j < i < k
    {
      for(int k=j+1;k<M;k++)(*pL)(j,j) += (*pL)(k,j)*(*pL)(j,k) ;
      for(int i=j+1;i<M;i++)
      {
        tmp =  0.;
        for(int k=i;k<M;k++) tmp += (*pL)(k,i)*(*pL)(j,k);
        (*pL)(j,i) = tmp;
      }
    }
    for(int k=M-2;k>=0;k--)
    {
      idxk = idx[k];
      if(idxk!=k)
      {
        for(int i=0     ;i<k   ;i++) __swap((*pL)( i , k ), (*pL)(  i ,idxk));
        for(int i=k+1   ;i<idxk;i++) __swap((*pL)( k , i ), (*pL)(  i ,idxk));
        for(int i=idxk+1;i<M   ;i++) __swap((*pL)( k , i ), (*pL)(idxk, i  ));
                                     __swap((*pL)( k , k ), (*pL)(idxk,idxk));
      }
    }
    #undef __swap
  }
}


//#####################################################################X
template<int N> inline
double invU(KUSymMat<double,N>* const pO)
{
  KNSforLDInverse::KDummySize<N> ds;
  int idx[N];
  double out = KNSforLDInverse::LD_Decompose_Here(pO,ds,idx);
  if(out<=0.0) return out;
  KNSforLDInverse::LD_inverse_and_Composite_Here(pO,idx,ds);
  return out;
}
//#####################################################################X
template<int N> inline
double iniU_invUO(KUSymMat<double,N>* const pOut, const KUSymMat<double,N>& A)
{
  *pOut = A;
  return invU(pOut);
}
//#####################################################################X

template<int N> inline
double iniU_invUO_11(KUSymMat<double,N>* const pOut , const KUSymMat<double,N>& A)
{
	double               aa =A(0,0);
	KVector <double,N-1> bb;
	KUSymMat<double,N-1> CC;
	for(int i=0;i<N-1;i++)                      bb(i)  =A(0,  i+1);
	for(int i=0;i<N-1;i++)for(int j=i;j<N-1;j++)CC(i,j)=A(i+1,j+1);
	addU_squVOSVT(&CC, bb, -1./aa );
	invU(&CC);
	KUSymMat<double,N> CCE; CCE.zero(); CCE(0,0) = 1./aa;
	for(int i=0;i<N-1;i++)for(int j=i;j<N-1;j++)CCE(i+1,j+1)=CC(i,j);
	KMatrix<double,N,N> MM; MM.zero();
	for(int i=0;i<N  ;i++)MM(i  ,i)=1;
	for(int i=0;i<N-1;i++)MM(i+1,0)=-bb(i)/aa;
	iniU_quadMTUOMO(pOut, CCE, MM );
	return 1;
}


#endif

