#ifndef Headder_KMatrix_Vector3D
#define Headder_KMatrix_Vector3D

/***************************************************/
inline void iniV_outer(KVector<double,3>* const pOut, const KVector<double, 3>& A, const KVector<double, 3>& B)
{
  const double* pA = ((double*) (&A)) +1 ;
  const double* pB = ((double*) (&B)) +2 ;
  double* pO = (double*) (pOut) ;
  (*pO) =  (*pA) *(*pB); pA++; pB--;
  (*pO) -= (*pA) *(*pB); pO++; pB--;  
  (*pO) =  (*pA) *(*pB); pA = ((double*) (&A))  ;pB = ((double*) (&B)) +2 ;
  (*pO) -= (*pA) *(*pB); pO++; pB--; 
  (*pO) =  (*pA) *(*pB); pA++; pB--;
  (*pO) -= (*pA) *(*pB); 
}
inline KVector<double,3> outer(const KVector<double,3>& A, const KVector<double,3>& B)
{
  KVector<double, 3> out;
  iniV_outer(&out, A, B);
  return out;
}
/***************************************************/

inline void iniS_triprod(double* const pOut, const KVector<double, 3>& A, const KVector<double,3>& B,const KVector<double, 3>& C)
{
  const double* pA = (double*) (&A) ;
  const double* pB = (double*) (&B) ;
  const double* pC = (double*) (&C) ;
  (*pOut)  = (*pA)*(B(1)*C(2)  -  B(2)*C(1) ); pA++; 
  (*pOut) += (*pA)*(B(2)*(*pC) - (*pB)*C(2) ); pA++;
  (*pOut) += (*pA)*((*pB)*C(1) - B(1)*(*pC) ); 
}

inline double vtriprod(const KVector<double, 3>& A, const KVector<double,3>& B,const KVector<double, 3>& C)
{
  double out;
  iniS_triprod(&out, A, B, C);
  return out;
}
/***************************************************/
inline void area_cal(double* const pOut, const KVector<double, 3>& A, const KVector<double,3>& B)
{
  const double* pA = (double*) (&A) ;
  const double* pB = (double*) (&B) ;
  (*pOut)  = ( A(1)*B(2)  -  A(2)*B(1) )*( A(1)*B(2)  -  A(2)*B(1) );
  (*pOut) += ( A(2)*(*pB) - (*pA)*B(2) )*( A(2)*(*pB) - (*pA)*B(2) );
  (*pOut) += ( (*pA)*B(1) - A(1)*(*pB) )*( (*pA)*B(1) - A(1)*(*pB) ); 
  (*pOut)  = sqrt(*pOut)/2;
}

inline double triangle_area(const KVector<double, 3>& A, const KVector<double,3>& B)
{
  double out;
  area_cal(&out, A, B);
  return out;
}

/***************************************************/
inline void addM_VPEquiv(KMatrix<double,3,3>* const pO, const KVector<double,3>& a)
{
  const double* pA = (double*)&a;
        double* pT = (double*)pO;
  pT[2*3+1] += (*pA);
  pT[1*3+2] -= (*pA);pA++;
  pT[0*3+2] += (*pA);
  pT[2*3+0] -= (*pA);pA++;
  pT[1*3+0] += (*pA);
  pT[0*3+1] -= (*pA);

}
inline void minM_VPEquiv(KMatrix<double,3,3>* const pO, const KVector<double,3>& a)
{
  const double* pA = (double*)&a;  
        double* pT = (double*)pO;
  pT[2*3+1] -= (*pA);
  pT[1*3+2] += (*pA);pA++;
  pT[0*3+2] -= (*pA);
  pT[2*3+0] += (*pA);pA++;
  pT[1*3+0] -= (*pA);
  pT[0*3+1] += (*pA);         
}
/***************************************************/
inline void iniS_outer(double* const pOut, const KVector<double,2>& A, const KVector<double,2>& B)
{
  (*pOut) = A(0)*B(1) - A(1)*B(0);
}
/***************************************************/

inline KVector<double,3> gsat(const double& a, const KVector<double,3>& b)
{
    double babs = vabs(b);
    if(babs>a) return b*(a/babs);
    return b;
}
/***************************************************/

inline 
void init_rotvec(KVector<double,3>* const pOut, const KVector<double,3>& a, const KVector<double,3>& b)
{
  iniV_outer(pOut, a, b);
  double tmp   = vabs(*pOut);
  double alpha ; iniS_mulVTVO(&alpha, a, b);
  if(tmp>0)
  {
    alpha = acos(max(-1,min(1,alpha)) );
    *pOut  *= (alpha/tmp);
    return;
  }
  if(alpha<0)
  {
    *pOut = a;
    int mni=-1;
    alpha = 2;
    for(int i=0;i<3;i++){ if(fabs(a(i))<alpha){alpha=fabs(a(i));mni=i;}} 
    (*pOut)(mni)+= 1;
    *pOut -= inner(*pOut,a)*a;
    alpha = vabs(*pOut);
    if(alpha>0) *pOut *= (M_PI/alpha);
    else     {pOut->zero();exit(1);}
    return;
  }
  pOut->zero();
}

inline
void init_RotationMatrix_Partial(KMatrix<double,3,2>* const pMatRot, const KVector<double,3>& aV)
{
  // [MatRot.vec(0) ,MatRot.vec(1), aV] constitutes a right-hand system.
  if( aV(1) != 0 || aV(2) != 0)
  {
    double  L  = vabs(aV);
    double va = aV(0)/L;
    double vb = aV(1)/L;
    double vc = aV(2)/L;
    double  tmp = sqrt( vb * vb + vc * vc );
    double  tmpinv = 1./tmp ;
    (*pMatRot)(0,0) = 0.            ;  (*pMatRot)(0,1) = tmp                ;
    (*pMatRot)(1,0) = - vc * tmpinv ;  (*pMatRot)(1,1) = - va * vb * tmpinv ;
    (*pMatRot)(2,0) =   vb * tmpinv ;  (*pMatRot)(2,1) = - va * vc * tmpinv ;
    return ;
  }
  if(aV(0)>0)
  {
    (*pMatRot)(0,0) = (double)0.;   (*pMatRot)(0,1) = (double)0.;
    (*pMatRot)(1,0) = (double)1.;    (*pMatRot)(1,1) = (double)0.;
    (*pMatRot)(2,0) = (double)0.;    (*pMatRot)(2,1) = (double)1.;
    return ;
  }
  if(aV(0)<0)
  {
    (*pMatRot)(0,0) = (double)0.;     (*pMatRot)(0,1) = (double)0.;
    (*pMatRot)(1,0) = (double)0.;     (*pMatRot)(1,1) = (double)1.;
    (*pMatRot)(2,0) = (double)1.;     (*pMatRot)(2,1) = (double)0.;
    return ;
  }
  pMatRot->zero();
  return ;
}

inline
void init_RotationMatrix(KMatrix<double,3,3>* const pMatRot, const KVector<double,3>& aV)
{
  // [MatRot.vec(0) ,MatRot.vec(1), aV] constitutes a right-hand system.
  if( aV(1) != 0 || aV(2) != 0)
  {
    double  L  = vabs(aV);
    double va = aV(0)/L;
    double vb = aV(1)/L;
    double vc = aV(2)/L;
    double  tmp = sqrt( vb * vb + vc * vc );
    double  tmpinv = 1./tmp ;
    (*pMatRot)(0,0) = 0.            ;  (*pMatRot)(0,1) = tmp                ;(*pMatRot)(0,2) = va ;
    (*pMatRot)(1,0) = - vc * tmpinv ;  (*pMatRot)(1,1) = - va * vb * tmpinv ;(*pMatRot)(1,2) = vb ;
    (*pMatRot)(2,0) =   vb * tmpinv ;  (*pMatRot)(2,1) = - va * vc * tmpinv ;(*pMatRot)(2,2) = vc ;
    return ;
  }
  if(aV(0)>0)
  {
    (*pMatRot)(0,0) = (double)0.;   (*pMatRot)(0,1) = (double)0.;   (*pMatRot)(0,2) = (double)1.;
    (*pMatRot)(1,0) = (double)1.;    (*pMatRot)(1,1) = (double)0.;    (*pMatRot)(1,2) = (double)0.;
    (*pMatRot)(2,0) = (double)0.;    (*pMatRot)(2,1) = (double)1.;    (*pMatRot)(2,2) = (double)0.;
    return ;
  }
  if(aV(0)<0)
  {
    (*pMatRot)(0,0) = (double)0.;   (*pMatRot)(0,1) = (double)0.;   (*pMatRot)(0,2) = (double)-1.;
    (*pMatRot)(1,0) = (double)0.;    (*pMatRot)(1,1) = (double)1.;    (*pMatRot)(1,2) = (double)0.;
    (*pMatRot)(2,0) = (double)1.;    (*pMatRot)(2,1) = (double)0.;    (*pMatRot)(2,2) = (double)0.;
    return ;
  }
  pMatRot->zero();
  return ;
}


inline
void init_rotate(KVector<double,3>* const pr2, const KVector<double,3>& r1, const KVector<double,3>& k, const double& angle)
{
  // SEE YOSHIKAWA'S TEXTBOOK pp.138-139, eq(5.10).,
  KVector<double,3> k_kr1 = k * inner(k,r1);
  *pr2 = k_kr1 + (r1 - k_kr1)* cos(angle) +  outer(k,r1)*sin(angle) ;
}

inline
KVector<double,3> rotate(const KVector<double,3>& r1, const KVector<double,3>& k, const double& angle)
{
  // SEE YOSHIKAWA'S TEXTBOOK pp.138-139, eq(5.10).,
  KVector<double,3> k_kr1 = k * inner(k,r1);
  return   k_kr1 + (r1 - k_kr1)* cos(angle) +  outer(k,r1)*sin(angle) ;
}


inline  KMatrix<double,3,3> axrotmat(const KVector<double,3>& delRa)
{
  KMatrix<double,3,3> out;
  double a = vabs(delRa);
  if(a==0)
  {
    out.zero();out(0,0)=out(1,1)=out(2,2)=1;
    return out;
  }
  KVector<double,3> delR = delRa/ a;
  out(0,0)= (1.-cos(a))*delR(0)*delR(0)+ cos(a) ;
  out(0,1)= (1.-cos(a))*delR(0)*delR(1)-(sin(a))*delR(2);
  out(0,2)= (1.-cos(a))*delR(0)*delR(2)+(sin(a))*delR(1);
  out(1,0)= (1.-cos(a))*delR(1)*delR(0)+(sin(a))*delR(2);
  out(1,1)= cos(a)+(1.-cos(a))*delR(1)*delR(1);
  out(1,2)= (1.-cos(a))*delR(1)*delR(2)-(sin(a))*delR(0);
  out(2,0)= (1.-cos(a))*delR(0)*delR(2)-(sin(a))*delR(1);
  out(2,1)= (1.-cos(a))*delR(2)*delR(1)+(sin(a))*delR(0);
  out(2,2)= cos(a)+(1.-cos(a))*delR(2)*delR(2);
  return out;
}

inline KVector<double,3> rotvec(const KMatrix<double,3,3>& A)
{
  KVector<double,3> out;
  out(0) = (A(2,1)-A(1,2))/2;
  out(1) = (A(0,2)-A(2,0))/2;
  out(2) = (A(1,0)-A(0,1))/2;
  double sa = vabs(out);
  double ca = (A(0,0)+A(1,1)+A(2,2)-1.)/2;
  double alpha = atan2(sa,ca);
  if(sa>0) out *= (alpha/sa) ;
  else out.zero();
  return out;
}
inline KVector<double,3> axrot(const KVector<double,3>& delRa,const KVector<double,3>& b)
{
  double a = vabs(delRa);
  if(a==0) return b;
  KVector<double,3> out;
  iniV_mulMOVO(&out,axrotmat(delRa),b);
  return out;
}


#endif

