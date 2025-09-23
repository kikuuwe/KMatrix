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
inline void iniS_outer(double* const pOut, const KVector<double,2>& A, const KVector<double,2>& B)
{
  (*pOut) = A(0)*B(1) - A(1)*B(0);
}
inline double outer(const KVector<double,2>& A, const KVector<double,2>& B)
{
  double out;
  iniS_outer(&out, A, B);
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
inline void ortho(KMatrix<double, 3, 3>* pmat)
{
  const double EPS = 1.0e-8;
  KMatrix<double,3,3>& m = *pmat;
  double a;
  a = 1./(std::max)(EPS,sqrt(m(0, 0)*m(0, 0)+m(1, 0)*m(1, 0)+m(2, 0)*m(2, 0)));
  m(0, 0) *= a;
  m(1, 0) *= a;
  m(2, 0) *= a;
  a = m(0, 0)*m(0, 1)+m(1, 0)*m(1, 1)+m(2, 0)*m(2, 1);
  m(0, 1) -= a*m(0, 0);
  m(1, 1) -= a*m(1, 0);
  m(2, 1) -= a*m(2, 0);
  a = 1./(std::max)(EPS,sqrt(m(0, 1)*m(0, 1)+m(1, 1)*m(1, 1)+m(2, 1)*m(2, 1)));
  m(0, 1) *= a;
  m(1, 1) *= a;
  m(2, 1) *= a;
  m(0, 2) = m(1, 0)*m(2, 1)-m(2, 0)*m(1, 1);
  m(1, 2) = m(2, 0)*m(0, 1)-m(0, 0)*m(2, 1);
  m(2, 2) = m(0, 0)*m(1, 1)-m(1, 0)*m(0, 1);
}
/***************************************************/
inline void iniM_wedgeV(KMatrix<double, 3, 3>* const pM, const KVector<double,3>& omg)
{
  KMatrix<double, 3, 3>& out = *pM;
  out(2, 1) = omg(0);
  out(1, 2) = -omg(0);
  out(0, 2) = omg(1);
  out(2, 0) = -omg(1);
  out(1, 0) = omg(2);
  out(0, 1) = -omg(2);
  out(0, 0) = 0;
  out(1, 1) = 0;
  out(2, 2) = 0;
}
inline KMatrix<double, 3, 3> wedgeV(const KVector<double,3>& omg)
{
  KMatrix<double, 3, 3> out ;
  iniM_wedgeV(&out,omg);
  return out;
}
inline void addM_wedgeV(KMatrix<double,3,3>* const pO, const KVector<double,3>& a)
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
inline void minM_wedgeV(KMatrix<double,3,3>* const pO, const KVector<double,3>& a)
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
inline KVector<double,3> sat(const double& a, const KVector<double,3>& b)
{
   const double EPS = 1.0e-8;
   double babs = (std::max)(EPS,vabs(b));
   if(babs>a) return b*(a/babs);
   return b;
}
inline KVector<double,3> gsat2(const KVector<double,3>& a,const KVector<double,3>& b)
{
  KVector<double,3> out;out.zero();
  out(0) = (std::max)(-a(0),(std::min)(a(0),b(0)));
  out(1) = (std::max)(-a(1),(std::min)(a(1),b(1)));
  out(2) = (std::max)(-a(2),(std::min)(a(2),b(2)));
  return out;
}

inline KVector<double,3> sat3E(double X, const KVector<double,3>& x)
{
  const double EPS = 1.0e-8;
  double ax0 = fabs(x(0));
  double ax1 = fabs(x(1));
  double ax2 = fabs(x(2));
  double axm = (std::max)(ax0,(std::max)(ax1,ax2));
  if(axm>X) return X/axm*x;
  else      return x;
}
/***************************************************/
inline double sinc(const double& a)
{
  const double EPS = 1.0e-10;
  if(fabs(a)<EPS) return 1.0-a*a*(1./6-a*a/120);
  return sin(a)/a;
}
inline double cosc(const double& a)
{
  const double EPS = 1.0e-6;
  if (fabs(a)<EPS) return 0.5-a*a*(1./24-a*a/720);
  return (1.-cos(a))/(a*a);
}
inline double sign(const double& A){ if(A>0)return 1.; if(A<0) return -1.; return 0; }
/***************************************************/

inline void iniM_expV(KMatrix<double,3,3>* const pO, const KVector<double,3>& V)
{
    KMatrix<double,3,3>& M = *pO;
    double a = vabs(V);
    KMatrix<double, 3, 3> K ; iniM_wedgeV (&K, V);
    KMatrix<double, 3, 3> KK; iniM_mulMOMO(&KK, K, K);
    iniX_mulXOS(&M, KK, cosc(a));
    addX_mulXOS(&M, K , sinc(a));
    addM_IS    (&M, 1.         );
    ortho(&M);
}
inline  KMatrix<double,3,3> expV(const KVector<double,3>& V){KMatrix<double,3,3> M;iniM_expV(&M,V);return M;}

/******************************************************************************************/

inline void iniV_logM(KVector<double,3>* const pV,const KMatrix<double,3,3>& R)
{
  KVector<double,3>& u = *pV ;
  double costheta = (trace(R)-1.)*0.5;
  double theta = acos((std::max)(-1.,(std::min)(1.,costheta)));
  if(theta<M_PI/3.) // costheta > 0.5
  {
    u(0) = R(2,1)-R(1,2);
    u(1) = R(0,2)-R(2,0);
    u(2) = R(1,0)-R(0,1);
    u /= (2. * sinc(theta));
  }
  else
  { // Using the idea of Kuo Kan Liang, https://arxiv.org/abs/1810.02999
    KMatrix<double,3,3> RI = R; minM_IS(&RI,1.);
    KVector<double,3>& v0 = *(KVector<double,3>*)&(RI(0,0));
    KVector<double,3>& v1 = *(KVector<double,3>*)&(RI(1,0));
    KVector<double,3>& v2 = *(KVector<double,3>*)&(RI(2,0));
    KVector<double,3> ua = outer(v1,v2); double vabs_ua = vabs(ua);
    KVector<double,3> ub = outer(v2,v0); double vabs_ub = vabs(ub);
    KVector<double,3> uc = outer(v0,v1); double vabs_uc = vabs(uc);
    double vabs_u ;
    if     (vabs_ua>=vabs_ub && vabs_ua>=vabs_uc){u = ua; vabs_u = vabs_ua;}
    else if(vabs_ub>=vabs_ua && vabs_ub>=vabs_uc){u = ub; vabs_u = vabs_ub;}
    else                                         {u = uc; vabs_u = vabs_uc;}
    u  *= (1./vabs_u) * atan2(-trace(wedgeV(u)*R),vabs_u*(trace(R)-1.)) ;
  }
}
  
inline KVector<double,3> logM(const KMatrix<double,3,3>& M){KVector<double,3> V;iniV_logM(&V,M);return V;}

/******************************************************************************************/

inline void init_mulExpVVO(KVector<double,3>* const pOut,const KVector<double,3>& a,const KVector<double,3>& V)
{
    iniV_mulMOVO(pOut,expV(a),V);
}
inline KVector<double,3> mulExpVVO(const KVector<double,3>& a, const KVector<double,3>& V)
{
  KVector<double,3> out;init_mulExpVVO(&out,a,V);return out;
}

/******************************************************************************************/

#endif
