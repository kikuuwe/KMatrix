#ifndef KQUARTERNION_H_INCLUDED
#define KQUARTERNION_H_INCLUDED
/**********************************************/
typedef KVector<double,3>   vec3;
typedef KMatrix<double,3,3> mat3;
/**********************************************/
struct quat
{
	double w,x,y,z; 
	void init(){w=1;x=y=z=0;};
	double norm()const{return std::sqrt(w*w+x*x+y*y+z*z);}
	void normalize()
	{
		double r = norm();
		if(r<=0) r=0; 
		else r= 1./r;
		w *= r;
		x *= r;
		y *= r;
		z *= r;
	}
};
struct vecvec
{
	double elm[6] ;
	void zero(){for(int i=0;i<6;i++)elm[i]=0;} 
	double& operator()(int i)const{return *(((double*)this)+i);} 
	vecvec& operator*=(const double& a){for(int i=0;i<6;i++)elm[i]*=a       ;return *this;} 
	vecvec& operator+=(const vecvec& a){for(int i=0;i<6;i++)elm[i]+=a.elm[i];return *this;} 
	vecvec& operator-=(const vecvec& a){for(int i=0;i<6;i++)elm[i]-=a.elm[i];return *this;} 
	vec3& sv0()const{return *((vec3*)this);}
	vec3& sv1()const{return *((vec3*)(&elm[3]));}
	KVector<double,6>& asvec(){return *((KVector<double,6>*)this);}
};
struct vecqua
{
	double elm[7]; 
	vec3& vec()const {return *((vec3*)this) ;}
	quat& qua()const {return *((quat*)(((double*)this)+3)) ;}
	double& operator()(int i)const{return *(((double*)this)+i);} 
	vecqua& operator +=(const vecvec& b);
	void  init(){vec().zero();qua().init();}
	void  normalize(){qua().normalize();}
}; 
/**********************************************/
inline quat normalize(const quat& q)
{
	quat out = q;
	out.normalize();
	return out;
}
/**********************************************/
inline quat inv(const quat& a)
{
	double r = sqrt(a.w*a.w+a.x*a.x+a.y*a.y+a.z*a.z);
	if(r==0){ return quat({1,0,0,0});}
	r = 1./r;
	return quat({a.w*r,-a.x*r,-a.y*r,-a.z*r});
}
/**********************************************/
inline void nearer(const quat& a, quat* pb)
{
	if(a.w*pb->w+a.x*pb->x+a.y*pb->y+a.z*pb->z>=0)return ;
	pb->w *=-1.;
	pb->x *=-1.;
	pb->y *=-1.;
	pb->z *=-1.;
}
inline void nearer(quat* pb)
{
	if(pb->w>=0)return ;
	pb->w *=-1.;
	pb->x *=-1.;
	pb->y *=-1.;
	pb->z *=-1.;
}
inline quat nearer(const quat& a)
{
	if(a.w>=0)return a ;
	return quat({-a.w,-a.x,-a.y,-a.z});
}

inline int nearer(const quat& a, const quat& b, const quat& c)
{
	double ab = (a.w-b.w)*(a.w-b.w)+(a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z);
	double bc = (b.w-c.w)*(b.w-c.w)+(b.x-c.x)*(b.x-c.x)+(b.y-c.y)*(b.y-c.y)+(b.z-c.z)*(b.z-c.z);
	if(ab<bc) return 1;
	else      return 0;
}
/**********************************************/
inline quat operator*(const quat& a, const quat& b)
{
	quat c; 
	c.w = a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z;
	c.x = a.x*b.w+a.w*b.x-a.z*b.y+a.y*b.z;
	c.y = a.y*b.w+a.z*b.x+a.w*b.y-a.x*b.z;
	c.z = a.z*b.w-a.y*b.x+a.x*b.y+a.w*b.z;
	return c;
}
inline vec3 rot(const quat&a, const vec3& b)
{
	quat c;
	c.w = -a.x*b(0)-a.y*b(1)-a.z*b(2);
	c.x = +a.w*b(0)-a.z*b(1)+a.y*b(2);
	c.y = +a.z*b(0)+a.w*b(1)-a.x*b(2);
	c.z = -a.y*b(0)+a.x*b(1)+a.w*b(2);
	vec3 d;
	d(0) = c.x*a.w-(+c.w*a.x-c.z*a.y+c.y*a.z);
	d(1) = c.y*a.w-(+c.z*a.x+c.w*a.y-c.x*a.z);
	d(2) = c.z*a.w-(-c.y*a.x+c.x*a.y+c.w*a.z);
	return d;
}
inline vec3 rotinv(const quat&a, const vec3& b)
{
	quat c;
	c.w = +a.x*b(0)+a.y*b(1)+a.z*b(2);
	c.x = +a.w*b(0)+a.z*b(1)-a.y*b(2);
	c.y = -a.z*b(0)+a.w*b(1)+a.x*b(2);
	c.z = +a.y*b(0)-a.x*b(1)+a.w*b(2);
	vec3 d;
	d(0) = c.x*a.w+(+c.w*a.x-c.z*a.y+c.y*a.z);
	d(1) = c.y*a.w+(+c.z*a.x+c.w*a.y-c.x*a.z);
	d(2) = c.z*a.w+(-c.y*a.x+c.x*a.y+c.w*a.z);
	return d;
}
/**********************************************/
inline quat pow(const quat& aa, const double& c)//assuming unit quaternion
{
	quat a =aa;
	a.normalize();
	quat b;
	double rx = sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
	double th = 2.*atan2(rx,a.w);
	if(rx==0 && a.w>0)
	{
		b.w=1;b.x=0;b.y=0;b.z=0;
		return b;
	}
	else 	if(rx==0 && a.w<=0)
	{
		b.w = cos(M_PI*c)   ;
		b.x = sin(M_PI*c)   ;
		b.y = 0.;
		b.z = 0.;
		return b;
	}
	else
	{
		double snrx = (rx==0)? 0. : (sin(th*c/2)/rx);
		b.w =     cos(th*c/2)   ;
		b.x = a.x*snrx;
		b.y = a.y*snrx;
		b.z = a.z*snrx;
		return b;
	}
}
/**********************************************/
inline vec3 qua2vec(const quat& a)
{
	vec3 b;
	double rx = sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
	if(rx<=0){b.zero();return b;} 
	double th = 2.*atan2(rx,a.w);
	if(th>=M_PI) th-=2.*M_PI;
	b(0) = (a.x/rx)*th;
	b(1) = (a.y/rx)*th;
	b(2) = (a.z/rx)*th;
	return b;
}
// inline double sinc(const double& x)
// {
// 	if(fabs(x)>1.0e-10) return sin(x)/x;
// 	else                return 1.0-x*x/6 + x*x*x*x/120;
// }
inline quat vec2qua(const vec3& a)
{
	double angh = 0.5* sqrt(a(0)*a(0)+a(1)*a(1)+a(2)*a(2));
	double sinc_angh = sinc(angh);
	quat b;
	b.w = cos(angh);
	b.x = (a(0)*0.5)*sinc_angh;
	b.y = (a(1)*0.5)*sinc_angh;
	b.z = (a(2)*0.5)*sinc_angh;
	return b;
}
/**********************************************/
inline void iniM_qua(mat3* const pmat, const quat& qua)
{
	const double& q0 = qua.w;
	const double& q1 = qua.x;
	const double& q2 = qua.y;
	const double& q3 = qua.z;
	mat3& mat = *pmat ;
	double r =qua.norm(); if(r>0) r=1./r;
	mat(0,0)=(q0*q0+q1*q1-q2*q2-q3*q3);
	mat(0,1)=(q1*q2-q0*q3)*2.;
	mat(0,2)=(q1*q3+q0*q2)*2.;
	mat(1,0)=(q1*q2+q0*q3)*2.;
	mat(1,1)=(q0*q0-q1*q1+q2*q2-q3*q3);
	mat(1,2)=(q2*q3-q0*q1)*2.;
	mat(2,0)=(q1*q3-q0*q2)*2.;
	mat(2,1)=(q2*q3+q0*q1)*2.;
	mat(2,2)=(q0*q0-q1*q1-q2*q2+q3*q3);
	mat *= r*r;
}

/**********************************************/
inline mat3 qua2mat(const quat& qua){mat3 out; iniM_qua(&out,qua);return out;}
inline quat mat2qua(const mat3& mat)
{
	quat q;
	double& m00 = mat(0,0);
	double& m01 = mat(0,1);
	double& m02 = mat(0,2);
	double& m10 = mat(1,0);
	double& m11 = mat(1,1);
	double& m12 = mat(1,2);
	double& m20 = mat(2,0);
	double& m21 = mat(2,1);
	double& m22 = mat(2,2);
	double tr = m00 + m11 + m22;
	if (tr > 0) { 
	  double S = sqrt(tr+1.0) * 2; // S=4*qw 
	  q.w = 0.25 * S;
	  q.x = (m21 - m12) / S;
	  q.y = (m02 - m20) / S; 
	  q.z = (m10 - m01) / S; 
	} else if ((m00 > m11)&(m00 > m22)) { 
	  double S = sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx 
	  q.w = (m21 - m12) / S;
	  q.x = 0.25 * S;
	  q.y = (m01 + m10) / S; 
	  q.z = (m02 + m20) / S; 
	} else if (m11 > m22) { 
	  double S = sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
	  q.w = (m02 - m20) / S;
	  q.x = (m01 + m10) / S; 
	  q.y = 0.25 * S;
	  q.z = (m12 + m21) / S; 
	} else { 
	  double S = sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
	  q.w = (m10 - m01) / S;
	  q.x = (m02 + m20) / S;
	  q.y = (m12 + m21) / S;
	  q.z = 0.25 * S;
	}
	double r = q.norm();
	q.w /= r;
	q.x /= r;
	q.y /= r;
	q.z /= r;
	return q;
}
/**********************************************/
inline vec3 operator -(const quat& a, const quat& b)
{
	/*   a = FqA , b = FqB , c = Fq{AB}      */
	return qua2vec(nearer(a*inv(b)));
}
inline vecvec operator -(const vecqua& a, const vecqua& b)
{
	vecvec c;
	c.sv0() = a.vec()-b.vec() ;
	c.sv1() = a.qua()-b.qua() ;
	return c;
}
/**********************************************/
inline quat operator +(const quat& a, const vec3& b)// 
{
	return normalize(nearer(vec2qua(b))*a);
}
inline quat qua_plus_vec3_nonearer(const quat& a, const vec3& b)
{
	return normalize(vec2qua(b)*a);
}
inline vecqua operator +(const vecqua& a, const vecvec& b)
{
	vecqua c;
	c.vec() = a.vec()+b.sv0() ;
	c.qua() = a.qua()+b.sv1() ;
	return c;
}
inline vecqua& vecqua::operator +=(const vecvec& b)
{
	this->vec() += b.sv0() ;
	this->qua() = this->qua()+b.sv1() ;
	return (*this);
}
/**********************************************/
inline vecqua tplus(const vecqua& a, const vecvec& b)
{
    vecqua c;
    c.vec() = a.vec() + qua2mat(a.qua())*b.sv0();
    c.qua() = normalize(a.qua()*nearer(vec2qua(b.sv1())));
    return c;
}
inline vecqua tplus(const vecqua& a, const vecvec& b, const vecvec& c){return tplus(tplus(a,b),c);}
inline vecqua tplus(const vecqua& a, const vecvec& b, const vecvec& c, const vecvec& d){return tplus(tplus(a,b,c),d);}
inline vecqua tplus(const vecqua& a, const vecvec& b, const vecvec& c, const vecvec& d, const vecvec& e){return tplus(tplus(a,b,c,d),e);}
/**********************************************/
inline vecvec   operator/(const vecvec& a, const double& b){vecvec c; for(int i=0;i<6;i++)c(i)=a(i)/b   ; return c;}
inline vecvec   operator*(const double& a, const vecvec& b){vecvec c; for(int i=0;i<6;i++)c(i)=a*b(i)   ; return c;}
inline vecvec   operator*(const vecvec& a, const double& b){vecvec c; for(int i=0;i<6;i++)c(i)=a(i)*b   ; return c;}
inline vecvec   operator+(const vecvec& a, const vecvec& b){vecvec c; for(int i=0;i<6;i++)c(i)=a(i)+b(i); return c;}
inline vecvec   operator-(const vecvec& a, const vecvec& b){vecvec c; for(int i=0;i<6;i++)c(i)=a(i)-b(i); return c;}
inline vecvec   operator-(const vecvec& a                 ){vecvec c; for(int i=0;i<6;i++)c(i)=-a(i)    ; return c;}
/**********************************************/

#endif
