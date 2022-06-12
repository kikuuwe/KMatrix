#ifndef KQUATERNION_H_INCLUDED
#define KQUATERNION_H_INCLUDED


class KQuaternion
{ 
public:
  double q0,q1,q2,q3;
  KVector<double,3>& im() const{ return *((KVector<double,3>*)(&q1)); }
	void init(){q0=1;q1=q2=q3=0;}
	void normalize()
	{
		double s = sqrt(q0*q0+q1*q1+q2*q2+q3*q3) ;
		if(s>0) {s=1./s; q0 *= s;  q1 *= s;  q2 *= s;  q3 *= s;}
		else    {q0=1;q1=q2=q3=0;}
	};
	KQuaternion adj()
	{
		KQuaternion q;
		q.q0   = q0 ;
		q.im() = -im();
		return q;
	}

	KQuaternion operator*(const double& a)
	{
		KQuaternion q;
		q.q0 = q0*a;
		q.q1 = q1*a;
		q.q2 = q2*a;
		q.q3 = q3*a;
		return q;
	}
	KQuaternion operator+(const KQuaternion& a)
	{
		KQuaternion q;
		q.q0 = q0 + a.q0;
		q.q1 = q1 + a.q1;
		q.q2 = q2 + a.q2;
		q.q3 = q3 + a.q3;
		return q;
	}
	KQuaternion& operator+=(const KQuaternion& a)
	{
		q0 += a.q0;
		q1 += a.q1;
		q2 += a.q2;
		q3 += a.q3;
		return *this;
	}
	void get_rotvec(KVector<double,3>* const pout) const
	{
		double s = sqrt(q1*q1+q2*q2+q3*q3);
		if(s<=0){pout->zero(); return;}
		*pout = im();
		*pout *= asin(   max(-1,min(1, 2.*s*q0 )) )/s;
	}
	double get_angle() const
	{
		double s = sqrt(q1*q1+q2*q2+q3*q3);
		return asin(   max(-1,min(1, 2.*s*q0 )) );
	}
	void put_rotvec(const KVector<double,3>& in) 
	{
		double s = vabs(in); // radian, 0 to PI
		q0 = cos(s*0.5);
		if(s>0){im() = in; im() *= sin(s*0.5)/s; }
		else    im().zero();
	}
};


/******************************************************************/

inline
void init_mulQ_QOQO(KQuaternion* const p_ans, const KQuaternion& q, const KQuaternion&  p) 
{ 
	p_ans->q0    = q.q0*p.q0 ;
	p_ans->q0 	-= inner(q.im(),p.im());
	iniV_outer(& p_ans->im() , p.im()  ,  q.im());
	addX_mulXOS(&(p_ans->im()), p.im()  ,  q.q0  ) ;
	addX_mulXOS(&(p_ans->im()), q.im()  ,  p.q0  ) ;
}
inline
void init_mulQ_QOQA(KQuaternion* const p_ans, const KQuaternion& q, const KQuaternion&  p) 
{ 
	p_ans->q0    = q.q0*p.q0 ;
	p_ans->q0 	+= inner(q.im(),p.im());
	iniV_outer(& p_ans->im() ,               q.im()  ,  p.im());
	           addX_mulXOS(&  p_ans->im() , p.im()  , -q.q0  ) ;
	           addX_mulXOS(&  p_ans->im() , q.im()  ,  p.q0  ) ;
} 
inline
void init_mulQ_QAQO(KQuaternion* const p_ans, const KQuaternion& q, const KQuaternion&  p) 
{ 
	p_ans->q0    = q.q0*p.q0 ;
	p_ans->q0 	+= inner(q.im(),p.im());
	iniV_outer (& p_ans->im() ,            q.im()  ,  p.im());
	addX_mulXOS(& p_ans->im() , p.im()  ,  q.q0  ) ;
	addX_mulXOS(& p_ans->im() , q.im()  , -p.q0  ) ;
} 

/******************************************************************/
inline
void init_mulQ_QOVO(KQuaternion* const p_ans, const KQuaternion& q, const KVector<double,3>& p_im) 
{ 
	p_ans->q0    =  - inner(q.im(),p_im);

	iniV_outer(& p_ans->im() ,               p_im   ,  q.im());
	 addX_mulXOS(&            p_ans->im() , p_im  ,  q.q0  ) ;
}
inline
void init_mulQ_QAVO(KQuaternion* const p_ans, const KQuaternion& q, const KVector<double,3>&  p_im) 
{ 
	p_ans->q0    =     inner(q.im(),p_im);
	iniV_outer(& p_ans->im() ,               q.im() ,  p_im  );
	     addX_mulXOS(&        p_ans->im() , p_im   ,  q.q0  ) ;

} 
/******************************************************************/
inline
void init_mulV_QOQO(KVector<double,3>* const p_ans, const KQuaternion& q, const KQuaternion& p) 
{ 
	iniV_outer( p_ans       ,               p.im()  ,  q.im());
	  addX_mulXOS(          p_ans,      p.im()  ,  q.q0  ) ;
	 addX_mulXOS(           p_ans,       q.im()  ,  p.q0  ) ;
} 
inline
void init_mulV_QOQA(KVector<double,3>* const p_ans, const KQuaternion& q, const KQuaternion& p) 
{ 
	iniV_outer( p_ans       ,               q.im()  ,  p.im());
	 addX_mulXOS (            p_ans,        p.im()  , -q.q0  ) ;
	  addX_mulXOS (           p_ans,     q.im()  ,  p.q0  ) ;
} 
/******************************************************************/



/*********************/
inline
void init_mulQ_QAQOQA(KQuaternion* const p_out, const KQuaternion& a, const KQuaternion& b, const KQuaternion& c) 
{
	KQuaternion ab;
	init_mulQ_QAQO( &ab , a, b);
	init_mulQ_QOQA( p_out,ab, c);
}
/*********************/
inline
void init_quadV_QAVOQO(KVector<double,3>* p_out, const KVector<double,3>& a, const KQuaternion& q) 
{
	KQuaternion qa;
	init_mulQ_QAVO( &qa , q, a);
	init_mulV_QOQO( p_out,qa, q);
}
/*********************/
inline
void init_quadV_QOVOQA(KVector<double,3>* p_out, const KVector<double,3>& a, const KQuaternion& q) 
{
	KQuaternion qa;
	init_mulQ_QOVO( &qa , q, a);
	init_mulV_QOQA( p_out,qa, q);
}
/*********************/
inline
void init_rotate_vec(KVector<double,3>* p_out, const KVector<double,3>& a, const KQuaternion& q) 
{
	init_quadV_QAVOQO( p_out, a, q);
}

/*********************/
inline
void init_qrate_from_omg(KQuaternion* const p_dq, const KVector<double,3>& omg, const KQuaternion& q)
{
	p_dq->q0  = - q.q1 * omg(0) ; 
	p_dq->q0 += - q.q2 * omg(1) ; 
	p_dq->q0 += - q.q3 * omg(2) ; 
                     
	p_dq->q1  =   q.q0 * omg(0) ; 
	p_dq->q1 +=   q.q3 * omg(1) ; 
	p_dq->q1 += - q.q2 * omg(2) ; 
                     
	p_dq->q2  = - q.q3 * omg(0) ; 
	p_dq->q2 +=   q.q0 * omg(1) ; 
	p_dq->q2 +=   q.q1 * omg(2) ; 
                     
	p_dq->q3  =   q.q2 * omg(0) ; 
	p_dq->q3 += - q.q1 * omg(1) ; 
	p_dq->q3 +=   q.q0 * omg(2) ; 

	p_dq->q0 *= 0.5;
	p_dq->q1 *= 0.5;
	p_dq->q2 *= 0.5;
	p_dq->q3 *= 0.5;
}




#endif
