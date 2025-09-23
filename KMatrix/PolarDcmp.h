/**** Decompose.c ****/
/* Ken Shoemake, 1993 */
#include <math.h>


namespace KNSforPolarDecomp
{
	typedef  double HMatrix[4][4]; /* Right-handed, for column vectors */
	double polar_decomp(HMatrix M, HMatrix Q, HMatrix S);
	/******* Matrix Preliminaries *******/
	/** Copy nxn matrix A to C using "gets" for assignment **/
	#define mat_copy(C,gets,A,n) {for(int i=0;i<n;i++) for(int j=0;j<n;j++) C[i][j] gets (A[i][j]);}
	/** Assign nxn matrix C the element-wise combination of A and B using "op" **/
	#define mat_binop(C,gets,A,op,B,n) {int i,j; for(i=0;i<n;i++) for(j=0;j<n;j++) C[i][j] gets (A[i][j]) op (B[i][j]);}
	/** Multiply the upper left 3x3 parts of A and B to get AB **/
	inline void mat_mult(const HMatrix A, const HMatrix B, HMatrix AB)
	{
		for (int i=0; i<3; i++) for (int j=0; j<3; j++) AB[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j];
	}

	/** Return dot product of length 3 vectors va and vb **/
	inline double vdot(const double *va, const double *vb)
	{
			return (va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2]);
	}

	/** Set v to cross product of length 3 vectors va and vb **/
	inline void vcross(const double *va, const double *vb, double *v)
	{
			v[0] = va[1]*vb[2] - va[2]*vb[1];
			v[1] = va[2]*vb[0] - va[0]*vb[2];
			v[2] = va[0]*vb[1] - va[1]*vb[0];
	}

	/** Set MadjT to transpose of inverse of M times determinant of M **/
	inline void adjoint_transpose(const HMatrix M, HMatrix MadjT)
	{
			vcross(M[1], M[2], MadjT[0]);
			vcross(M[2], M[0], MadjT[1]);
			vcross(M[0], M[1], MadjT[2]);
	}


	/******* Decomp Auxiliaries *******/

	static HMatrix mat_id = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

	/** Compute either the 1 or infinity norm of M, depending on tpose **/
	inline double mat_norm(HMatrix M, int tpose)
	{
		double sum, max;
		max = 0.0;
		for (int i=0; i<3; i++)
		{
			if (tpose) sum = fabs(M[0][i])+fabs(M[1][i])+fabs(M[2][i]);
			else	     sum = fabs(M[i][0])+fabs(M[i][1])+fabs(M[i][2]);
			if (max<sum) max = sum;
		}
		return max;
	}

	inline double norm_inf(HMatrix M) {return mat_norm(M, 0);}
	inline double norm_one(HMatrix M) {return mat_norm(M, 1);}

	/** Return index of column of M containing maximum abs entry, or -1 if M=0 **/
	inline int find_max_col(HMatrix M)
	{
			double abs, max;
			int col;
			max = 0.0; col = -1;
			for (int i=0; i<3; i++) for (int j=0; j<3; j++)
			{
				abs = M[i][j]; if (abs<0.0) abs = -abs;
				if (abs>max) {max = abs; col = j;}
			}
			return col;
	}

	/** Setup u for Household reflection to zero all v components but first **/
	inline void make_reflector(const double *v, double *u)
	{
			double s = sqrt(vdot(v, v));
			u[0] = v[0]; u[1] = v[1];
			u[2] = v[2] + ((v[2]<0.0) ? -s : s);
			s = sqrt(2.0/vdot(u, u));
			u[0] = u[0]*s; u[1] = u[1]*s; u[2] = u[2]*s;
	}

	/** Apply Householder reflection represented by u to column vectors of M **/
	inline void reflect_cols(HMatrix M, double *u)
	{
			for (int i=0; i<3; i++) 
			{
				double s = u[0]*M[0][i] + u[1]*M[1][i] + u[2]*M[2][i];
				for (int j=0; j<3; j++) M[j][i] -= u[j]*s;
			}
	}
	/** Apply Householder reflection represented by u to row vectors of M **/
	inline void reflect_rows(HMatrix M, double *u)
	{
			int i, j;
			for (i=0; i<3; i++) {
		double s = vdot(u, M[i]);
		for (j=0; j<3; j++) M[i][j] -= u[j]*s;
			}
	}

	/** Find orthogonal factor Q of rank 1 (or less) M **/
	inline void do_rank1(HMatrix M, HMatrix Q)
	{
			double v1[3], v2[3], s;
			int col;
			mat_copy(Q,=,mat_id,4);
			/* If rank(M) is 1, we should find a non-zero column in M */
			col = find_max_col(M);
			if (col<0) return; /* Rank is 0 */
			v1[0] = M[0][col]; v1[1] = M[1][col]; v1[2] = M[2][col];
			make_reflector(v1, v1); reflect_cols(M, v1);
			v2[0] = M[2][0]; v2[1] = M[2][1]; v2[2] = M[2][2];
			make_reflector(v2, v2); reflect_rows(M, v2);
			s = M[2][2];
			if (s<0.0) Q[2][2] = -1.0;
			reflect_cols(Q, v1); reflect_rows(Q, v2);
	}

	/** Find orthogonal factor Q of rank 2 (or less) M using adjoint transpose **/
	inline void do_rank2(HMatrix M, HMatrix MadjT, HMatrix Q)
	{
			double v1[3], v2[3];
			double w, x, y, z, c, s, d;
			int col;
			/* If rank(M) is 2, we should find a non-zero column in MadjT */
			col = find_max_col(MadjT);
			if (col<0) {do_rank1(M, Q); return;} /* Rank<2 */
			v1[0] = MadjT[0][col]; v1[1] = MadjT[1][col]; v1[2] = MadjT[2][col];
			make_reflector(v1, v1); reflect_cols(M, v1);
			vcross(M[0], M[1], v2);
			make_reflector(v2, v2); reflect_rows(M, v2);
			w = M[0][0]; x = M[0][1]; y = M[1][0]; z = M[1][1];
			if (w*z>x*y) 
			{
				c = z+w; s = y-x; d = sqrt(c*c+s*s); c = c/d; s = s/d;
				Q[0][0] = Q[1][1] = c; Q[0][1] = -(Q[1][0] = s);
			}
			else
			{
				c = z-w; s = y+x; d = sqrt(c*c+s*s); c = c/d; s = s/d;
				Q[0][0] = -(Q[1][1] = c); Q[0][1] = Q[1][0] = s;
			}
			Q[0][2] = Q[2][0] = Q[1][2] = Q[2][1] = 0.0; Q[2][2] = 1.0;
			reflect_cols(Q, v1);
			reflect_rows(Q, v2);
	}


	/******* Polar Decomposition *******/

	/* Polar Decomposition of 3x3 matrix in 4x4,
	 * M = QS.  See Nicholas Higham and Robert S. Schreiber,
	 * Fast Polar Decomposition of An Arbitrary Matrix,
	 * Technical Report 88-942, October 1988,
	 * Department of Computer Science, Cornell University.
	 */
	const double TOL = 1.0e-6;

};
	// A = U H where U is unitary (rotational) and H is hermitian



inline double iniM_PolarDcmp(KMatrix<double,3,3>* const pU, const KMatrix<double,3,3>& A)
{
	KNSforPolarDecomp::HMatrix Mk, MadjTk, Ek;
	double det, M_one, M_inf, MadjT_one, MadjT_inf, mu, gamma, g1, g2;
	for(int i=0;i<3;i++) for(int j=0;j<3;j++) Mk[i][j] = A(j,i);
	M_one = KNSforPolarDecomp::norm_one(Mk);  
	M_inf = KNSforPolarDecomp::norm_inf(Mk);
	do
	{
		KNSforPolarDecomp::adjoint_transpose(Mk, MadjTk);
		det = KNSforPolarDecomp::vdot(Mk[0], MadjTk[0]);
		if (det==0.0) {KNSforPolarDecomp::do_rank2(Mk, MadjTk, Mk); break;}
		MadjT_one = KNSforPolarDecomp::norm_one(MadjTk);
		MadjT_inf = KNSforPolarDecomp::norm_inf(MadjTk);
		gamma = sqrt(sqrt((MadjT_one*MadjT_inf)/(M_one*M_inf))/fabs(det));
		g1 = gamma*0.5;
		g2 = 0.5/(gamma*det);
		mat_copy(Ek,=,Mk,3);
		mat_binop(Mk,=,g1*Mk,+,g2*MadjTk,3);
		mat_copy(Ek,-=,Mk,3);
		mu = KNSforPolarDecomp::norm_one(Ek);
		M_one = KNSforPolarDecomp::norm_one(Mk);
		M_inf = KNSforPolarDecomp::norm_inf(Mk);
	}
	while (mu > (M_one*KNSforPolarDecomp::TOL ));
	for(int i=0;i<3;i++)for(int j=0;j<3;j++) (*pU)(i,j) = (Mk[j][i]);
	return (det);
}

inline double iniMU_PolarDcmp(KMatrix<double,3,3>* const pU, KUSymMat<double,3>* const pS, const KMatrix<double,3,3>& A)
{
	double det = iniM_PolarDcmp(pU,A);
	iniM_mulMTMO(&(pS->asMatrix()),*pU,A);
	for(int i=0;i<3;i++)for(int j=i+1;j<3;j++) (*pS)(i,j)=0.5*((*pS)(i,j)+(*pS)(j,i));
	return (det);
}
