#ifndef Headder_KMatrix_LUInverse
#define Headder_KMatrix_LUInverse

template<int N> inline
void iniM_invMO(KMatrix<double,N,N>* const pOut, const KMatrix<double,N,N>& A)		// LU����
{
	double piv;
	double max;
	double big, big0;
	double em;
	double sum;
	double scale[N];

	int i,j,k;
	int piv_i;				// �s�{�b�g�̔ԍ�
	int ip, kp, ib;
	int ips[N];				// �s�{�b�g�œ���ւ�����̕���

	KMatrix<double,N,N> Qmat = A;
	KMatrix<double,N,N>& Qinv = *pOut;
// ---- LU�ɕ��� ---------------------------
	for(i=0;i<N;i++){						// �s�{�b�g�I��
		ips[i] = i;
		max = 0.0;
		for(j=0;j<N;j++){
			if( max<fabs( Qmat(i,j) ) )	max = fabs( Qmat(i,j) );
		}
		scale[i] = 1.0/max;
	}

	for(k=0;k<N-1;k++){
		big = 0.0;
		for(i=k;i<N;i++){
			ip = ips[i];
			big0 = fabs( Qmat(ip,k) )*scale[ip];
			if( big0>big ){
				big = big0;
				piv_i = i;
			}
		}
		if( piv_i!=k ){
			j = ips[k];						// �s�{�b�g�̓���ւ�
			ips[k] = ips[piv_i];
			ips[piv_i] = j;
		}
		kp = ips[k];
		piv = Qmat(kp,k);
		for(i=k+1;i<N;i++){
			ip = ips[i];
			em = -Qmat(ip,k)/piv;
			Qmat(ip,k) = -em;					// L�̗v�f
			for(j=k+1;j<N;j++){
				Qmat(ip,j) += em*Qmat(kp,j);	// U�̗v�f
			}
		}
	}
	kp = ips[N-1];

// ----- �t�s��̓��o -------------------
	KMatrix<double,N,N> Id; 
	Qinv.zero();
	Id.zero();
	for(i=0;i<N;i++)	Id(i,i) = 1.0;

	for(k=0;k<N;k++){
		ip = ips[0];
		Qinv(0,k) = Id(ip,k);

		for(i=1;i<N;i++){
			ip = ips[i];
			sum = 0.0;
			for(j=0;j<i;j++){
				sum += Qmat(ip,j)*Qinv(j,k);
			}
			Qinv(i,k) = Id(ip,k) - sum;
		}

		ip = ips[N-1];
		Qinv(N-1,k) = Qinv(N-1,k)/Qmat(ip,N-1);

		for(ib=2;ib<=N;ib++){
			i = N - ib;
			ip = ips[i];
			sum = 0.0;
			for(j=i+1;j<N;j++){
				sum += Qmat(ip,j)*Qinv(j,k);
			}
			Qinv(i,k) = ( Qinv(i,k) - sum )/Qmat(ip,i);
		}
	}
}

#endif


