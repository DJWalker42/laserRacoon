#ifndef KRYLOV_ACCN_HPP
#define KRYLOC_ACCN_HPP

#include <adv/SparseMat.h>
#include <adv/Multigrid.h>

namespace phys{

	/*
		Preconditioned Conjugate Gradient method
		Conjugate Gradient Squared method
		Tranpose-free Quasi-minimal Residual
		General Minimal Residual
	*/

	template<class T>
	void PCG(	const Multigrid<T>& MG,
				const sparseMat<T>& A,
				const vector<T>& f,
				vector<T> x	)
	{
		const double eps = 1.e-15, thresh = 1.e-12;
		const int iter_num = 1000;
		vector<T> zero(x.size(), 0.); 
		vector<T> keep(x);
		vector<T> rr(MG.Vcycle(f,keep) - x);
		vector<T> r = f - A * x;
		vector<T> pp(rr);
		vector<T> p = r;
		double gamma = r * rr;
		double gamma_init = gamma;
		int count = 0; 
		while((abs(gamma/gamma_init) >= thresh * thresh) && (count++ <= iter_num))
		{
			keep = pp;
			vector<T> ww = pp - MG.Vcycle(zero,keep);
			vector<T> w = A * pp;
			double alpha = gamma / (pp * w);
			x += alpha * pp;
			rr -= alpha * ww;
			r -= alpha * w;
			double gamma_old = gamma;
			gamma = r * rr;
			double beta = gamma / gamma_old;
			pp = rr + beta * pp;
			std::cout << "At MG it. " << count << " in PCG, (r,Pr) = " << gamma << "\n";
		}
		std::cout << "Total MG it. in PCG = " << count << "\n";
	}// The PCG acceleration method

	template<class T>
	void CGS(	const Multigrid<T>& MG,
				const sparseMat<T>& A,
				int TFQMR,
				const vector<T>& f,
				const vector<T>& x	)
	{
		const double eps = 1.e-15, thresh = 1.e-12;
		const int iter_num = 1000;
		T omega_init, omega_0, omega_1, tau, vv;
		T rho_0, sigma, eta, alpha, rho_1, beta;
		vector<T> keep(x); 
		vector<T> rr(MG.Vcycle(f,keep) - x);
		vector<T> rbar(rr);
		vector<T> rcgs(rr), u(rr), pp(rr);
		tau = omega_init = omega_1 = omega_0 = sqrt(rcgs * rcgs);
		std::cout << "res_0 = " << omega_init << "\n";
		eta = 0.;
		vv = 0.;
		if(fabs(rho_0 = rbar * rr) < eps)
			std::cout << "rho_0 = " << fabs(rho_0) << "\n";
		vector<T> zero(x.size(), 0.);
		vector<T> d(zero), v(zero), q(zero);
		int count = 1;
		do
		{
			keep = pp;
			v = pp - MG.Vcycle(zero,keep);
			if(fabs(sigma = rbar * v) < eps)
				std::cout << "sigma = " << fabs(sigma) << " mg it. = " << 2 * count << "\n";
			if(fabs(alpha = rho_0 / sigma) < eps)
				std::cout << "alpha = " << fabs(alpha) << " mg it. = " << 2 * count << "\n";
			q = u - alpha * v;
			vector<T> uq = u + q;
			keep = uq;
			rcgs -= alpha * (uq - MG.Vcycle(zero,keep)); 
			omega_1 = sqrt(rcgs * rcgs);
			if(!TFMQR)
			{
				x += alpha * uq;
				std::cout << "res = " << fabs(omega_1) << " mg it. = " << 2 * count + 1 << "\n";
			}
			else
			{
				int m = 0;
				while( m++ < 2 ) //loop twice with m == 1 then m == 2 (see description of postfix ++)
				{
					T omega;
					if(m == 1)
					{
						omega = sqrt( omega_0 * omega_1 );
						keep = u;
					}
					else
					{
						omega = omega_1;
						keep = q;
					}
					T scalar = vv * vv * eta / alpha;
					d = keep + scalar * d;
					vv = omega / tau;
					T c = 1. / sqrt(1. + vv * vv);
					tau *= vv * c;
					eta = c * c * alpha;
					x += eta * d;
					std::cout << "res = " << sqrt( (A*x - f) * (A*x - f) ) 
						<< "it. = " << 2 * count + 1 << "\n";
				}
			}
			omega_0 = omega_1;
			if(fabs(rho_1 = rbar * rcgs) < eps)
				std::cout << "rho_1 = " << fabs(rho_1) << " it. = " 
				<< 2 * count + 1 << "\n";
			beta = rho_1 / rho_0;
			rho_0 = rho_1;
			u = rcgs + beta * q;
			pp = u + beta * (q + beta * pp);
		}while((fabs(omega_1 / omega_init) >= thresh) && (++count <= iter_num)); 
		std::cout << "Total MG it. in CGS = " << 2 * count + 1 << "\n";
	} // the CGS or the TFQMR acceleration method


	template<class T>
	const vector<T>& GMRES(	const Multigrid<T>&MG, 
							const sparseMat<T>&A,
							int preIterations, 
							int K,
							const vector<T>&f,
							vector<T>&x)
	{
		for(int i = 0; i < preIterations; ++i)
			MG.Vcycle(f,x);
		vector<T> s = x;
		vector<T> r = x;
		T R[K+1][K+1];
		for(int i=0;i<=K;i++)
			for(int j=0;j<=K;j++)
				R[i][j]=0.0;
		T Givens[2][K];
		T xi[K];
		for(int i = 0; i < K; i++)
		{
			Givens[0][i] = 0.0;
			Givens[1][i] = 0.0;
			xi[i] = 0.0;
		}
		vector<T>* Q[K+1];
		vector<T> zero(x.dim(),0.0);
		vector<T> keep(x);
		double res = sqrt((A*x-f) * (A*x-f));
		for(int k = 0;k <= K; k++)
		{
			if(k) keep = *Q[k-1];
			Q[k] = k ? new vector<T>(*Q[k-1]-MG.Vcycle(zero,keep)) : new vector<T>(MG.Vcycle(f,keep)-x);
			for(int j = 0; j < k; j++)
				*Q[k] -= (R[j][k] = *Q[j] * *Q[k]) * *Q[j];
			*Q[k] *= 1.0/(R[k][k] = sqrt(*Q[k] * *Q[k]));
			T givensa;
			T givensb;
			if(k)
			{
				for(int j = 1; j < k; j++)
				{
					givensa = R[j-1][k];
					givensb = R[j][k];
					R[j-1][k] = givensa*Givens[0][j-1] + givensb*Givens[1][j-1];
					R[j][k] = -givensa*Givens[1][j-1] + givensb*Givens[0][j-1];
				}
				T ab = sqrt( R[k-1][k]*R[k-1][k] + R[k][k]*R[k][k] );
				Givens[0][k-1] = R[k-1][k] / ab;
				Givens[1][k-1] = R[k][k] / ab;
				givensa = R[k-1][k];
				givensb = R[k][k];
				R[k-1][k] = givensa * Givens[0][k-1] + givensb*Givens[1][k-1];
				R[k][k] = 0.0;
				R[k][0] = -R[k-1][0] * Givens[1][k-1];
				R[k-1][0] = R[k-1][0] * Givens[0][k-1];
			}
			for(int i = k-1; i >= 0; i--)
			{
				xi[i] = R[i][0];
				for(int j = i+2; j <= k; j++)
					xi[i] -= R[i][j] * xi[j-1];
				xi[i] /= R[i][i+1];
			}
			s = x;
			for(int j = 0; j < k; j++)
				s += xi[j] * *Q[j];
			std::cout << "res at step k = " << k << " is " << res = sqrt((r = A*s - f)*r) << "\n";
		}
		return x = s;
	} // GMRES (with initial iterations) method

}

#endif