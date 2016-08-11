#include <Fourier.h>
#include <ODESolvers.h>
#include <Visualise.h>

/*
	***PROGRAM SOLVING THE VAN DER POL OSCILLATOR EQUATION***

	Uses the Runge-Kutta algorithm to compute the solution to the Van Der Pol oscillator.
	Passes the stable part of that solution to a fast fourier transform to have a look at the
	frequencies making up the solution.
*/

int main()
{
	const double M = 10;
	const size_t N = size_t( pow(2,M) ); //1024

	double h = 120.0/double(3 * N); //0.0390625

	phys::ode::state initial_sys(0.0, 1.0, 0.0);

	phys::diffs::Van_Der_Pol* VDP = new phys::diffs::Van_Der_Pol();

	phys::ode::RK4 runge_kutta(VDP, initial_sys, h);

	phys::storage::ODE_Storage solution = runge_kutta.fullSolve(double(3*N) * h);
    solution.set_x_name("Time");
    solution.set_y_name("Amplitude");

	phys::visual::Viewer viewerA;
	viewerA.set_x_range(0., 20.);
	viewerA.withLines();
	viewerA.plot(solution);

	phys::stdVec_d time = solution.get_independent();
	phys::stdVec_d posn = solution.get_dependent();

	phys::stdVec_d sub_pos = phys::sub_vector(posn, N, posn.size() - 1);
	phys::stdVec_d sub_tim = phys::sub_vector(time, N, time.size() - 1);

	phys::stdVec_c F(N);

	for(size_t i = 0; i < N; ++i)
		F[i] = phys::complex (sub_pos[i], 0.0);

	phys::stdVec_c G = phys::fourier::FFT(F);

	phys::stdVec_d normG(N);

	//The fourier coefficients of the VDP solution have a large dynamic range
	//To compare the coefficient peaks we need to take the log of the values to plot.
	for(size_t i = 0; i < N; ++i){
		G[i] /= double(N); 
		normG[i] = log10(phys::norm(G[i]));
	}

	//FFT solution contains both positive and negative frequencies
	//take only the postivie frequencies (first half of FFT solution).
	normG.resize(512);
	phys::stdVec_d x_axis(512);

	//create x-axis in frequency units
	for(size_t i = 0; i < 512; ++i)
		x_axis[i] = 2 * 3.14159 * double(i)/N/h;

	phys::visual::Viewer viewerB;

	viewerB.set_x_range(0., 10.);
    viewerB.set_x_name("Frequency");
    viewerB.set_y_name("Power");
	viewerB.withLines();

	viewerB.plot(x_axis, normG); 

	phys::storage::Storage<double> transform;

	transform.copy(x_axis, normG);
	//transform.write(/*supply a location*/);

	return 0;
}
