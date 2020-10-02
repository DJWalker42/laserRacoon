#include <ODESolvers.h>
#include <RootSearch.h>
#include <Storage.h>
#include <Visualise.h>

#include <algorithm>

/*
	***COMPUTES THE ELECTRON ENERGY STATES OF A GIVEN POTENTIAL FUNCTION USING THE NUMEROV ALGORITHM** 
	
	For a user defined potential computes the electron energy states of the well. It guesses at an
	energy level and integrates the resulting wavefunction from the classically forbidden region
	up to a matching point within the well. Then it integrates the wavefunction from the other
	direction again starting from the classically forbidden region, up to the matching point. It 
	use the information of the guessed energy level and how close the wavefunctions get at the matching
	point to improve the energy level guess; in essence it root searches for the precise energy level.

	The matching point can be any arbitrary point in the well however it convenient to take the
	matching point where the guessed energy level crosses the RHS boundary of the potential function.
	We integrate from the classically forbidden region into the well to suppress the unwanted exponential
	growth term in the solution to Schrodinger's equation that tends to blow up due to floating point
	number precision.

*/


double E;				//!< Energy eigenvalues to find in eV
double omega = 1.0;		//!< oscillator frequency (Hz)
/* 
	The size of the frequency is arbitrary but will depend on choice of energy and length units
	The natural choice for these is eV and nm making omega ~order 10^14 Hz, given the harmonic potential.
*/
double N = 500;						//!< Total number of steps to take in the integration
double x_lft = -3.0;				//!< Left extreme
double x_rht =  3.0;				//!< Right extreme
double step = (x_rht - x_lft)/N;	//!< step size
double u0 = 0.0;					//!< value of wavefunction at first grid point
double u1 = 1.e-10;					//!< value of wavefunction at second grid point

phys::stdVec_d phi;					//!< Storage for the wavefunction(s)
size_t i_match_p1;					//!< Index location where phi switches from left to right integration.

phys::ode::state u_lft(x_lft, u0, u1);	//!< Constructs for the Numerov solver Left initial
phys::ode::state u_rht(x_rht, u0, u1);	//!< Right initial

phys::diffs::Diff_eqn* q_func = new phys::diffs::User_eqn(); //auto deleted at end of program

//!< harmonic oscillator potential / square potential well (uncomment choice)
double V(double x)      
{
	using namespace phys::constants;
	return (0.5 * electron_mass * omega * omega * x * x);// harmonic potential
	//return fabs(x) > 1.0 ? 10.0 : 0.0; //square potential Vo = 10 eV, width = 2 nm
}

double find_match( double start, double step, double end )
{
	if(step < 0)
		step = (start < end) ? -step : step;
	else
		step = (start < end) ? step : -step;
	double match = start;         
    while (V(match) > E)  {  // in forbidden region
        match += step;
        if (fabs(match) > fabs(end)) {
            std::cerr << "can't find the turning point" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
	return match;
}

//!< q(x) function for the Numerov method.
double q (double x)
{
	using namespace phys::constants;
	return 2 * 10.0 * electron_mass / (hbar_ev * hbar_ev) * (E - V(x));
}

//interface to differential equation object -- compiler optimisations *may* remove redundant arguments here - but they are all cheap and quick to pass.
double phys::diffs::User_eqn::differential_function(double x, const stdVec_d& y, int N, int i)
{
	return q(x);
}

//!< Search function - input an energy level guess, returns with the difference between wavefunctions at matching point.
double F(double En)
{

	E = En; //assign search value to global E required by q(x) and find_match()

	//variables to ensure function is continuous over energy.
	//local static variables are preserved from one function call to the next.
	static int sign = 1;
	static int nodes = 0;

	// find right-hand crossing point in the potential function
	double match = find_match(x_rht, -step, x_lft); 

	//Storage for the wavefunction so we can conveniently apply the matching condition.
	//These will both contain:
	//[match - step, match, match + step] in that order.
	phys::stdVec_d phi_lft(3, 0.0), phi_rht(3, 0.0); 

	//integrate from the left up to one step past the matching point.(-5 inc.-> match + step)
	phys::ode::Numerov numerov_lft(q_func, u_lft, step);
	phys::storage::ODEStorage temp = numerov_lft.fullSolve(match + step);
	phys::stdVec_d lft_vals = temp.get_dependent();

	//lft_vals size - 1 gives switching index i_match + 1
	i_match_p1 = lft_vals.size() - 1;

	phi_lft[0] = lft_vals[i_match_p1 - 2];	//match - step
	phi_lft[1] = lft_vals[i_match_p1 - 1];	//match
	phi_lft[2] = lft_vals[i_match_p1];		//match + step

	//integrate from the right up to one step before the matching point (+5 dec.-> match - step)
	phys::ode::Numerov numerov_rht(q_func, u_rht, -step);
	temp = numerov_rht.fullSolve(match - step);
	phys::stdVec_d rht_vals = temp.get_dependent();

	size_t i_temp = rht_vals.size() - 1;

	phi_rht[0] = rht_vals[i_temp];		//match - step
	phi_rht[1] = rht_vals[i_temp - 1];	//match
	phi_rht[2] = rht_vals[i_temp - 2];	//match + step

	//scale the matching condition left values, using the matching point value
	double scale = phi_rht[1]/phi_lft[1];
	//scale all left values in phi and monitor nodes.

	{
		using phys::operator*=; 
		phi_lft *= scale;	
		lft_vals *= scale;
	}

	int n = 0;
	for(size_t i = 1; i <= i_match_p1; ++i){
		if(lft_vals[i -1] * lft_vals[i] < 0) ++n;
	}

	// flip the sign when a new node has developed
	if (n != nodes) {
		nodes = n;
		sign = -sign;
	}

	//merge data into phi (global wavefunction)
	lft_vals.insert(lft_vals.end(), rht_vals.begin(), rht_vals.end());
	phi = lft_vals;

	return sign * ( phi_rht[0] - phi_rht[2]- phi_lft[0] + phi_lft[2] )
		/ (2 * step * phi_rht[1]);
}


void normalize() 
{
	//We can erase the final three elements of phi as these 
	//contain the overlapped values from the right integration (match + 1, match, match - 1).
	phi.erase(phi.end() - 3, phi.end() );

	// manipulate phi so that it is in x order i.e. [x_lft -> x_rht]
	// current order is [x_lft -> match + 1, x_rht -> match + 2]
	// so we reverse the back "half" of the vector at the matching-point-plus-one plus one
    std::reverse(phi.begin() + i_match_p1 + 1, phi.end());

	double norm = phys::squaredNorm(phi);
	norm /= phi.size();
	norm = sqrt(norm);
	{
		using phys::operator/=; //this will do an elementwise division
		phi /= norm;
		//alternate syntax without the using
		/*phi = phys::operator/=(phi, norm);*/ 
	}
}

int main()
{
    //change this to where you want it saved
    //- make sure it exists first, the write/save functions don't make directories.
    std::string rootDir = "/path/to/data/directory/quantum_well";
    
	std::string pot_filename = rootDir + "/potential.txt";
	std::string lvl_filename = rootDir + "/levels.txt";
	std::string wvf_filename = rootDir + "/wavefunction.txt";

	std::cout	<< " Eigenvalues of the Schroedinger equation\n"
				<< " for the harmonic oscillator V(x) = 0.5 x^2\n"
				<< " ------------------------------------------\n";
	double E_max;

	do{
		std::cout << "Enter maximum energy E (must be less than " << V(x_rht) << ")\n";
		std::cin >> E_max;
	}while(E_max > V(x_rht));

	phys::storage::Storage<double> potential, levels, wavefunction;

	for(size_t i = 0; i <= size_t(N) ; ++i)
	{
		double x = x_lft + i*step;
		potential.store( x, V(x) );
		wavefunction.store(x);
	}

	potential.write(pot_filename);

    // find the energy levels
    size_t level = 0;           // energy level number
    double E_old = 0;           // previous energy eigenvalue attempted
    E = 0.01;                    // guess an initial E below the ground state	

	do
	{
        // estimate next E and dE
        double dE = 0.5 * (E - E_old);
		//monitor dE
        E_old = E;
        E += dE;
		phys::roots::Hybrid_B_S secant(F); //default tolerance 1.e-8
		if( secant.find_brackets(E, dE, E_max) )
		{
			E = secant.find_root();
			//store x coordinates and energy values of the levels (for drawing)
			secant.set_function(q);
			secant.find_brackets(x_lft, step, (x_rht + x_lft)/2); //use mid-point for end		
			levels.store(secant.lft_bracket(), E); 
			//normalise the wavefunction
			normalize();
			//store current wavefunction
			wavefunction.copy(phi);
			std::ostringstream convert;
			convert << level++;
			wavefunction.add_to_key( "E" + convert.str());
			convert.str("");

			std::cout << "Found a level at E = " << E << "\n";
			//std::cout << "constant = " << E/(2*level - 1) << "\n";//this should be a constant for the harmonic potential
		}
	}while(E < E_max);

	std::cout << "Found " << level << " levels\n";
	
	if(level == 0) return 0;
	
	levels.write(lvl_filename);
	wavefunction.write(wvf_filename);

	//compute square of the wavefunction (probability distribution) and
	//add on the particular E associated with the wavefunction for visuals
	std::vector<phys::stdVec_d> data = wavefunction.get_multi();
	phys::stdVec_d energy = levels.get_dependent();

	for(size_t j = 0; j < data.size(); ++j)
		for(size_t i = 0; i < data[j].size(); ++i){
			data[j][i] *= data[j][i];
			data[j][i] += energy[j];
		}

	phys::visual::Viewer viewer;
	viewer.set_key_name( wavefunction.get_key_names() );
	viewer.set_x_name("x/nm");
	viewer.set_y_name("Energy/eV");
	viewer.set_shape(phys::visual::CROSS, 5);
	viewer.set_x_range(-3,3);
	viewer.animation();
	viewer.plot(wavefunction.get_independent(), data);
	viewer.set_shape(phys::visual::SQUARE, 5);
	viewer.noAnimation();
	viewer.withLines();
	viewer.add_data(potential.get_independent(), potential.get_dependent(), true);
	viewer.save(rootDir + "/harmonic.png");

	return 0;
}
