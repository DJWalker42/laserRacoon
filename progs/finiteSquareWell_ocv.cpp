#include "RootSearch.h"
#include "PhysicalUnits.h"
#include "Storage.h"

#include "DynVector.h"
#include "Visualise.h"

#include <cmath>
#include <sstream>

/*
 *  Can this be rewritten to remove the use of global variables?
 *  One sticking point is that the RootSearch classes expect a function
 *  prototype that excepts a single, double argument. One possible
 *  solution is to make the matching functions static member functions
 *  of a class that has V0, half_width, and hbar2 as its data members.
 */


//convenience globals to compute hbar2
double h = phys::constants::planck;
double eV = phys::constants::electron_charge;
double me = phys::constants::electron_mass;
double pi = phys::constants::PI;

// using these units for hbar^2 keeps numbers computer friendly
// and we can say the mass of the electron is unity in calculations
double hbar2 = h * h * 1.e2 / 4 / pi / pi/ eV / me;

/*
 * Note that if these change the root brackets are no longer valid.
 * Is there a way of programming a "root look-up" automatically?
 */

double V0 = 10.; //eV
double half_width = 5.; //Angstroms


double even(double energy) {

	double alpha = sqrt(2. * energy/hbar2);
	double beta = sqrt(2 * (V0 - energy)/ hbar2);

	return beta * cos(alpha * half_width) - alpha * sin(alpha * half_width);
}

double odd(double energy) {
	double alpha = sqrt(2. * energy/hbar2);
	double beta = sqrt(2 * (V0 - energy)/ hbar2);

	return alpha * cos(alpha * half_width) + beta * sin(alpha * half_width);
}



int main (int argc, char ** argv) {

	using Bisecant = phys::roots::Hybrid_B_S;
	using Storage = phys::storage::Storage<double>;

	//you may want to change the location of the data file
	std::string wavefunction_filename {"./wavefunctions.log"};

	//properties of the linespace to plot data
	const double xmin {-10.}, xmax {10.}, xstep {0.1};

	//the following brackets are only valid for V0 = 10, half_width = 5

	//left and right brackets of the even parity functions
	phys::stdVec_d even_left {0.1, 2.5, 7.0};
	phys::stdVec_d even_right {0.5, 3.0, 7.5};

	//left and right brackets of the odd parity functions
	phys::stdVec_d odd_left {1.0, 4.5, 9.5};
	phys::stdVec_d odd_right {1.5, 5.0, 9.9};


	Bisecant even_search(&even);
	Bisecant odd_search(&odd);

	phys::stdVec_d energies;

	//find each energy level (3 even, 3 odd)
	for (int i = 0; i < 3; ++i) {
		even_search.set_brackets(even_left[i], even_right[i]);
		odd_search.set_brackets(odd_left[i], odd_right[i]);

		energies.push_back(even_search.find_root());
		energies.push_back(odd_search.find_root());
	}

	std::cout << "hbar2: " << hbar2 << std::endl;
	std::cout << "energy roots: " << energies << std::endl;


	//compute the wavefunction over a line-space x for regions I, II and III
	//and store the potential energy function for plotting
	phys::stdVec_d x_linespace;
	Storage wavefunctions; //for visualisation
	double x_val {xmin};
	while (x_val <= xmax) {

		if (x_val < -half_width) {
			wavefunctions.store(x_val, V0); //region I
		} else if (x_val < half_width) {
			wavefunctions.store(x_val, 0.); //region II
		} else {
			wavefunctions.store(x_val, V0); //region III
		}

		x_linespace.push_back(x_val);
		x_val += xstep;
	}

	std::cout << "x min: " << *(x_linespace.begin()) << " x max: " << *(x_linespace.end() - 1) << std::endl;

	//add a header for the potential function
	wavefunctions.add_to_key("V(x)");

	// We compute the wavefunctions for these energies, normalised such that
	// integral [-inf, inf] psi^2 dx = 1
	for (int i = 0; i < 6; ++i) {
		double alpha {sqrt(2 * energies[i] / hbar2)};
		double beta {sqrt(2 * (V0 - energies[i]) / hbar2)};
		double D {sqrt(beta / (1 + beta * half_width))};
		double C;
		if ((i % 2) == 0) {
			//even parity
			C = D * cos(alpha * half_width) * exp(beta * half_width);
		} else {
			//odd parity
			C = -D * sin(alpha * half_width) * exp(beta * half_width);
		}


		phys::stdVec_d psi;

		for (auto& x : x_linespace) {
			if (x < -half_width) {
				//region I
				psi.push_back(C * exp(beta * x));

			} else if ( x < half_width) {
				//region II
				if ((i % 2) == 0) {
					//even parity
					psi.push_back(D * cos(alpha * x));
				} else {
					//odd parity
					psi.push_back(D * sin(alpha * x));
				}
			} else {
				//region III - for odd parity solutions reverse sense of constant C
				if ((i % 2) == 0) {
					//even parity
					psi.push_back(C * exp(-beta * x));
				} else {
					//odd parity
					psi.push_back(-C * exp(-beta * x));
				}
			} //region condition end
		} //x_linespace loop end

		{
			using phys::operator+=;
			psi += energies[i]; //element-wise addition of the corresponding energy for visuals
		}

		//wavefunctions and potential share the same x axis
		wavefunctions.copy(psi); //for plotting each wavefunction
		std::stringstream ss;
		ss << i;
		wavefunctions.add_to_key("E" + ss.str()); //add a unique key label
	}


	wavefunctions.write(wavefunction_filename, true); //2nd arg. == print headers

	phys::visual::Viewer viewer;
	viewer.set_key_name(wavefunctions.get_key_names());
	viewer.set_x_name("x/Angs.");
	viewer.set_y_name("energy/eV");
	viewer.plot(wavefunctions.get_independent(), wavefunctions.get_multi());

	//press space-bar to progress through the data

	return 0;
}
