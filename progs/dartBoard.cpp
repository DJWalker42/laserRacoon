#include <RandomNumGen.h>
#include <Storage.h>
#include <Visualise.h>
#include <Maths.h>
#include <iostream>


int main(){
	size_t repeat, darts_to_throw;

	std::cout << "Input the number of darts to throw: ";
	std::cin >> darts_to_throw;

	std::cout << "Input the number of repeats: "; 
	std::cin >> repeat;

	phys::stdVec_d bin(100,0.);
	phys::stdVec_d pi_est;
	
	for(size_t k = 0; k < repeat; ++k)
	{
		phys::RNG<> rand_0_1(0.,1.);//,k);
		std::vector<double> v_rand = rand_0_1.random_vector(2*darts_to_throw);

		std::vector<double> result;

		for(size_t i = 0; i < 2*darts_to_throw; i += 2)
			result.push_back( sqrt(v_rand[i] * v_rand[i] + v_rand[i+1] * v_rand[i+1]) );

		int hit = 0, miss = 0;
		std::vector<double>::iterator it = result.begin();
		for(;it != result.end(); ++it){
			if(*it > 1.0) 
				++miss;
			else
				++hit;
		}
		result.clear();

		pi_est.push_back(4.0*double(hit)/double((hit+miss)));

		size_t idx = size_t((pi_est[k]*100.0/0.3) - 1000);// bin range [3.0, 3.3) steps 0.003

		if(idx < 100) ++bin[idx]; //else out-of-range of bins
		//std::cout << "Done " << k << "/" << repeat << "\n"; 
	}

	std::pair<double, double> stats = phys::maths::mean_stdDev(pi_est);

	std::cout << "No darts thrown: " << darts_to_throw << "\n";
	std::cout << "Mean value: " << stats.first << " Sigma:" << stats.second << "\n";

	phys::stdVec_d xbin(100);

	for(size_t i = 0; i < 100; ++i)
		xbin[i] = 3.0 + 0.003 * i;

	phys::storage::Storage<double> pi_dart("Pi estimate", "Frequency");
	pi_dart.copy(xbin, bin);

	phys::visual::Viewer viewer;

	viewer.set_x_range(3., 3.3);

	viewer.plot(pi_dart);

	return 0;

}
