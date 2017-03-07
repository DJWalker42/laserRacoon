#include <RandomNumGen.h>
#include <Storage.h>
#include <Visualise.h>
#include <Maths.h>
#include <iostream>

/*
	********** THE MONTE-CARLO DARTBOARD ***********

	Monte-Carlo dart board experiment to estimate a value for pi.
	User is asked for the number of darts to throw to get an estimate of pi.
	Then asked the number of times to repeat the experiment for statistical analysis.
	
	Pi is estimated by throwing darts at a unit square board defined by (0,0), (1,0), (0,1), and (1,1). 
	All darts thrown hit this unit square by design.
	By comparing the ratio of darts to land within a unit quarter circle (centred at the origin) to those that
	land outside said quarter circle we can obtain an estimate for pi.To explain, the probability that a randomly thrown
	dart will land within a given area is proportional to the size of that given area. Thus, the probability that a
	dart lands within the unit quarter circle is given by k.pi/4, where k is the constant of proportionality. 
	In this case we can state that k is unity as we have designed all darts thrown to hit the unit square. Thus, the ratio
	of "hits" (landed in the quarter circle) to "misses" (not landed in the quarter circle) should equal pi/4.
	
	Experiment with the number of darts thrown vs. the number of repeats paying specific attention to the statistics.
	
*/


int main(){
	size_t repeat, darts_to_throw;

	std::cout << "Input the number of darts to throw: ";
	std::cin >> darts_to_throw;

	std::cout << "Input the number of repeats: "; 
	std::cin >> repeat;

	phys::stdVec_d bin(100,0.); //vector to in-order to bin results to plot as a histogram
	phys::stdVec_d pi_est; //value to store the current pi estimate for a given experiment
	
	
	/** main loop **/
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
