#include <RandomNumGen.h>
#include <Visualise.h>
#include <Storage.h>
/*
	***SIMULATED NUCLEAR DECAY OF AN ATOM***

	Uses an pseudo-random number generator to simulate the spontaneous emission of radiation
	from an atom due to nuclear decay. Exhibits the statistical average exponential decay for
	large numbers of nuclei but which breaks down into a more stochastic process as the number
	of remaining nuclei drop below a certain number.
*/


int main(){

	/*	Instantiate a (pseudo)random number generator with uniform distribution on the
		interval [0,1). The "seeding" is handled by the RNG constructor */
	phys::RNG<> rand_0_1(0.0, 1.0);

	int N = 1000;					//number of nuclei
	size_t M = 1000;				//total number of time steps
	size_t N0 = N;

	double lambda = 0.01;			//decay constant
	double dt = 1.;				//time step
	
	std::vector<double> time;		//container for the time
	std::vector<double> nuclei;		//container for the number of nuclei remaining

	/* store initial values */
	size_t j = 0;
	time.push_back(double(dt*j));
	nuclei.push_back(double(N));


	/*	For each time step loop through nuclei remaining to see if they decay or not. 
		Early exit if all nuclei decay before reaching last time step  */
	while(j++ < M && N > 0){
		for(size_t i = 0; i < N0; ++i){
			if(rand_0_1.random_number() < lambda*dt) --N; 
		}
		N0 = N;
		if(j%10 == 0 && N > 0){ //store every 10th point and avoid taking log of zero
			time.push_back(double(dt*j));
			nuclei.push_back(double(N));
		}
	}

	phys::Storage<double> nuclear_container("time", "N");
	nuclear_container.copy(time, nuclei);

	//nuclear_container.write("./nuclear_decay.txt", true);

	phys::Viewer viewer;
	viewer.plot(nuclear_container);

	return 0;

}
