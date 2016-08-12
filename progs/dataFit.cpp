#include <DataFit.h>
#include <Storage.h>
#include <Visualise.h>

#include <string>

/*
	***SHOW CASE OF THE DATA FITTING ROUTINES***
	
	Here we use a sample of Millikan's oil drop experimental data values to show how the data fitting algorithms can be used.
*/

int main(int argc, char** argv)
{
	std::string filename = "../../../../laserRacoon/progs/resource/millikanData.txt";

	//... or use an argument in the call to the binary
	if (argc > 1) {
		filename = std::string (argv[1]); 
	}

	phys::storage::Storage<double> data;

	data.read(filename); 

	std::pair<double,double> llsq_result = phys::dfit::llsq(data.get_independent(), data.get_dependent()); 

	std::cout << "Gradient = " << llsq_result.first << ", y intercept = " << llsq_result.second  << "\n";

	std::pair<double,double> direct_res = phys::dfit::linear_fit(data.get_independent(), data.get_dependent());

	std::cout << "Gradient = " << direct_res.first << ", y intercept = " << direct_res.second  << "\n";

	std::pair<phys::stdVec_d,phys::mat> poly_fit = phys::dfit::polynomial_fit(size_t(6), data.get_independent(), data.get_dependent());

	std::cout << "coeff = \t" << poly_fit.first << "\n";

	phys::mat U = poly_fit.second;

	size_t n = U[0].size();

	phys::stdVec_d x = data.get_independent();

	double sum_0 = 0.0;
	double sum_t = 0.0;
	for(size_t i = 0; i < n; ++i){
		sum_0 += U[0][i] * U[0][i];
		sum_t += x[i]*U[0][i] * U[0][i];
	}

	double e0 = poly_fit.first[1];
	double de = poly_fit.first[0] - poly_fit.first[1]*sum_t/sum_0; 

	std::cout << " electron charge = " << e0 << "\t +/- " << de << "\n";

	phys::visual::Viewer viewer;

	viewer.set_shape(phys::visual::CROSS, 10);

	viewer.plot(data); 

	viewer.draw_line(llsq_result.first, llsq_result.second); 

	return 0;
}
