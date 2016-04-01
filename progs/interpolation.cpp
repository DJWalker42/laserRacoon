#include <Interpolation.h>
#include <Visualise.h>
#include <iomanip>

using namespace phys::interp;

double f (double x)
{
	return (x == 0)? 1.0 : sin(x)/x;
}


int main()
{
	size_t dim = 10;

	phys::stdVec_d x;
	phys::stdVec_d y;

	for(size_t i = 0; i <= dim; ++i)
		x.push_back(2*i); 

	for(size_t i = 0; i < x.size(); ++i)
		y.push_back(f(x[i]));

	data my_data(x,y);

	std::vector<Interpolator*> interpolators(5);

	interpolators[0] = new Linear(my_data);
	interpolators[1] = new Lagrange(my_data);
	interpolators[2] = new Aitken(my_data);
	interpolators[3] = new UpDown(my_data); 
	interpolators[4] = new Cubic(my_data);

	std::vector<std::string> key_names(5); 

	key_names[0] = "Linear";
	key_names[1] = "Lagrange";
	key_names[2] = "Aitken";
	key_names[3] = "UpDown";
	key_names[4] = "Spline";

	phys::stdVec_d x_interp;

	for(size_t i = 0; i < 199; ++i)
		x_interp.push_back((i+1)*0.1);

	std::vector<phys::stdVec_d> y_interp(5); 

	for(size_t i = 0; i < y_interp.size(); ++i)
		y_interp[i] = interpolators[i]->interpolate(x_interp);

	phys::visual::Viewer viewer;

	viewer.set_shape(phys::visual::CROSS, 10);

	viewer.set_y_range(-0.3, 1.1);

	//plotter.animation(true);
	viewer.set_key_name( key_names ); 	
	viewer.plot(x, y);

	viewer.set_shape(phys::visual::CIRCLE, 1);
	

	viewer.add_data(x_interp, y_interp);

	return 0;
}
