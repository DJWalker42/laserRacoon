#include "Storage.h"
#include "Visualise.h"


int main() {

	//file created and written by odeSolvers program
	phys::ODEStorage ode_store("./euler_shm.log", 2, 1);

	phys::Viewer viewer;
	viewer.plot(ode_store, phys::Viewer::PHASE);


	return 0;
}
