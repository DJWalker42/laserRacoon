#include "Storage.h"
#include "Visualise.h"


int main() {

	phys::ODEStorage ode_store("./euler_shm.log", 2, 1);

	phys::Viewer viewer;
	viewer.plot(ode_store, phys::Viewer::PHASE);


	return 0;
}
