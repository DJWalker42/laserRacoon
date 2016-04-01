#ifndef MESH_HPP
#define MESH_HPP

#include <adv/physList.hpp>
#include <physVector.hpp>
#include <physMaths.hpp>
#include <adv/FiniteElement.hpp>

namespace phys{

	template <class T>
	class mesh : public connectedList<T>{
	public:
		mesh();
		mesh(T& e);

		// indexing the nodes in the mesh
		size_t indexing();

		/*The following functions are currently implemented for triangle element type only*/

		//used with mesh<triangle> only; refinement step
		void refine(const std::vector<double>&, double);

		//used with mesh<triangle> only; refine the neighbour of a refined triangle
		void refine_neighbour(node<point2>&, node<point2>&, node<point2>&);

		//used with mesh<triangle> only; the maximum modulus a nodes away from (>0.01) a singularity; 
		double maxNorm(const std::vector<double>&);

		//used with mesh<triangle> only; mesh refinement at a circular domain
		void refineBoundary(size_t);
	};
}

#include <adv/mesh.inl>

#endif