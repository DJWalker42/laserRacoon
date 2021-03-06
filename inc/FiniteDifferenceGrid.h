#ifndef DJW_FINITEDIFFERENCEGRID_H
#define DJW_FINITEDIFFERENCEGRID_H

#include <vector>
#include "Node.h"

#include <cstdint>

namespace phys {

//forward declaration of the BvpODE class so we can declare it as a friend
class BvpODE;

class FiniteDifferenceGrid {
	friend class BvpODE;
public:
	/* Use this constructor for a uniform 1D grid */
	FiniteDifferenceGrid(uint32_t numNodes, double xMin, double xMax);
	/* Use this constructor for a non-uniform 1D grid, you have to specify the nodes */
	FiniteDifferenceGrid(const std::vector<Node> &customGrid);

	std::vector<double> getGrid1D() const;
private:
	std::vector<Node> m_Nodes;
};

/* Suggested class for a more general description of a finite grid that should be able to deal with one, two, and three dimensions */
class FiniteDifferenceGridZ {
	friend class BvpODE;
public:
	/* Constructor for a uniform 1D grid */
	FiniteDifferenceGridZ(uint32_t numNodesX, double xMin, double xMax);

	/* Constructor for a uniform 2D grid */
	FiniteDifferenceGridZ(uint32_t numNodesX, double xMin, double xMax,
			uint32_t numNodesY, double yMin, double yMax);

	/* Constructor for a uniform 3D grid */
	FiniteDifferenceGridZ(uint32_t numNodesX, double xMin, double xMax,
			uint32_t numNodesY, double yMin, double yMax, uint32_t numNodesZ,
			double zMin, double zMax);

	/* How might you write a constructor for a user defined grid that can be used for one, two, or three dimensions? */

private:
	std::vector<Node> m_boundaryNodes;
	std::vector<Node> m_interiorNodes;
	enum tDimensions {
		D1, D2, D3
	} m_dimensions;
};

}

#endif
