#include <cassert>
#include <iostream>
#include <stdexcept>

#include "FiniteDifferenceGrid.h"

namespace phys{

	FiniteDifferenceGrid::FiniteDifferenceGrid(	uint numNodes, 
												double xMin, 
												double xMax) :
												m_Nodes(numNodes)
	{
		if(numNodes < 3)
		{
			throw std::runtime_error("A finite difference grid requires at least three nodes\n");
		}

		if( xMin > xMax)
		{
			std::cout << "Limits wrong way round -- swapping\n";
			double tmp = xMin;
			xMin = xMax;
			xMax = tmp;
		}

		double step = (xMax - xMin) / double (numNodes - 1);
		for(uint i = 0; i < m_Nodes.size(); ++i)
		{
			m_Nodes[i] = Node(xMin + i * step);
		}
	}


	FiniteDifferenceGrid::FiniteDifferenceGrid(const std::vector<Node>& customGrid) : m_Nodes(customGrid)
	{}

	std::vector<double> FiniteDifferenceGrid::getGrid1D() const
	{
		uint n = static_cast<uint>(m_Nodes.size());
		std::vector<double> retval(n);
		for (uint i = 0; i < n; ++i){
			retval[i] = m_Nodes[i].getX();
		}		
		return retval;
	}


	FiniteDifferenceGridZ::FiniteDifferenceGridZ(	uint numNodesX,
													double xMin,
													double xMax ) :
		m_boundaryNodes(),
		m_interiorNodes()
	{
		if (numNodesX < 3)
		{
			throw std::runtime_error("1D Finite difference grid requires at least three nodes\n");
		}

		if (xMin > xMax)
		{
			std::cout << "Limits wrong way round -- swapping\n";
			double tmp = xMin;
			xMin = xMax;
			xMax = tmp;
		}

		double step_x = (xMax - xMin) / double(numNodesX - 1); //uniform grid spacing in x
		
		m_boundaryNodes.push_back(Node(xMin));
		m_boundaryNodes.push_back(Node(xMax)); 

		m_interiorNodes.reserve(numNodesX - 2);

		for(uint i = 1; i < numNodesX - 1; ++i) //note the range
		{
			m_interiorNodes.push_back(Node(xMin + i * step_x));
		}

		m_dimensions = D1;
	}

	FiniteDifferenceGridZ::FiniteDifferenceGridZ(	uint numNodesX,
													double xMin,
													double xMax,
													uint numNodesY,
													double yMin,
													double yMax) :
													m_boundaryNodes(),
													m_interiorNodes()
	{
		if (numNodesX < 3 || numNodesY < 3)
		{
			throw std::runtime_error("2D Finite difference grid requires at least three nodes per dimension\n");
		}

		if (xMin > xMax)
		{
			std::cout << "X limits wrong way round -- swapping\n";
			double tmp = xMin;
			xMin = xMax;
			xMax = tmp;
		}

		if (yMin > yMax)
		{
			std::cout << "Y limits wrong way round -- swapping\n";
			double tmp = yMin;
			yMin = yMax;
			yMax = tmp;
		}

		double step_x = (xMax - xMin) / double(numNodesX - 1); //uniform grid spacing in x
		double step_y = (yMax - yMin) / double(numNodesY - 1); //uniform grid spacing in y

		m_boundaryNodes.reserve(2 * (numNodesX + numNodesY) - 4);
		m_interiorNodes.reserve((numNodesX - 2) * (numNodesY - 2));

		/********Boundary Nodes********/

		//first row of nodes
		for (int i = 0; i < numNodesX; ++i)
		{
			m_boundaryNodes.push_back(Node(xMin + i * step_x, yMin, 0)); //x coord, y coord, globalIdx (row index)
		}

		//boundaries across the interior per row (two either side of the domain)
		for (int j = 1; j < numNodesY - 1; ++j)
		{
			m_boundaryNodes.push_back(Node(xMin, yMin + j * step_y, j)); //left boundary
			m_boundaryNodes.push_back(Node(xMax, yMin + j * step_y, j)); //right boundary
		}

		//last row of nodes
		for (int i = 0; i < numNodesX; ++i)
		{
			m_boundaryNodes.push_back(Node(xMin + i * step_x, yMax, numNodesY - 1)); //x coord, y coord, globalIdx (row index)
		}

		/****Interior Nodes*****/

		for (uint j = 1; j < numNodesY - 1; ++j)
		{
			for (uint i = 1; i < numNodesX - 1; ++i)
			{
				m_interiorNodes.push_back(Node(xMin + i * step_x, yMin + j * step_y, j));
			}
		}

		m_dimensions = D2;
	}


	FiniteDifferenceGridZ::FiniteDifferenceGridZ(	uint numNodesX,
													double xMin,
													double xMax,
													uint numNodesY,
													double yMin,
													double yMax,
													uint numNodesZ,
													double zMin,
													double zMax) :
													m_boundaryNodes(),
													m_interiorNodes()
	{	
		/* Body of constructor left as an exercise */
	}


}

