#ifndef DJW_NODE_H
#define DJW_NODE_H

namespace phys{

	/*
		Node class: describes a node or point in a 1D or 2D space for use in the FiniteDifferenceGrid class
	*/
	class Node{
	public:

		/* Constructor for a Node in 1D space */
		Node(double xCoord = 0.) : m_Xcoord(xCoord), m_Ycoord(0.), m_globalIndex(0)
		{}
		/* Constructor for a Node in 2D space */
		Node(	double xCoord, 
				double yCoord, 
				int globalIdx) : 
				m_Xcoord(xCoord), 
				m_Ycoord(yCoord),
				m_globalIndex(globalIdx)
		{}

		double getX() const {return m_Xcoord;}
		double getY() const {return m_Ycoord;}

	private:
		double m_Xcoord;
		double m_Ycoord;
		int m_globalIndex;
	};

	/*
		The code below is in a developmental stage and the comments that follow are my thoughts on how to design code 
		for a 2D (or higher dimensional) general finite difference grid. 

		First assume we cannot add or remove nodes from the grid once it is created - requires reordering/refactoring 
		of indicies. However, we can move individual nodes within the grid, up to a limit of the nearest neighbour 
		location (indicies would have to swap otherwise). The Boundary Value ODE solver class has been written to 
		solve on 1D grids of general spacing i.e. not necessarily grids of uniformly spaced nodes - note that a 2D
		implementation of the BvpODE solver class remains to be written. 

		We use the Node class in the FiniteDifferenceGrid class to construct the grid; the grid is simply a container of 
		nodes with a description of the connections between them. The nodes are either interior nodes or boundary nodes 
		(any corner nodes **must** be boundary nodes by defintion) and they must not be confused (an argument perhaps for
		avoiding inheritence/polymorphism?) 

		Note that if the grid is a uniform array (or its indicies can be described as a uniform array) then the
		FiniteDifferenceGrid class should be made responsible for computing the connections between nodes; it's a trival 
		task. For an arbritrary grid the responsibility of the connections between nodes lies with the user.  

		As written, the BoundaryNode class currently shares the same internal description as a Node class. However, 
		it will require data members to descripe which boundary it is on and/or to which other nodes it is connected. 

		As written, the InteriorNode class is a Node so should inherit from the Node class rather than repeat code. 
		An InteriorNode stores the locations of its four nearest neighbours for reference; we're assuming a 2D grid.

		To how many other nodes are BoundaryNodes connected? Three; two boundary nodes, one interior node

		What about nodes at the corners of the grid, assuming it has corners? Two boundary nodes, different boundaries

		How would you handle circular or arbritrary grid shapes? - not an easy answer to this question.
  
		Would you consider a BoundaryNode to be a special case of an InteriorNode, given the internal 
		description of an InteriorNode? Would it better to keep theses classes distinct? A Node is either an InteriorNode
		or a BoundaryNode, it cannot be both. 

		As an alternative to seperate classes we could give the Node class an attribute flag (set up as an internal enum
		say) that identifies it as either a interior node or a boundary node. 
	*/

	enum {NORTH, EAST, SOUTH, WEST};

	class BoundaryNode{
	public:
		BoundaryNode(double xCoord, double yCoord = 0., int globalIdx = 0) :
			m_Xcoord(xCoord), m_Ycoord(yCoord), m_globalIdx(globalIdx)
		{}
	private:
		double m_Xcoord;
		double m_Ycoord;
		int m_globalIdx;

	};


	class InteriorNode{
	public:
		InteriorNode(	double xCoord, 
						double yCoord = 0., 
						int globalIdx = 0, 
						int northIdx = 0, 
						int eastIdx = 0, 
						int southIdx = 0, 
						int westIdx = 0	) :
						m_Xcoord(xCoord),
						m_Ycoord(yCoord),
						m_globalIdx(globalIdx)
		{
			m_neighbour[NORTH] = northIdx;
			m_neighbour[EAST] = eastIdx;
			m_neighbour[SOUTH] = southIdx;
			m_neighbour[WEST] = westIdx;
		}
	private:
		double m_Xcoord;
		double m_Ycoord;
		int m_globalIdx;
		int m_neighbour[4];
	};

}


#endif
