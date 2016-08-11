#ifndef DJW_BVPODE_H
#define DJW_BVPODE_H

#include "DynVector.h"
#include "DynMatrix.h"
#include "LinearSolvers.h"
#include "FiniteDifferenceGrid.h"
#include "SecondOrderODE.h"
#include "BoundaryConditions.h"

#ifdef _HAVE_EIGEN
#include <Eigen/Sparse>
#endif

namespace phys{

#ifdef _HAVE_EIGEN
	typedef Eigen::SparseMatrix<double> sparseMat;
#endif

	class BvpODE{
	public:

		BvpODE(SecondOrderODE* pODE, BoundaryConditions* pBCs, uint numNodes);
		~BvpODE();

	private:
		BvpODE(const BvpODE&);
		const BvpODE& operator=(const BvpODE&);

	public:

		/* solves the boundary value ODE problem and writes the solution to the file if specified */
		stdVec_d solve(const std::string& filename = std::string());
		void plot() const; 

	private:
		void setupMatrix();
		void setupVector();
		void applyBCs();
		void appendBCs();
		void writeToFile(const std::string& fn);

	private:
		SecondOrderODE* m_pODE;
		BoundaryConditions* m_pBCs;
		uint m_numNodes;

		FiniteDifferenceGrid* m_pGrid;
				
		mat* m_pLhsMat; 
		stdVec_d* m_pRhsVec;
		stdVec_d m_SolVec;

#ifdef _HAVE_EIGEN
		sparseMat m_LHSsparse;
#endif

	};
}

#endif
