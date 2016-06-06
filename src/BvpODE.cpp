#include <iostream>
#include <fstream>
#include <cassert>

#include <BvpODE.h>
#include <Visualise.h>

namespace phys{
	
	/*	Note that we are embedding the boundary conditions into the rhs vector of knowns, and
		the lhs matrix of coefficients such that the dimensions of the problem are 2 fewer than
		the number of nodes used in the 1D grid */

	BvpODE::BvpODE(	SecondOrderODE* pODE, 
					BoundaryConditions* pBCs, 
					size_t numNodes) :
					m_pODE(pODE),
					m_pBCs(pBCs),
					m_numNodes(numNodes),
					m_pGrid(new FiniteDifferenceGrid(m_numNodes, m_pODE->m_xMin, m_pODE->m_xMax)),
					m_pLhsMat(new mat(numNodes - 2, numNodes - 2)),
					m_pRhsVec(new stdVec_d(numNodes - 2)),					
					m_SolVec(numNodes - 2)					
#ifdef _HAVE_EIGEN
					, m_LHSsparse(numNodes - 2, numNodes - 2)
#endif
	{}

	BvpODE::~BvpODE()
	{
		delete m_pLhsMat;
		delete m_pRhsVec;
		delete m_pGrid;
	}

	stdVec_d BvpODE::solve(const std::string& fn)
	{
		setupMatrix(); //sets up m_pLhsMat or the sparseMat with Eigen
		setupVector(); //sets up m_pRhsVec
		applyBCs();
#ifdef _HAVE_EIGEN
		{
			using namespace Eigen;
			SparseLU<sparseMat> solver(m_LHSsparse);
			Map<VectorXd> b(m_pRhsVec->data(), m_pRhsVec->size()); //copy m_pRhsVec (std::vector<double>) to b 
			VectorXd x = solver.solve(b);
			VectorXd::Map(m_SolVec.data(), x.size()) = x; //copy x (Eigen::VectorXd) to m_SolVec (std::vector<double>)
		}
#else
		LinearSolver* solver = new LU(*m_pLhsMat);
		m_SolVec = solver->solve(*m_pRhsVec);
		delete solver;
#endif
		appendBCs();
		if (fn.empty() == false){
			writeToFile(fn);	
		}
		return m_SolVec; 
	}

	void BvpODE::plot() const
	{
		visual::Viewer viewer;
		viewer.plot(m_pGrid->getGrid1D(), m_SolVec);
	}


	void BvpODE::setupMatrix()
	{
#ifdef _HAVE_EIGEN
		/*	Eigen::Triplet in this context contains the row index, the column index, 
			and the value of the element at that location, respectively 
		*/

		typedef Eigen::Triplet<double> triplet;

		std::vector<triplet> triplet_list;
		triplet_list.reserve(m_numNodes); 

		//special cases for boundary conditions

		//LHS boundary
		double xm = m_pGrid->m_Nodes[0].getX();
		double x = m_pGrid->m_Nodes[1].getX();
		double xp = m_pGrid->m_Nodes[2].getX();

		double alpha = 2. / (xp - xm) / (x - xm);
		double beta = -2. / (xp - x) / (x - xm);
		double gamma = 2. / (xp - xm) / (xp - x);

		triplet_list.push_back(triplet(0, 0, m_pODE->m_Uxx * beta + m_pODE->m_U));
		triplet_list.push_back(triplet(0, 1, m_pODE->m_Uxx * gamma + m_pODE->m_Ux / (xp - xm)));

		//RHS boundary - might be more efficient to put this after the interior loop to maintain 
		//order in the triplet list before calling the member function .setFromTriplets()
		xm = m_pGrid->m_Nodes[m_numNodes - 3].getX();
		x = m_pGrid->m_Nodes[m_numNodes - 2].getX();
		xp = m_pGrid->m_Nodes[m_numNodes - 1].getX();

		alpha = 2. / (xp - xm) / (x - xm);
		beta = -2. / (xp - x) / (x - xm);
		gamma = 2. / (xp - xm) / (xp - x);

		triplet_list.push_back(triplet(m_numNodes - 3, m_numNodes - 4, m_pODE->m_Uxx * alpha - m_pODE->m_Ux / (xp - xm)));
		triplet_list.push_back(triplet(m_numNodes - 3, m_numNodes - 3, m_pODE->m_Uxx * beta + m_pODE->m_U));

		//interior
		for (size_t i = 1; i < m_numNodes - 3; ++i)
		{
			xm = m_pGrid->m_Nodes[i - 1].getX();
			x = m_pGrid->m_Nodes[i].getX();
			xp = m_pGrid->m_Nodes[i + 1].getX();

			alpha = 2. / (xp - xm) / (x - xm);
			beta = -2. / (xp - x) / (x - xm);
			gamma = 2. / (xp - xm) / (xp - x);

			triplet_list.push_back(triplet(i, i - 1, m_pODE->m_Uxx * alpha - m_pODE->m_Ux / (xp - xm)));
			triplet_list.push_back(triplet(i, i, m_pODE->m_Uxx * beta + m_pODE->m_U));
			triplet_list.push_back(triplet(i, i + 1, m_pODE->m_Uxx * gamma + m_pODE->m_Ux / (xp - xm)));
		}

		m_LHSsparse.setFromTriplets(triplet_list.begin(), triplet_list.end());
#else
		//special cases for boundary conditions

		//LHS boundary
		double xm = m_pGrid->m_Nodes[0].getX();
		double x = m_pGrid->m_Nodes[1].getX();
		double xp = m_pGrid->m_Nodes[2].getX();

		double alpha = 2. / (xp - xm) / (x - xm);
		double beta = -2. / (xp - x) / (x - xm);
		double gamma = 2. / (xp - xm) / (xp - x);

		//note that alpha is not used in the following two lines
		(*m_pLhsMat)[0][0] = m_pODE->m_Uxx * beta + m_pODE->m_U;
		(*m_pLhsMat)[0][1] = m_pODE->m_Uxx * gamma + m_pODE->m_Ux / (xp - xm); 

		//RHS boundary - same notation as before with i = m_numNodes - 2
		xm = m_pGrid->m_Nodes[m_numNodes - 3].getX();
		x = m_pGrid->m_Nodes[m_numNodes - 2].getX();
		xp = m_pGrid->m_Nodes[m_numNodes - 1].getX();

		alpha = 2. / (xp - xm) / (x - xm);
		beta = -2. / (xp - x) / (x - xm);
		gamma = 2. / (xp - xm) / (xp - x); 

		//note that gamma is not used in the following two lines
		(*m_pLhsMat)[m_numNodes - 3][m_numNodes - 4] = m_pODE->m_Uxx * alpha - m_pODE->m_Ux / (xp - xm);
		(*m_pLhsMat)[m_numNodes - 3][m_numNodes - 3] = m_pODE->m_Uxx * beta + m_pODE->m_U;

		//Interior nodes

		for (size_t i = 1; i < m_numNodes - 3; ++i)
		{
			xm = m_pGrid->m_Nodes[i - 1].getX();
			x = m_pGrid->m_Nodes[i].getX();
			xp = m_pGrid->m_Nodes[i + 1].getX();

			alpha = 2. / (xp - xm) / (x - xm);
			beta = -2. / (xp - x) / (x - xm);
			gamma = 2. / (xp - xm) / (xp - x);

			(*m_pLhsMat)[i][i - 1]	= m_pODE->m_Uxx * alpha - m_pODE->m_Ux / (xp - xm);
			(*m_pLhsMat)[i][i]		= m_pODE->m_Uxx * beta + m_pODE->m_U;
			(*m_pLhsMat)[i][i + 1]	= m_pODE->m_Uxx * gamma + m_pODE->m_Ux / (xp - xm);
		}

#endif

	}

	void BvpODE::setupVector()
	{
		for (size_t i = 1; i < m_numNodes - 1; ++i) 
		{
			double x = m_pGrid->m_Nodes[i].getX();
			(*m_pRhsVec)[i - 1] = m_pODE->m_pDiff(x);
		}
	}

	void BvpODE::applyBCs()
	{
		/*
			As we specify the boundary types and values in the sole constructor for the
			BoundaryConditions class then the data memembers of that class are gaurenteed 
			to be initialised. 
		*/

		//LHS boundary - xm == x[i-1] xp == x[i+1] here i = 1
		double xm = m_pGrid->m_Nodes[0].getX();
		double x = m_pGrid->m_Nodes[1].getX();
		double xp = m_pGrid->m_Nodes[2].getX();

		double alpha = 2. / (xp - xm) / (x - xm);
		double gamma = 2. / (xp - xm) / (xp - x);

		if (m_pBCs->m_LHStype == DIRICHLET)
		{
			(*m_pRhsVec)[0] -= ((m_pODE->m_Uxx * alpha - m_pODE->m_Ux / (xp - xm)) * m_pBCs->m_LHSvalue); 
		}

		if (m_pBCs->m_LHStype == NEUMANN)
		{	

#ifdef _HAVE_EIGEN
			m_LHSsparse.coeffRef(0, 0) += (m_pODE->m_Uxx * alpha - m_pODE->m_Ux / (xp - xm));
#else
			(*m_pLhsMat)[0][0] += (m_pODE->m_Uxx * alpha - m_pODE->m_Ux / (xp - xm));
#endif
			(*m_pRhsVec)[0] += ((m_pODE->m_Uxx * alpha - m_pODE->m_Ux / (xp - xm)) * m_pBCs->m_LHSvalue * (x - xm));
		}

		//RHS boundary
		xm = m_pGrid->m_Nodes[m_numNodes - 3].getX();
		x = m_pGrid->m_Nodes[m_numNodes - 2].getX();
		xp = m_pGrid->m_Nodes[m_numNodes - 1].getX();

		alpha = 2. / (xp - xm) / (x - xm);
		gamma = 2. / (xp - xm) / (xp - x);


		if (m_pBCs->m_RHStype == DIRICHLET)
		{
			(*m_pRhsVec)[m_numNodes - 3] -= ((m_pODE->m_Uxx * gamma + m_pODE->m_Ux / (xp - xm)) * m_pBCs->m_RHSvalue);
		}

		if (m_pBCs->m_RHStype == NEUMANN)
		{
#ifdef _HAVE_EIGEN
			m_LHSsparse.coeffRef(m_numNodes - 3, m_numNodes - 3) += (m_pODE->m_Uxx * gamma + m_pODE->m_Ux / (xp - xm));
#else
			(*m_pLhsMat)[m_numNodes - 3][m_numNodes - 3] += (m_pODE->m_Uxx * gamma + m_pODE->m_Ux / (xp - xm));
#endif
			(*m_pRhsVec)[m_numNodes - 3] -= ((m_pODE->m_Uxx * gamma + m_pODE->m_Ux / (xp - xm)) * m_pBCs->m_RHSvalue * (xp - x));
		}
	}

	void BvpODE::appendBCs()
	{
		double lhs_boundary = 0., rhs_boundary = 0.;

		if (m_pBCs->m_LHStype == DIRICHLET)
		{
			lhs_boundary = m_pBCs->m_LHSvalue;
		}

		if (m_pBCs->m_RHStype == DIRICHLET)
		{
			rhs_boundary = m_pBCs->m_RHSvalue;
		}

		if (m_pBCs->m_LHStype == NEUMANN)
		{
			double xm = m_pGrid->m_Nodes[0].getX();
			double x = m_pGrid->m_Nodes[1].getX(); 

			//extrapolate using the boundary gradient
			lhs_boundary = m_SolVec[0] - m_pBCs->m_LHSvalue * (x - xm); //eqn for a straight line
		}

		if (m_pBCs->m_RHStype == NEUMANN)
		{
			double x = m_pGrid->m_Nodes[m_numNodes - 2].getX();
			double xp = m_pGrid->m_Nodes[m_numNodes - 1].getX();

			//extrapolate using the boundary gradient
			rhs_boundary = m_SolVec[m_numNodes - 3] +  m_pBCs->m_RHSvalue * (xp - x);
		}

		m_SolVec.emplace(m_SolVec.begin(), lhs_boundary);
		m_SolVec.push_back(rhs_boundary);
	}



	void BvpODE::writeToFile(const std::string& fn)
	{
		std::ofstream output_file(fn);
		if (output_file.is_open())
		{
			for (size_t i = 0; i < m_numNodes; ++i)
			{
				output_file << m_pGrid->m_Nodes[i + 1].getX() << "\t" << m_SolVec[i] << "\n";
			}

			output_file.flush();
			output_file.close();
			std::cout << "Solution written to " << fn << "\n";
		}
		else //fix by asking the user to supply a new path and filename then calling function recursively; q to quit
		{
			std::cout << "Unable to open " << fn << " for writing\n";
			std::cout << "Please supply a new location and filename or enter q to quit (no data will be saved):\n";
			std::string alt_fn;
			std::cin >> alt_fn;

			if (alt_fn == "q")
			{
				std::cout << "Quitting -- no output written\n";
			}
			else
			{
				writeToFile(alt_fn); //using new filename
			}
		}

		return;
	}

}
