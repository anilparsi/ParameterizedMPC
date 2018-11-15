#pragma once
#include <string>
#include "DefineSettings.h"
#include "ActiveConstraints.h"

/*! \class QPSolver
 * \brief Solve quadratic programming problems using active set approach.
 *
 * QP formulation:
 *	min		0.5 z^T G z + g^T z
 *	s.t		lbineq <= AiZ z <= ubineq
 * 
 * The matrix G is specified in terms of the inverse of its Cholesky decomposition (Li):
 * G = LL^T
 * Li =  inv(L)
 */
class QPSolver{

public:
	/*! 
	 * \brief Constructor to use MATLAB interface
	 * 
	 * Constructor to use the QPSolver object with MATLAB
	 * \param dir contains the address of the directory with the matrices Li, g, lbineq,
	 * AiZ, ubineq in .txt files.
	 */
	QPSolver(std::string dir);

	/*! 
	 * \brief Constructor for C++ interface
	 * 
	 * Constructor to use the QPSolver object with C++
	 * \param Li_i is the matrix containing Li
	 * \param g_i is the matrix containing g
	 * \param Aineq_i is the matrix containing AiZ
	 * \param lbineq_i is the matrix containing lbineq
	 * \param ubineq_i is the matrix containing ubineq
	 * \param nz_i is the number of variables nz
	 * \param nc_i is the number of constraints nc
	 * \param tolMin_i is the minimum tolerance used in the active set method. This is the default tolerance initially used by the active set solver.
	 * \param tolMax_i is the maximum tolerance used in the active set method. The tolerance is increased to this value after iterRelax_i iterations are performed.
	 * \param MAXITER_i is the maximum number of iterations in the active set method
	 * \param iterRelax_i is the number of iterations after which the tolerance is relaxed
	 */
	QPSolver(const real_t*const Li_i, const real_t*const g_i, const real_t*const Aineq_i,
		const real_t*const lbineq_i, const real_t*const ubineq_i, const int_t nz_i, const int_t nc_i,
		const real_t tolMin_i = 1e-9, const real_t tolMax_i = 1e-5, const int_t MAXITER_i = 50, const int_t iterRelax_i = 25);
	
	/// destructor
	~QPSolver(); 
	
	/// function to solve QP using incremental active set approach
	void	solve();

	/*! \brief number of active set iterations
	 */
	int_t	getIterNumber() const {return iter;}
	
	/*! \brief get exit flag
	 */
	int_t	getExitFlag() const { return exitFlag; }

	/*! \brief returns the solution to QP
	 */
	void	getSolutionCopy(real_t *z_out) const;
private:
	/// perform initialization 
	void	initialize();
	/// perform standard active set approach (when lambda>0)
	void	activeSetIterations(const int_t extra_idx=0);

	/// calculate Lagrange multipliers for the current active set
	void	calcLambda();
	
	/// add constraint to the active set
	void	addConstraint(const  int_t viol_idx);
		
	/// calculate the error for a particular constraint
	void	calculateError(const int_t idx,const real_t *const x, real_t *const err) const;

	// params
	int_t	MAXITER;			///< Max iterations for Active Set approach
			
	// dimensions
	int_t	iter;				///< number of active set iterations
	
	// variables for active set iterations
	real_t	*z_sol,				///< solution of previous iteration
			*delta,				///< difference between current and previous solutions
			*a_del;				///< AiZ(inactive)*delta

	int_t	*indices;			///< inactive indices which move in direction of constraint

	int_t	exitFlag;			///< 0 for solved
								///< -1 for infeasible IC (comes from kickout constraint)
								///< -2 for maxIter in Primal active set method
								///< -3 for maxIter in solver

protected:
	real_t	*Li;				///< inverse of Cholesky decomposition of G
	
	int_t	nc,					///< total number of inequality constraints;
			nz;					///< number of decision variables in the QP

	real_t	*z,					///< values of decision variables at each iteration
			*g,					///< linear part of cost function in QP
			*AiZ,				///< inequality constraints	
			*lbineq,			///< lower bound of inequality constraints			
			*ubineq,			///< upper bound of inequality constraints		

			*temp_nz,			///< temporary variable
			*temp_nz2,			///< temporary variable

			*LiTLi,				///< LiTLi = Li^T * Li	
			*lambda,			///< Lagrange multipliers of active set
		
			tolMin,				///< minimum tolerance used for checking constraints
			tolMax,				///< maximum tolerance used for checking constraints
			TOL;				///< absolute tolerance for active set approach 

	int_t	iterRelax,			///< iteration at which tolerance is relaxed to maxTol
			viol_idx;			///< index of maximum violation; negative index for lb
								// add 1 to absolute index because 0 and -0 are same
	
	/// class containing active constraint coefficients
	ActiveConstraints *activeCons;
								
	bool	viol;				///< violation status

	/// check constraints of the QP for violations
	virtual void	checkConstraints();

	/*! \brief calculates the solution for the current active set, where the constraints in the active set are 
	 * considered as equality constraints, and all other constriants are neglected.
	 */
	virtual void	calc_z();

};

