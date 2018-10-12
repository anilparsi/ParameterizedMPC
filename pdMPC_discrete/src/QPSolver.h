#pragma once
#include <string>
#include "DefineSettings.h"
#include "ActiveConstraints.h"

class QPSolver{
	/* Solve quadratic programming problems using active set approach.
	 *
	 * QP formulation:  
	 *	min		0.5 z^T G z + g^T z
	 *	s.t		lbineq <= AiZ z <= ubineq
	 */

public:
	//Default constructor
	QPSolver(std::string dir);

	// Constructor for C++ interface
	QPSolver(real_t * Li_i, real_t * g_i, real_t * Aineq_i, real_t * lbineq_i,
			real_t * ubineq_i, int_t nz_i, int_t nc_i,
			real_t tolMin_i = 1e-9, real_t tolMax_i = 1e-5, int_t MAXITER_i = 50, int_t iterRelax_i = 25);

	~QPSolver(); // destructor
	
	// solve QP with inequality constraints using incremental Active Set approach
	void	solve();

	// number of active set iterations
	int_t	getIterNumber() const {return iter;}
	
	// get exit flag
	int_t	getExitFlag() const { return exitFlag; }
	
	// update problem for warmstart
	void	updateProblem();
private:
	// Perform initialization
	void	initialize();
	// perform standard active set approach (when lambda>0)
	void	activeSetIterations(const int_t extra_idx=0);

	// calculate Lagrange multipliers for the current active set
	void	calcLambda();
	
	// add constraint to the active set
	void	addConstraint(const  int_t viol_idx);
		
	// calculate the error for a particular constraint
	void	calculateError(const int_t idx,const real_t *const x, real_t *const err) const;

	// params
	int_t	MAXITER;			// Max iterations for Active Set approach

			
	// dimensions
	int_t	iter;				// number of active set iterations
	
	// variables for active set iterations
	real_t	*z_sol,				// solution of previous iteration
			*delta,				// difference between current and previous solutions
			*a_del;				// AiZ(inactive)*delta

	int_t	*indices;			// inactive indices which move in direction of constraint

	int_t	exitFlag;			// 0 for solved
								// -1 for infeasible IC (comes from kickout constraint)
								// -2 for maxIter in Primal active set method
								// -3 for maxIter in solver

protected:
	real_t	*Li;				// Li = inv(chol(G,lower))
	
	int_t	nc,					// Total number of inequality constraints;
			nz;					// number of variables in z formulation = m*s-n 

	real_t	*z,					// variable in reduced problem
			*g,					// linear part of cost function in formulation 2
			*AiZ,				// inequality constraints	
			*lbineq,			// lower bound of inequality constraints:			
			*ubineq,			// upper bound of inequality constraints		

			*temp_nz,			// temporary variables
			*temp_nz2,

			*LiTLi,				// LiTLi = Li^T * Li	(repeatedly used, precomputed)
			*lambda,			// lagrange multipliers	of active set
		
			tolMin,				// minimum tolerance used for checking constraints
			tolMax,				// maximum tolerance used for checking constraints
			TOL;				// Absolute tolerance for active set approach 

	int_t	iterRelax,			// iteration at which tolerance is relaxed to maxTol
			viol_idx;			// index of maximum violation; negative index for lb
								// add 1 to absolute index because 0 and -0 are same
	
	// Matrix with active constraints 
	ActiveConstraints *activeCons;
								
	bool	viol;				// violation status
	// check constraints for the reduced problem
	virtual void	checkConstraints();

	// calculate the solution for the current active set
	virtual void	calc_z();

};

