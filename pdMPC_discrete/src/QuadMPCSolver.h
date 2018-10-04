#pragma once

#include "QPSolver.h"
class QuadMPCSolver: public QPSolver{
public:
	//Default constructor
	QuadMPCSolver(std::string dir);

	~QuadMPCSolver(); // destructor
	
	// update implementation of solve
	void	solve(const real_t *const x_IC, const real_t *const u_old);
	
	// gives control inputs
	void	getSolutionCopy(real_t *const u_out) const;

	// gives size of control inputs
	int_t	getNumberOfOutputs() const {return m;}
protected:
	virtual void	updateMPCProblem(const real_t *const x_IC, const real_t *const u_old);

	// update implementation of check constraints
	virtual void	checkConstraints() override;

	// to implement skip constraints method
	virtual void	checkConstraints_skip();


	real_t	*Z,					// from qr decomposition of Aeq
			*C,					// C = inv(Y)*R*D;: constant for the problem
			*AiC,				// AiC = Aineq*C
			*F,					// F = Z'*H*C

			*Mi,				// inv(M), where M is the matrix which transforms
								// propeller force vector into [omega_dot;f_tot-f0]
			// Matrices for checkConstraintsSkip
			// Cs = kron(Cons, eye(s)); (Cons is constraint matrix at each time step)
			*C0,				// C0	 = Cs*C
			*C1;				// C1	 = Cs*Z
	


	const real_t  *x0;				// initial condition for solving optimization

	real_t 	*u,					// control input
			*tauk,				// tau vector with size t_star*s
			*b_l,				// fixed bounds on Cxu for one time step.
			*b_u,
			*eta_u,				// parameter vector for input variables
								// eta_z = [eta_x; eta_u];
			*eta2u,				// conversion matrix from eta to u: kron(eye(m),tau0d')
			
			*temp_nc,			// temporary variable
			*temp_nw,			// temporary variable
			
			*lbineq_c,			// lower bound of inequality constraints			
			*ubineq_c,			// upper bound of inequality constraints	

			*eta_w,				// eta_w = kron(Cxu,eye(s))*eta_z;
			*norm_w,			// norm of each part of eta_w
			*norms;				// norms of tauk^T*(Md-I);	(k from 0 to t_star)

			
	
	int_t	n,					// number of states
			m,					// number of inputs
			s,					// number of basis funcs
			m_np,				// number of constraints for single time step (Cxu)
			m_nw;				// nw = np*s: number of variables in eta_w formulation

	real_t	*est_lbErr,			// estimates of error
			*est_ubErr;

	int_t	t_star,				// Number of time steps used in maximal output admissible set
			OL_count,			// Exception handling: applies previous input when problem not solved
			*time_indices;		// Time based indices of active constraints: list of active 
								// constraints at each time step from 0 to t_star

};