#pragma once

#include "QPSolver.h"
/*! \class MPCSolver
 * \brief This class is be used to solve the QP problems encountered in parameterized 
 * model predictive control (pdMPC).
 * 
 */
class MPCSolver: public QPSolver{
public:  
	/*!
	 * Constructor to implement pdMPC solver with MATLAB interface
	 * \param dir contains the address of the directory with required matrices to solve the pdMPC problem
	 * in .txt files.
	 */
	MPCSolver(std::string dir);

	/// destructor
	~MPCSolver(); 
	
	/// update implementation of solve
	void	solve(const real_t *const x_IC);
	
	/// returns control inputs
	void	getControlInputs(real_t *u_out) const;

	/// returns size of control inputs
	int_t	getNumberOfOutputs() const {return m;}
private:
	/*!
	 * Updates the QP which has to be solved based on the current state of the system.
	 */
	void	updateMPCProblem(const real_t *const x_IC);

	/// update implementation of check constraints
	virtual void checkConstraints() override;

	/// to implement skip constraints method
	void checkConstraints_skip();


	real_t	*Z,					///< from qr decomposition of Aeq
			*C,					///< C = inv(Y)*R*D;: constant for the problem
			*AiC,				///< AiC = Aineq*C
			*F,					///< F = Z'*H*C

			/// Matrices for checkConstraintsSkip
			/// Cs = kron(Cons, eye(s)); (Cons is constraint matrix at each time step)
			*C0,				///< C0	 = Cs*C
			*C1;				///< C1	 = Cs*Z
	


	const real_t  *x0;				///< initial condition for solving optimization

	real_t 	*u,					///< control input
			*tauk,				///< tau vector with size t_star*s
			*b_l,				///< fixed bounds on Cxu for one time step.
			*b_u,
			*eta_u,				///< parameter vector for input variables
								///< eta_z = [eta_x; eta_u];
			*eta2u,				///< conversion matrix from eta to u: kron(eye(m),tau0d')
			
			*temp_nc,			///< temporary variable
			*temp_nw,			///< temporary variable
			
			*lbineq_c,			///< lower bound of inequality constraints			
			*ubineq_c,			///< upper bound of inequality constraints	

			*eta_w,				///< eta_w = kron(Cxu,eye(s))*eta_z;
			*norm_w,			///< norm of each part of eta_w
			*norms;				///< norms of tauk^T*(Md-I);	(k from 0 to t_star)

			
	
	int_t	n,					///< number of states
			m,					///< number of inputs
			s,					///< number of basis funcs
			m_np,				///< number of constraints for single time step (Cxu)
			m_nw;				///< nw = np*s: number of variables in eta_w formulation

	real_t	*est_lbErr,			///< estimates of error
			*est_ubErr;

	int_t	t_star,				///< Number of time steps used in maximal output admissible set
			*time_indices;		///< Time based indices of active constraints: list of active 
								///< constraints at each time step from 0 to t_star

};