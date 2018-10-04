#pragma once
#include "QuadMPCSolver.h"

class QuadPendMPCSolver : public QuadMPCSolver {
public: 
	// Default constructor
	QuadPendMPCSolver(std::string dir);

	// destructor
	~QuadPendMPCSolver();

private:

	real_t	ep,				// slack variable 
			mu,				// cost term for slack
			*w_b,			// vector denoting softness of active constraints

			*temp_wb;		// temporary storage of w_b


	int_t	*soft_indices;	// indices of soft constraints in AiZ matrix (loaded)
			
	void calc_z() override;
	
	void checkConstraints() override;

	void checkConstraints_skip() override;

	void updateMPCProblem(const real_t *const x_IC, const real_t *const u_old) override;
};


