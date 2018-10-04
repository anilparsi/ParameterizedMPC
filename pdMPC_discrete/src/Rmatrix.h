#pragma once

#include "DefineSettings.h"

class Rmatrix{
public:
	// constructor
	Rmatrix(const int_t nz, const int_t *const nac);

	//destructor
	~Rmatrix();

	void updateR(real_t *const vec1);
	void downdateR(const int_t idx);
	void performRTRSubstitution(real_t *const vec1);

	// return flag to indicate if the constraint set is linearly dependent
	bool getLD_Flag();
	// updates Q matrix: Q_new = Q_old*GqT
	void multiplyGqT(real_t *const Q);
private:
	// calculate givens coefficients
	void givens(const real_t x, const real_t y);
	
	// update the current givens matrix at positions r1 and r2
	void givensMatUpdate(const int_t r1, const int_t r2);
	
	real_t	*m_R,										// R matrix
			
			giv_c,										// givens cos
			giv_s,										// givens sin
			*Gq,										// givens matrix product

			// temporary variables
			*temp_nz,
			*temp_nz2,
			*temp_nznz;
									
	const int_t m_nz;
	const int_t *const m_nac;
	const real_t TOL;									// Tolerance to check linear dependence of active set: set to large values for safety (1e-5)
};