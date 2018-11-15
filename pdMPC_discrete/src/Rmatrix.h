#pragma once

#include "DefineSettings.h"
/*!
 * \brief This class is used to store the R matrix in the QR decomposition of active set, 
 * and provides the functions to perform matrix operations on it. 
 */
class Rmatrix{
public:
	/// constructor
	Rmatrix(const int_t nz, const int_t *const nac);

	/// destructor
	~Rmatrix();

	/* \brief add a column to R matrix
	 *
	 * \param vec1 is a pointer to the column to be added
	 */
	void updateR(real_t *const vec1);

	/* \brief remove a column from R matrix
	 *
	 * \param idx is the index of the column to be removed
	 */
	void downdateR(const int_t idx);
	
	/// returns the value inv(R^T*R)*vec1 in the same vector vec1.
	void performRTRSubstitution(real_t *const vec1);

	/// return flag to indicate if the constraint set is linearly dependent
	bool getLD_Flag();

	/* \brief multiplies the input matrix with Gq^T
	 * 
	 * This is used to update the Q matrix in the QR decomposition. Whenever
	 * a column is added or removed from the active set, this is called so that
	 * the Q matrix can be updated along with the R matrix
	 */
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