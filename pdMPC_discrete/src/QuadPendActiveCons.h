#pragma once
#include "ActiveConstraints.h"

class QuadPendActiveCons : public ActiveConstraints {
public:
	QuadPendActiveCons(const real_t *const AiZ, const real_t *const lb,
		const real_t *const ub, const real_t *const Li, const int_t nz, const int_t* const soft, const real_t &mu);
	
	// destructor
	~QuadPendActiveCons();

	void addConstraint(const int_t viol_idx) override;

	void updateActiveSet() override;

	void dotProductFunction(const real_t * const lambda, real_t &val) override;

	// update function to use m_nz_ori
	void multiplyW_vector(const real_t * const vec1, real_t * const vec2) const override;

	// update function to use m_nz_ori
	void multiplyWT_vector(const real_t * const vec1, real_t * const vec2) const override;

private:
	int_t	m_nz_ori;				// original size of the variables with slack not appended

	real_t	*w_b;					// variable with slack coefficiens for active constraints

	const int_t	*soft_indices;	// soft indices

	void update_wb(int_t viol_idx);

	const real_t sq_mu;
};