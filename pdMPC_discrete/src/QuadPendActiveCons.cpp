#include "QuadPendActiveCons.h"

QuadPendActiveCons::QuadPendActiveCons(const real_t *const AiZ, const real_t *const lb,
	const real_t *const ub, const real_t *const Li, const int_t nz, const int_t* const soft, const real_t &mu) :
	ActiveConstraints(AiZ, lb, ub, Li, nz), soft_indices(soft), sq_mu(Utils::getSqrt(mu))
{
	m_nz_ori = nz - 1;
	w_b = new real_t[m_nz_ori]();	// initialize with max size = original variable size

}

QuadPendActiveCons::~QuadPendActiveCons() {
	delete[] w_b;
}


void QuadPendActiveCons::addConstraint(const int_t viol_idx) {
	// Update Q,R,active set:
	// Adds QT*Li*viol_lhs' at the end of matrix R (as a column)

	// update w_b
	update_wb(viol_idx);

	// col = QT*Li*viol_lhs';
	if (viol_idx>0)
		Utils::MatVecMult(m_Li, &m_AiZ[(viol_idx - 1)*m_nz_ori], m_temp_nz, m_nz_ori, m_nz_ori);		// temp = Li*viol'
	else {
		Utils::MatVecMult(m_Li, &m_AiZ[(-viol_idx - 1)*m_nz_ori], m_temp_nz, m_nz_ori, m_nz_ori);		// temp = Li*viol'
		Utils::ScalarVectorMult(m_temp_nz, -1, m_nz_ori);
	}
	m_temp_nz[m_nz-1] = w_b[getActiveSetSize()] / sq_mu;
	Utils::MatTVecMult(m_Q, m_temp_nz, m_temp_nz2, m_nz, m_nz);							// temp2 = QT*temp

																						// update R matrix
	Rmat->updateR(m_temp_nz2);

	// update m_active
	active.incrementSet(viol_idx);

	// update Q: Qnew = Qold*GqT
	Rmat->multiplyGqT(m_Q);

}

void QuadPendActiveCons::update_wb(int_t viol_idx) {
	// w_b = soft*abs(w);
	if (viol_idx> 0)
	{	// ub
		if (m_ub[viol_idx - 1] > 0) { w_b[getActiveSetSize()] = soft_indices[viol_idx - 1] * m_ub[viol_idx - 1]; }
		else						{ w_b[getActiveSetSize()] = -soft_indices[viol_idx - 1] * m_ub[viol_idx - 1]; }
	}
	else
	{	// lb
		if (m_lb[-viol_idx - 1] > 0) { w_b[getActiveSetSize()] = soft_indices[-viol_idx - 1] * m_lb[-viol_idx - 1]; }
		else						 { w_b[getActiveSetSize()] = -soft_indices[-viol_idx - 1] * m_lb[-viol_idx - 1]; }
	}
}

void QuadPendActiveCons::updateActiveSet() {
	// updates the w_b vector when problem is updated
	// w_b = soft*abs(w);

	for (int_t i = 0; i <getActiveSetSize(); ++i) {
		int_t cons_idx = getActiveIndex(i);

		if (cons_idx> 0)
		{	// ub
			if (m_ub[cons_idx - 1] > 0) { w_b[i] = soft_indices[cons_idx - 1] * m_ub[cons_idx - 1]; }
			else						{ w_b[i] = -soft_indices[cons_idx - 1] * m_ub[cons_idx - 1]; }
		}
		else
		{	// lb
			if (m_lb[-cons_idx - 1] > 0) { w_b[i] = soft_indices[-cons_idx - 1] * m_lb[-cons_idx - 1]; }
			else						 { w_b[i] = -soft_indices[-cons_idx - 1] * m_lb[-cons_idx - 1]; }
		}
	}
}

void QuadPendActiveCons::dotProductFunction(const real_t *const lambda, real_t &val) {
	// use to multiply lambda with w_b
	Utils::DotProduct(lambda, w_b, getActiveSetSize(), val);
}

void QuadPendActiveCons::multiplyW_vector(const real_t* const vec1, real_t* const vec2) const {

	for (int_t i = 0; i<active.getSize(); ++i) {
		vec2[i] = 0.0;

		if (active.getIndex(i)>0) { // upper bound, no sign inversion
			for (int_t k = 0; k<m_nz_ori; ++k) {
				vec2[i] += m_AiZ[(active.getIndex(i) - 1)*m_nz_ori + k] * vec1[k];
			}
		}
		else {				// lower bound, sign inversion
			for (int_t k = 0; k<m_nz_ori; ++k) {
				vec2[i] += -m_AiZ[(-active.getIndex(i) - 1)*m_nz_ori + k] * vec1[k];
			}

		}
	}

}

void QuadPendActiveCons::multiplyWT_vector(const real_t* const vec1, real_t* const vec2) const {
	// vec2 must be a zero vector
	for (int_t k = 0; k<m_nz_ori; ++k) {
		vec2[k] = 0.0;
	}

	for (int_t i = 0; i<active.getSize(); ++i) {
		if (active.getIndex(i)>0) { // upper bound, no sign inversion
			for (int_t k = 0; k<m_nz_ori; ++k) {
				vec2[k] += m_AiZ[(active.getIndex(i) - 1)*m_nz_ori + k] * vec1[i];
			}
		}
		else {				// lower bound, sign inversion
			for (int_t k = 0; k<m_nz_ori; ++k) {
				vec2[k] += -m_AiZ[(-active.getIndex(i) - 1)*m_nz_ori + k] * vec1[i];
			}

		}
	}

}