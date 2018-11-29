#include "ActiveConstraints.h"
#include "DefineSettings.h"

#include "Utils.h"
#include "Rmatrix.h"

#include <vector>
#include <cassert>
#include <algorithm>

ActiveConstraints::ActiveConstraints(const real_t *const AiZ, const real_t *const lb,
									 const real_t *const ub, const real_t *const Li, const int_t nz):
									 m_AiZ(AiZ), m_lb(lb), m_ub(ub), m_Li(Li),m_nz(nz)
{
	m_Q			= new real_t [m_nz*m_nz]();
	
	for (int i=0; i<m_nz; ++i){
		m_Q[i*m_nz+i] = 1.0;				// Initialize to identity
	}
	Rmat		= new Rmatrix (m_nz,active.getSizePtr());

	m_temp_nz	= new real_t [m_nz];
	m_temp_nz2	= new real_t [m_nz];
};

ActiveConstraints::~ActiveConstraints(){
	delete[] m_Q;
	delete[] m_temp_nz;
	delete[] m_temp_nz2;
	delete	 Rmat;
}


void ActiveConstraints::addConstraint(const int_t viol_idx){
	// Update Q,R,active set:
	// Adds QT*Li*viol_lhs' at the end of matrix R (as a column)

	// col = QT*Li*viol_lhs';
	if(viol_idx>0)
		Utils::MatVecMult(m_Li,&m_AiZ[(viol_idx-1)*m_nz],m_temp_nz,m_nz,m_nz);		// temp = Li*viol'
	else{
		Utils::MatVecMult(m_Li,&m_AiZ[(-viol_idx-1)*m_nz],m_temp_nz,m_nz,m_nz);		// temp = Li*viol'
		Utils::ScalarVectorMult(m_temp_nz,-1,m_nz);
	}
	Utils::MatTVecMult(m_Q,m_temp_nz,m_temp_nz2,m_nz,m_nz);							// temp2 = QT*temp
	
	// update R matrix
	Rmat->updateR(m_temp_nz2);									
	
	// update m_active
	active.incrementSet(viol_idx);

	// update Q: Qnew = Qold*GqT
	Rmat->multiplyGqT(m_Q);

}

void ActiveConstraints::resetActiveSet()
{
	// remove elements from active set
	for (int_t i = getActiveSetSize(); i > 0; --i) {
		active.decrementSet(i-1);
	}

	// set Q matrix to identity
	for (int i = 0; i<m_nz; ++i) {
		for (int j = 0; j<m_nz; ++j) {
			m_Q[i*m_nz + j] = 0.0;
		}
		m_Q[i*m_nz + i] = 1.0;
	}
}

void ActiveConstraints::removeConstraint(const int_t idx){
	// downdate Q, R, active set

	//printf("Index to be removed is %d and # of cons are %d \n",idx,active.getSize(););
	assert(idx<active.getSize() && "Constraint to remove is out of range");
	assert(idx>=0 && "Constraint to remove is out of range");
	
	Rmat->downdateR(idx);

	// update m_active
	active.decrementSet(idx);

	// update Q
	if (active.getSize()==0){ // Initialize to identity
		for (int i=0; i<m_nz; ++i){
			for (int j=0; j<m_nz; ++j){
				m_Q[i*m_nz+j] = 0.0;
			}
			m_Q[i*m_nz+i] = 1.0;				
		}
	}else{
		// update Q: Qnew = Qold*GqT
		Rmat->multiplyGqT(m_Q);
	}
	
}

void ActiveConstraints::multiplyW_vector(const real_t* const vec1, real_t* const vec2) const{
	
	for (int_t i = 0; i<active.getSize(); ++i){
		vec2[i] = 0.0;
		
		if (active.getIndex(i)>0){ // upper bound, no sign inversion
			for (int_t k = 0; k<m_nz; ++k){
				vec2[i] += m_AiZ[(active.getIndex(i)-1)*m_nz+k]*vec1[k];
			}
		}else{				// lower bound, sign inversion
			for (int_t k = 0; k<m_nz; ++k){
				vec2[i] += -m_AiZ[(-active.getIndex(i)-1)*m_nz+k]*vec1[k];
			}
		
		}
	}

}

void ActiveConstraints::multiplyWT_vector(const real_t* const vec1, real_t* const vec2) const{
	// vec2 must be a zero vector
		for (int_t k = 0; k<m_nz; ++k){
			vec2[k] = 0.0;
		}
	
	for (int_t i = 0; i<active.getSize(); ++i){		
		if (active.getIndex(i)>0){ // upper bound, no sign inversion
			for (int_t k = 0; k<m_nz; ++k){
				vec2[k] += m_AiZ[(active.getIndex(i)-1)*m_nz+k]*vec1[i];
			}
		}else{				// lower bound, sign inversion
			for (int_t k = 0; k<m_nz; ++k){
				vec2[k] += -m_AiZ[(-active.getIndex(i)-1)*m_nz+k]*vec1[i];
			}
		
		}
	}

}

void ActiveConstraints::add_w_vector(const real_t *vec1, real_t *const vec2) const{
	for (int_t i = 0; i<active.getSize(); ++i){
		vec2[i] = vec1[i];
		
		if (active.getIndex(i)>0){ // upper bound, no sign inversion
			vec2[i] += m_ub[(active.getIndex(i)-1)];
		}else{				// lower bound, sign inversion
			vec2[i] += -m_lb[(-active.getIndex(i)-1)];
		}
	}
}
