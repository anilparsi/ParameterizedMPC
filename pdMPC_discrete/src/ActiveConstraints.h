#pragma once
#include "DefineSettings.h"
#include "Utils.h"
#include "Rmatrix.h"

class ConstraintSet{
public:
	// constructor
	ConstraintSet(){
		// initialize number of elements to 0
		n = 0;
	}

	// destructor
	~ConstraintSet(){
		
	}

	void incrementSet(const int_t viol_idx){
		indices[n] = viol_idx;
		++n;
	};

	void decrementSet(const int_t idx){
		// delete element
		for(int_t i=idx;i<n-1;++i){
			indices[i] = indices[i+1];
		}
		// update set size
		--n;
	};


	const int_t& getSize()const	{ return n;} ;

	const int_t* getSizePtr() const {return &n;};

	const int_t& getIndex(const int_t idx) const{
		return indices[idx];};
private:
	// size of set
	int_t n;

	// indices in the set
	int_t indices[MAX_VARS];						// need nz to get exact size


};

class ActiveConstraints{
	// W is the matrix containing the active constraints' coefficients
	// Li*W^T = Q*R is the way in which it is stored

public:
	// constructor
	ActiveConstraints(const real_t *const AiZ, const real_t *const lb, 
				const real_t *const ub, const real_t *const Li, const int_t nz); 
	
	// destructor
	~ActiveConstraints(); 

	// multiply matrix W with vec1
	virtual void multiplyW_vector(const real_t *const vec1, real_t *const vec2) const;
	
	// multiply matrix WT with vec1
	virtual void multiplyWT_vector(const real_t *const vec1, real_t *const vec2) const;

	// add vector w to vec1
	void add_w_vector(const real_t *vec1, real_t *const vec2) const;
	
	// vec1 = (R'*R)\vec1;
	void performRTRSub(real_t *const vec1){
		Rmat->performRTRSubstitution(vec1);
	};
	
	// remove one constraint from Q and R
	void removeConstraint(const int_t idx);
	
	// add new constraint to Q and R
	virtual void addConstraint(const int_t viol_idx);

	const int_t& getActiveSetSize() const{
		return active.getSize();
	};

	// Used for check: must be deleted later
	const int_t& getActiveIndex(const int_t idx) const{
		return active.getIndex(idx);
	};

	// returns flag to indicate if the constraint set is linearly dependent
	bool getLD_Flag() {
		return Rmat->getLD_Flag();
	}

	// reset active set: error handling
	void resetActiveSet();

	// update active set for warmstart
	virtual void updateActiveSet() {
		// used in virtual class for slack constraints
	}

	// empty function, used for slack constraints
	virtual void dotProductFunction(const real_t *const lambda, real_t &val) {};

protected:
	
	// member variables
	real_t	*m_Q;							// current Q matrix

	Rmatrix *Rmat;							// R matrix from QR decomposition

	
	// const member variables 
	
	const real_t	*m_AiZ,					// constraint matrix with all constraints
					*m_lb,					// lower and upper bounds
					*m_ub,
					*m_Li;					// from the Çholesky decomposition of the Hessian

	const int_t		m_nz;

	ConstraintSet	active;					// active set

	
	real_t			*m_temp_nz,				// temporary variables
					*m_temp_nz2;
	
};
