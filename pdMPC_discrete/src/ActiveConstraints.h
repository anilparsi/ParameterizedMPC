#pragma once
#include "DefineSettings.h"
#include "Utils.h"
#include "Rmatrix.h"
/// This class contains the indices of constraints which are active.
/// It is updated whenever the active set is changed
class ConstraintSet{
public:
	/// constructor
	ConstraintSet(){
		// initialize number of elements to 0
		n = 0;
	}

	/// destructor
	~ConstraintSet(){
		
	}

	/// add a constraint index to the set
	void incrementSet(const int_t viol_idx){
		indices[n] = viol_idx;
		++n;
	};

	/// remove a constraint index from the set
	void decrementSet(const int_t idx){
		// delete element
		for(int_t i=idx;i<n-1;++i){
			indices[i] = indices[i+1];
		}
		// update set size
		--n;
	};

	/// get the size of the constraint set
	const int_t& getSize()const	{ return n;} ;

	/// get a pointer to the size of the constraint set
	const int_t* getSizePtr() const {return &n;};

	/// get the index of the constraint in the constraint set
	const int_t& getIndex(const int_t idx) const{
		return indices[idx];};
private:
	// size of set
	int_t n;

	// indices in the set
	int_t indices[MAX_VARS];						// need nz to get exact size


};

/// This class performs matrix operations involving the active set. 
/// Constraints can be added and removed from the active set, and the active set matrix can be multiplied with vectors.
/// W is the matrix containing the active constraints' coefficients
/// w is the vector containing the active bounds on these constraints
class ActiveConstraints{
	// Li*W^T = Q*R is the way in which it is stored

public:
	/// constructor
	ActiveConstraints(const real_t *const AiZ, const real_t *const lb, 
				const real_t *const ub, const real_t *const Li, const int_t nz); 
	
	/// destructor
	~ActiveConstraints(); 

	/// multiply matrix W with vec1
	virtual void multiplyW_vector(const real_t *const vec1, real_t *const vec2) const;
	
	/// multiply matrix WT with vec1
	virtual void multiplyWT_vector(const real_t *const vec1, real_t *const vec2) const;

	/// add vector w to vec1
	void add_w_vector(const real_t *vec1, real_t *const vec2) const;
	
	/// vec1 = (R'*R)\vec1;
	void performRTRSub(real_t *const vec1){
		Rmat->performRTRSubstitution(vec1);
	};
	
	/// remove one constraint from Q and R
	void removeConstraint(const int_t idx);
	
	/// add new constraint to Q and R
	virtual void addConstraint(const int_t viol_idx);

	const int_t& getActiveSetSize() const{
		return active.getSize();
	};

	/// returns the index of the active constraint from the original constraint matrix
	const int_t& getActiveIndex(const int_t idx) const{
		return active.getIndex(idx);
	};

	/// returns flag to indicate if the constraint set is linearly dependent
	bool getLD_Flag() {
		return Rmat->getLD_Flag();
	}

	/// reset active set: error handling
	void resetActiveSet();

	/// update active set for warmstart
	virtual void updateActiveSet() {
		/// used in virtual class for slack constraints
	}

	/// empty function, used for slack constraints
	virtual void dotProductFunction(const real_t *const lambda, real_t &val) {};

protected:
	
	// member variables
	real_t	*m_Q;							///< current Q matrix

	Rmatrix *Rmat;							///< R matrix from QR decomposition

	
	// const member variables 
	
	const real_t	*m_AiZ,					///< constraint matrix with all constraints
					*m_lb,					///< lower and upper bounds
					*m_ub,
					*m_Li;					///< from the �holesky decomposition of the Hessian

	const int_t		m_nz;

	ConstraintSet	active;					///< active set

	
	real_t			*m_temp_nz,				// temporary variables
					*m_temp_nz2;
	
};
