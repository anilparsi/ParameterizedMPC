#pragma once
#include "DefineSettings.h"
#include "Utils.h"
#include "Rmatrix.h"
/*!
 * \brief This class contains the indices of constraints which are active.
 * It is updated whenever the active set is changed
 */
class ConstraintSet{
public:
	/// constructor
	ConstraintSet(){		
		// initialize number of elements to 0
		n = 0;
	}

	/// destructor
	~ConstraintSet(){}

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

/*! 
 * \brief This class performs matrix operations involving the active set. 
 * 
 * Constraints can be added and removed from the active set, and the active set matrix can be multiplied with vectors. <br>
 * W is the matrix containing the active constraints' coefficients. <br>
 * w is the vector containing the active bounds on these constraints. <br>
 * Matrix operations to be performed on the matrix W and the vector w are performed through this class. 
 * However, the matrix W is stored in the form of the QR decomposition of Li*W^T. The variable Rmat is used to
 * store the R matrix from the QR decomposition. 
 * 
 */
class ActiveConstraints{
	

public:
	/*! \brief constructor
	 * 
	 * \param AiZ contains the coefficients of inequality constraints
	 * \param lb contains the lower bounds of inequality constraints
	 * \param ub contains the containing upper bounds of inequality constraints
	 * \param 
	 */
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
	
	/// vec1 = inv(R'*R)*vec1;
	void performRTRSub(real_t *const vec1){
		Rmat->performRTRSubstitution(vec1);
	};
	
	/// remove one constraint from Q and R
	void removeConstraint(const int_t idx);
	
	/// add new constraint to Q and R
	virtual void addConstraint(const int_t viol_idx);

	/// returns the current active set size
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


protected:
	
	// member variables
	real_t	*m_Q;							///< current Q matrix

	Rmatrix *Rmat;							///< R matrix from QR decomposition

	
	// const member variables 
	
	const real_t	*m_AiZ,					///< constraint matrix with all constraints
					*m_lb,					///< lower and upper bounds
					*m_ub,
					*m_Li;					///< from the Çholesky decomposition of the Hessian

	const int_t		m_nz;

	ConstraintSet	active;					///< active set

	
	real_t			*m_temp_nz,				// temporary variables
					*m_temp_nz2;
	
};
