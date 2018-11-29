#include "Rmatrix.h"
#include "Utils.h"
#include "DefineSettings.h"
#include <cmath>

Rmatrix::Rmatrix(const int_t nz, const int_t *const nac):
		m_nz(nz), m_nac(nac),TOL(1e-15){	// hardcoded tolerance for linear dependency

	temp_nz		= new real_t [m_nz];
	temp_nz2	= new real_t [m_nz];
	temp_nznz	= new real_t [m_nz*m_nz];
	m_R			= new real_t [static_cast<int_t>(m_nz*m_nz+m_nz)/2];
	Gq			= new real_t [m_nz*m_nz];
};

Rmatrix::~Rmatrix(){
	delete[] temp_nz;
	delete[] temp_nz2;
	delete[] temp_nznz;
	delete[] m_R;
	delete[] Gq;

};
void Rmatrix::updateR(real_t *const vec1){

	for (int i=0; i<m_nz; ++i){
		// convert Gq into an identity matrix
		for(int j = 0; j<m_nz; ++j){
			Gq[i*m_nz+j] = 0.0;	
		}
		Gq[i*m_nz+i] = 1.0;			
	}

	// convert vec1 to have *nac elements
	for (int i=m_nz-2;i>*m_nac-1;--i){
		givens(vec1[i], vec1[i+1]);
		givensMatUpdate(i,i+1);									
		vec1[i] = giv_c*vec1[i] + giv_s*vec1[i+1];				// update ìth element in vec1
		vec1[(i+1)] = 0.0;						// set the i+1th element to zero
	}

	// add updated vec to R	   
	for (int i = 0; i<*m_nac + 1; ++i) {
		m_R[(*m_nac* *m_nac + *m_nac) / 2 + i] = vec1[i];
	}

	
}

void Rmatrix::downdateR(const int_t idx){
	
	// convert Gq into an identity matrix
	for (int i=0; i<m_nz; ++i){
		for(int j = 0; j<m_nz; ++j){
			Gq[i*m_nz+j] = 0.0;	
		}
		Gq[i*m_nz+i] = 1.0;			
	}
	
	// downdate R (edit columns to right of idx: remove last elements)
	for(int i = idx+1; i<*m_nac; ++i){		// each column
		int_t i1 = (i*i + 3*i)/2;
		givens(m_R[i1-1],m_R[i1]);		// R[i-1,i], R[i,i]
		givensMatUpdate(i-1,i);									
		
		for(int j=i;j<*m_nac;++j){			// columns affected by givens 
			int_t k = (j*j + j)/2;
			m_R[k+i-1] = giv_c*m_R[(j*j + j)/2+i-1] + giv_s*m_R[k+i];					// R[i-1,j]
			m_R[k+i] = -giv_s/giv_c*m_R[k+i-1] + (giv_c + giv_s*giv_s/giv_c)*m_R[k+i];	// R[i,j]
		}
		
	}
	
	// remove the column at idx
	for(int j=idx;j<*m_nac-1;++j){	// each column
		int_t k = (j*j + j)/2;
		for(int i = 0;i<j+1;++i){		//each row
			m_R[k+i] = m_R[k+j+1+i];			// swap elements column wise R[i,j] = R[i,j+1]
		}
	}

	for(int i = 0;i<*m_nac;++i){	//each row
		m_R[(*m_nac* *m_nac - *m_nac)/2 + i] = 0.0;						// set elements in last active column to zero	
	}
	
}

void Rmatrix::givens(const real_t x, const real_t y){
	real_t absx = x>=0?x:-x;
	if(x==0.0){
		giv_c = 0.0;
		giv_s = 1.0;
	}else{
		real_t nrm = sqrt(x*x+y*y);
		giv_c = absx/nrm;
		giv_s = x/absx*y/nrm;
	}
}

void Rmatrix::givensMatUpdate(const int_t r1, const int_t r2){

	// Updtes GqT using the givens coefficients
	for(int_t i=0; i<m_nz; ++i){
		temp_nz[i] = giv_c*Gq[r1*m_nz+i] + giv_s*Gq[r2*m_nz+i];
		temp_nz2[i] = -giv_s*Gq[r1*m_nz+i] + giv_c*Gq[r2*m_nz+i];
	}

	Utils::VectorCopy(temp_nz, &Gq[r1*m_nz],m_nz);	// copy row 1
	Utils::VectorCopy(temp_nz2, &Gq[r2*m_nz],m_nz);	// copy row 2
	
}

void Rmatrix::multiplyGqT(real_t *const Q){
	// Qnew = Qold*GqT
	real_t val = 0.0;
	for(int_t i =0; i<m_nz; ++i){					// each row
		for(int_t k =0; k<m_nz; ++k){				// each column
			Utils::DotProduct(&Q[i*m_nz],&Gq[k*m_nz],m_nz,temp_nznz[i*m_nz+k]);
		}
	}

	// need this extra copy of the whole matrix
	Utils::MatrixCopy(temp_nznz,Q,m_nz,m_nz);

}


void Rmatrix::performRTRSubstitution(real_t *const vec1){
	// vec1 = (R'*R)\vec1;

	//Utils::ForwardSubstitution(m_RT,vec1,*m_nac,m_nz);					// lambda = (RT)^-1 * lambda
	for (int i=0; i<*m_nac; ++i){	//each column
		int_t k = (i*i+i)/2;
		for (int j=0;j<i; ++j){		// each row until j
			vec1[i] -= m_R[k+j]*vec1[j];
		}
		vec1[i] = vec1[i]/m_R[k+i];
	}

	//Utils::BackwardSubstitution(m_R,vec1,m_nz,*m_nac);					// lambda = R^-1 * (RT)^-1 * lambda = (R'*R)\lambda
	for (int i=*m_nac-1; i>=0; --i){  // each row
		for (int j=*m_nac-1;j>i; --j){ // each column
			vec1[i] -= m_R[(j*j+j)/2+i]*vec1[j];
		}
		vec1[i] = vec1[i]/m_R[(i*i+3*i)/2];
	}
}

bool Rmatrix::getLD_Flag()
{	// returns flag to indicate if the constraint set is linearly dependent
	/* implemented as a calculation every time it is called than a member variable, because
	 * this is efficient compared to evaluating this at every update and downdate for the 
	 * QP solver algorithm. If number of calls are high, member variable is better.
	 */
	int_t idx = (*m_nac* *m_nac + *m_nac) / 2 - 1;	// index of last diagonal element
	for (int_t i = *m_nac-1; i >= 0; --i) {
		if (m_R[idx]<TOL && -m_R[idx]<TOL) {
			return true;
		}
		idx = idx - i-1;							// update index to previous diagonal element
	}

	return false;
}
