#include "QuadMPCSolver.h"
#include "Utils.h"
#include <cassert>

QuadMPCSolver::QuadMPCSolver(std::string dir): QPSolver(dir){
	std::string tmp;
	int_t tmp2;
			
	// Basic paramters
	{
		tmp=dir+"/params";
		int_t nparams;
		real_t *tmpvec;
		
		Utils::LoadVec(tmp.c_str(),&tmpvec,nparams);
		
		n = (int_t)tmpvec[5];
		m = (int_t)tmpvec[6];
		s = (int_t)tmpvec[7];

		delete [] tmpvec;
	}
	
	tmp=dir+"/AiC";
	Utils::LoadVec(tmp.c_str(),&AiC,tmp2);			// tmp2 = nc*n

	tmp=dir+"/C";
	Utils::LoadVec(tmp.c_str(),&C,tmp2);			

	tmp=dir+"/eta2u";
	Utils::LoadVec(tmp.c_str(),&eta2u,tmp2);		
	
	tmp=dir+"/Z";
	Utils::LoadVec(tmp.c_str(),&Z,tmp2);			
	
	tmp=dir+"/F";
	Utils::LoadVec(tmp.c_str(),&F,tmp2);
	
	tmp=dir+"/C0";
	Utils::LoadVec(tmp.c_str(),&C0,tmp2);
	
	tmp=dir+"/C1";
	Utils::LoadVec(tmp.c_str(),&C1,tmp2);
	
	tmp=dir+"/time_indices";
	Utils::LoadVec(tmp.c_str(),&time_indices,tmp2);

	tmp=dir+"/norms";
	Utils::LoadVec(tmp.c_str(),&norms,tmp2);
	
	tmp=dir+"/Mi";
	Utils::LoadVec(tmp.c_str(),&Mi,tmp2);

	tmp=dir+"/tauk";
	Utils::LoadVec(tmp.c_str(),&tauk,tmp2);
	
	t_star = static_cast<int_t>(tmp2/s);

	tmp=dir+"/b_u";
	Utils::LoadVec(tmp.c_str(),&b_u,tmp2);

	tmp=dir+"/b_l";
	Utils::LoadVec(tmp.c_str(),&b_l,tmp2);
	
	m_np = tmp2;
	m_nw = m_np*s;
	
	eta_u		= new real_t[m*s]();
	est_lbErr	= new real_t[m_np]();	
	est_ubErr	= new real_t[m_np]();
	eta_w		= new real_t[m_nw]();
	norm_w		= new real_t[m_np]();
	temp_nw		= new real_t[m_nw];
	temp_nc		= new real_t[nc]();
	lbineq_c	= new real_t[nc]();
	ubineq_c	= new real_t[nc]();
	u			= new real_t[m]();
	
	Utils::VectorCopy(lbineq,lbineq_c,nc);
	Utils::VectorCopy(ubineq,ubineq_c,nc);

	OL_count = 0;
}

QuadMPCSolver::~QuadMPCSolver(){
	delete[] eta_u;
	delete[] est_lbErr;
	delete[] est_ubErr;
	delete[] eta_w;
	delete[] norm_w;
	delete[] temp_nw;
	delete[] temp_nc;
	delete[] lbineq_c;
	delete[] ubineq_c;
	delete[] u;


	delete[] AiC;
	delete[] C;
	delete[] eta2u;
	delete[] Z;
	delete[] F;
	delete[] C0;
	delete[] C1;
	delete[] time_indices;
	delete[] norms;
	delete[] Mi;
	delete[] tauk;
	delete[] b_u;
	delete[] b_l;
}

void QuadMPCSolver::updateMPCProblem(const real_t *const x_IC , const real_t *const u_old){
	// update linear cost g = F*x0;
	Utils::MatVecMult(F,x_IC,g,nz,n);
	
	// update bounds on inequality constraint lbineq = lbineq_c - AiC*x0;
	Utils::MatVecMult(AiC,x_IC,temp_nc,nc,n);
	for(int_t k=0;k<nc;++k){				// for each constraint
		
		// update bounds
		lbineq[k] = lbineq_c[k] - temp_nc[k];
		ubineq[k] = ubineq_c[k] - temp_nc[k];
		
	}
	if(u_old){
		// update the first m constraints: input constraints due to omega rate
		// lbineq = lbineq + Mi * [omega[-1]; 0]
			
		for (int_t i = 0; i<m; ++i){
			temp_nc[i] = 0;
			for (int_t k = 0; k<m-1; ++k){
				temp_nc[i] += Mi[i*m+k]*u_old[k];
			}
		}
		Utils::VectorAdd(lbineq,temp_nc,lbineq,m);
		Utils::VectorAdd(ubineq,temp_nc,ubineq,m);
	}
}

void QuadMPCSolver::solve(const real_t *const x_IC, const real_t *const u_old){
	x0 = x_IC;
	// update the parameters depending on x0
	updateMPCProblem(x_IC,u_old);
	QPSolver::solve();
	if (QPSolver::getExitFlag()<0)
	{
		// problem not solved
		// copy previous solution
		// Utils::VectorCopy(u_old, u, m);

		// increment open loop count
		++OL_count;

		// get previous open loop solution
		for(int i = 0; i < m; ++i) {
			Utils::DotProduct(&eta_u[i*s], &tauk[s*OL_count], s, u[i]);			
		}
		
			
		return;
	}else{
		// QP solved: convert z to u
		real_t *temp = new real_t[m*s];
		
		// reset open loop count
		OL_count = 0;			
			
		// eta_z = C*x0 + Z*zk_sol;
		// eta_u = C(ns+1:end,:)*x0 + Z(ns+1:end,:)*zk_sol

		Utils::MatVecMult(&C[n*s*n],x_IC,temp,m*s,n);
		Utils::MatVecMult(&Z[n*s*nz],z,eta_u,m*s,nz);
		Utils::VectorAdd(eta_u,temp,eta_u,m*s);

		// u = eta2u*eta_z(ns+1:end) = eta2u * eta_u;
		Utils::MatVecMult(eta2u,eta_u,u,m,m*s);

		delete[] temp;
	}
}

void QuadMPCSolver::checkConstraints(){
	// choose method with less number of variables
	if(s<nz){	//
		QPSolver::checkConstraints();
	}else{
		checkConstraints_skip();
	}
}

void QuadMPCSolver::checkConstraints_skip(){
	//eta_w = C0*x0 + C1*z
	Utils::MatVecMult(C0,x0,temp_nw,m_nw,n);
	Utils::MatVecMult(C1,z,eta_w,m_nw,nz);
	Utils::VectorAdd(eta_w,temp_nw,eta_w,m_nw);
	
	for (int i=0;i<m_np;++i){
		norm_w[i] = Utils::VectorNorm(&eta_w[i*s],s);
	}

	// val = tau0d'*eta_w;
	
#ifdef DEBUG_SKIP
	eval = 0;
#endif
	
	viol_idx = 0;
	real_t max_error = -INFVAL;
	// calculate the upper and lower bound errors for omega_rate
	for(int_t i=0;i<m;++i){
		real_t prod = 0.0;
		
		for(int j=0;j<nz;++j){
			prod += AiZ[i*nz+j]*z[j];
		}
		
		// errors
		real_t	e1 = prod-ubineq[i];
		
		if (e1>max_error)
		{	// update max error
			viol_idx = i+1;
			max_error = e1;
		}
		 
		e1 = lbineq[i] - prod;
		if (e1>max_error)
		{	// update max error
			viol_idx = -i-1;
			max_error = e1;
		}
	}

	
	real_t val = 0.0;
	int_t idx1 = m;							// index in constraint matrix
	
	// calculate the upper and lower bound errors for t=0
	for(int_t i =0; i<m_np; ++i){
		++idx1;

		Utils::DotProduct(tauk,&eta_w[i*s],s,val);
		est_ubErr[i] = val - b_u[i];
		est_lbErr[i] = b_l[i] - val;

		if(est_ubErr[i]>max_error){
			max_error = est_ubErr[i];
			viol_idx = idx1;
		}

		if(est_lbErr[i]>max_error){
			max_error = est_lbErr[i];
			viol_idx = -idx1;
		}
	}
	

	// skip constraints loop
	for(int i = 0; i < t_star; ++i){
		for(int k = 0; k<m_np; ++k){
				
			if(time_indices[(i+1)*m_np+k]>0){ // constraint is in non-redundant set
				++idx1;
				
				// update estimates
				est_ubErr[k] += norms[i]*norm_w[k];
				est_lbErr[k] += norms[i]*norm_w[k];

				if (est_ubErr[k] > max_error || est_lbErr[k] > max_error){
					// estimate crosses bound: find exact value
					Utils::DotProduct(&tauk[(i+1)*s],&eta_w[k*s],s,val);
					est_ubErr[k] = val - b_u[k]; 
					est_lbErr[k] = b_l[k] - val;

					if(est_ubErr[k]>max_error){
						max_error = est_ubErr[k];
						viol_idx = idx1;
					}

					if(est_lbErr[k]>max_error){
						max_error = est_lbErr[k];
						viol_idx = -idx1;
					}
#ifdef DEBUG_SKIP
	++eval;
#endif
				}
							
			}
		}
	}
#ifdef DEBUG_SKIP
	printf("%d constraints were evaluated \n", eval);
#endif
	
	viol = (max_error>TOL);
		
}
void QuadMPCSolver::getSolutionCopy(real_t *const u_out) const{
	for(int i=0;i<m;++i){
	u_out[i] = u[i];
	}
	//printf("input is %f.\n",u[0]);
}