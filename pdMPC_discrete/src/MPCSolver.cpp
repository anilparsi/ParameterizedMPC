#include "MPCSolver.h"
#include "Utils.h"

MPCSolver::MPCSolver(std::string dir): QPSolver(dir){
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

	tmp=dir+"/tauk";
	Utils::LoadVec(tmp.c_str(),&tauk,tmp2);
	
	t_star = static_cast<int_t>(tmp2/s);

	tmp=dir+"/b_u";
	Utils::LoadVec(tmp.c_str(),&b_u,tmp2);

	tmp=dir+"/b_l";
	Utils::LoadVec(tmp.c_str(),&b_l,tmp2);
	
	m_np = tmp2;
	m_nw = m_np*s;
	

	eta_u = new real_t[m*s]();
	est_lbErr = new real_t[m_np]();	
	est_ubErr = new real_t[m_np]();
	eta_w = new real_t[m_nw]();
	norm_w = new real_t[m_np]();
	temp_nw = new real_t[m_nw];
	temp_nc = new real_t[nc]();
	lbineq_c = new real_t[nc]();
	ubineq_c = new real_t[nc]();
	u = new real_t[m]();
	
	Utils::VectorCopy(lbineq,lbineq_c,nc);
	Utils::VectorCopy(ubineq,ubineq_c,nc);
}

MPCSolver::~MPCSolver(){
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
	delete[] tauk;
	delete[] b_u;
	delete[] b_l;
}

void MPCSolver::updateMPCProblem(const real_t *const x_IC ){
	// update linear cost g = F*x0;
	Utils::MatVecMult(F,x_IC,g,nz,n);
	
	// update bounds on inequality constraint	 lbineq = lbineq_c - AiC*x0;
	
	Utils::MatVecMult(AiC,x_IC,temp_nc,nc,n);
	for(int_t k=0;k<nc;++k){				// for each constraint
		
		// update bounds
		lbineq[k] = lbineq_c[k] - temp_nc[k];
		ubineq[k] = ubineq_c[k] - temp_nc[k];
		
	}
}

void MPCSolver::solve(const real_t *const x_IC){
	x0 = x_IC;
	// update the parameters depending on x0
	updateMPCProblem(x_IC);

	QPSolver::solve();

	if (!viol){
	// QP solved: convert z to u
	real_t *temp = new real_t[m*s];
			
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

void MPCSolver::checkConstraints(){
	// choose method with less number of variables
	//if(s>nz){	
		QPSolver::checkConstraints();
	//}else{
		//checkConstraints_skip();
	//}
}

void MPCSolver::checkConstraints_skip(){
	//eta_w = C0*x0 + C1*z
	Utils::MatVecMult(C0,x0,temp_nw,m_nw,n);
	Utils::MatVecMult(C1,z,eta_w,m_nw,nz);
	Utils::VectorAdd(eta_w,temp_nw,eta_w,m_nw);
	
	for (int i=0;i<m_np;++i){
		norm_w[i] = Utils::VectorNorm(&eta_w[i*s],s);
	}

	// val = tau0d'*eta_w;
	
	real_t val = 0.0;
	int_t idx1 = 0;							// index in constraint matrix
	real_t max_error = -INFVAL;

	// calculate the upper and lower bound errors for state constraints at t=0 
	for (int_t i = 0; i < m_np -m; ++i) {	
		// do not increment index because these are not in inequality matrix
		est_ubErr[i] = INFVAL;
		est_lbErr[i] = INFVAL;
	}

	// calculate the upper and lower bound errors for input constraints at t=0 
	for (int_t i = m_np - m; i<m_np; ++i) {
		++idx1;

		Utils::DotProduct(tauk, &eta_w[i*s], s, val);
		est_ubErr[i] = val - b_u[i];
		est_lbErr[i] = b_l[i] - val;

		if (est_ubErr[i]>max_error) {
			max_error = est_ubErr[i];
			viol_idx = idx1;
		}

		if (est_lbErr[i]>max_error) {
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

				}
							
			}
		}
	}

	
	viol = (max_error>TOL);
		
}
void MPCSolver::getControlInputs(real_t *u_out) const{
	for(int i=0;i<m;++i){
	u_out[i] = u[i];
	}
	//printf("input is %f.\n",u[0]);
}