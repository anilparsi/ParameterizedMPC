
#include "QPSolver.h"
#include "DefineSettings.h"

#include "Utils.h"
#include "ActiveConstraints.h"

#include <direct.h>
#include <cassert>
#include <iostream>
#include <algorithm>

#ifdef TIME
	real_t	tchkStart, tchkTot;
	real_t titStart, titTot=0.0;
	int_t it_count;
#endif

	
#ifdef DEBUG_SKIP
	int_t eval;
#endif

QPSolver::QPSolver(std::string dir){
	// Constructor: Load matrices from directory
	{
		char cCurrentPath[FILENAME_MAX];

		if (!_getcwd(cCurrentPath, sizeof(cCurrentPath))){
			assert(false);
		}

		cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

		printf ("The current working directory is %s \n", cCurrentPath);
	}

	std::string tmp;

	// Basic paramters
	{
		tmp=dir+"/params";
		int_t nparams;
		real_t *tmpvec;
		
		Utils::LoadVec(tmp.c_str(),&tmpvec,nparams);
		
		tolMin = tmpvec[0];
		tolMax = tmpvec[1];
		
		// start with minimum tolerance
		TOL = tolMin;		

		MAXITER = (uint_t)tmpvec[2];
		iterRelax = (int_t)MAXITER / 2;

		nz = (int_t)tmpvec[3];
		nc = (int_t)tmpvec[4];
		delete [] tmpvec;
	}
	// LOAD MPC DATA
	
	int_t tmp2;

	tmp=dir+"/AiZ";
	Utils::LoadVec(tmp.c_str(),&AiZ,tmp2);

	tmp=dir+"/Li";
	Utils::LoadVec(tmp.c_str(),&Li,tmp2);
	
	tmp=dir+"/g";
	Utils::LoadVec(tmp.c_str(),&g,tmp2);
	
	tmp=dir+"/lbineq";
	Utils::LoadVec(tmp.c_str(),&lbineq,tmp2);

	tmp=dir+"/ubineq";
	Utils::LoadVec(tmp.c_str(),&ubineq,tmp2);			

	initialize();
}

QPSolver::QPSolver(real_t * Li_i, real_t * g_i, real_t * Aineq_i, real_t * lbineq_i,
	real_t * ubineq_i, int_t nz_i, int_t nc_i,
	real_t tolMin_i, real_t tolMax_i, int_t MAXITER_i, int_t iterRelax_i)
{
	Li = Li_i;
	g = g_i;
	AiZ = Aineq_i;
	lbineq = lbineq_i;
	ubineq = ubineq_i;

	// get parameters
	nz = nz_i;
	nc = nc_i;
	tolMin = tolMin_i;
	tolMax = tolMax_i;
	MAXITER = MAXITER_i;
	iterRelax = iterRelax_i;

	// start with minimum tolerance
	TOL = tolMin;

	initialize();
}

void QPSolver::initialize()
{
	z = new real_t[nz]();

	lambda = new real_t[nz]();
	activeCons = new ActiveConstraints(AiZ, lbineq, ubineq, Li, nz);

	assert(nz <= MAX_VARS && "nz is less than MAX_VARS");


	temp_nz = new real_t[nz];
	temp_nz2 = new real_t[nz];

	z_sol = new real_t[nz];
	delta = new real_t[nz];
	a_del = new real_t[nz];
	LiTLi = new real_t[nz*nz];
	indices = new int_t[nz + 1];

	// construct LiTLi matrix;
	real_t *temp_nznz = new real_t[nz*nz];
	Utils::MatrixTranspose(Li, temp_nznz, nz, nz);
	Utils::MatrixMult(temp_nznz, Li, LiTLi, nz, nz, nz);
	delete[] temp_nznz;
}

QPSolver::~QPSolver(){
	delete[] z;
	delete[] lambda;
	delete activeCons;	
	delete[] temp_nz;
	delete[] temp_nz2;
	delete[] z_sol;
	delete[] delta;
	delete[] a_del;
	delete[] indices;

	delete[] AiZ;
	delete[] Li;
	delete[] LiTLi;
	delete[] g;
	delete[] lbineq;
	delete[] ubineq;
}



void QPSolver::solve(){
	
	iter = 1;
	exitFlag = 0;
	bool reRunFlag = true;
	while (iter<MAXITER)
	{
		calcLambda();
		calc_z();

		if(Utils::anyPositive(lambda,activeCons->getActiveSetSize())){
			// positive lagrange multiplier: must be removed from active set
#ifdef TIME
	titStart = Utils::getCPUtime();
#endif
			activeSetIterations();
#ifdef TIME
	titTot += Utils::getCPUtime()-titStart;
#endif
		}

#ifdef TIME
	tchkStart = Utils::getCPUtime();
#endif	
		checkConstraints();
#ifdef TIME
	tchkTot += Utils::getCPUtime()-tchkStart;
#endif
		if (viol)
		{	
			addConstraint(viol_idx);

			if (exitFlag < 0) {
				// add exception handler: in derived classes
				// Increase tolerance for the next run. 
				TOL = tolMax;
				activeCons->resetActiveSet();
				return;
			}
			++iter;
		}
		else
		{	//QP solved
			exitFlag = 0;	
			TOL = tolMin;			// Reset tolerance to tight
			return;
		}

		if (iter == iterRelax) {
			// relax constraint tolerance and try to solve the problem
			TOL = tolMax;
			activeCons->resetActiveSet();			
		}
	}

	// maximum iterations reached in activeSet 
	exitFlag = -3;
	// simplify problem
	TOL = tolMax;
	activeCons->resetActiveSet();
}

void QPSolver::calcLambda(){
	if (activeCons->getActiveSetSize()>0){
		
		// lambda =  (R'*R)\(bineq(active) + AiZ(active,:)*(Li'*Li*g));
		Utils::MatVecMult(LiTLi,g,temp_nz2,nz,nz);				// temp = (Li'*Li*g)
		activeCons->multiplyW_vector(temp_nz2,temp_nz);			// temp_nz = AiZ(active,:)*(Li'*Li*g)
		activeCons->add_w_vector(temp_nz,lambda);				// lambda = bineq(active) + AiZ(active,:)*(Li'*Li*g)
	
		//	lambda = (R'*R)\lambda;
		activeCons->performRTRSub(lambda);
		
	}
}

void QPSolver::calc_z(){
	if (activeCons->getActiveSetSize()==0){
		// no active constraints: z = -LiTLig;
		for (int_t i = 0; i<nz; ++i){		
			z[i] = 0;
			for (int_t k = 0; k<nz; ++k){
				z[i] -= LiTLi[i*nz+k]*g[k];				
			}
		}

	}else{
		// xk = LiTLi*(WT*lambda-g)	
		activeCons->multiplyWT_vector(lambda,temp_nz);			// temp = WT*lambda
		Utils::VectorSubstract(temp_nz,g,temp_nz,nz);			// xk = WT*lambda -g
		Utils::MatVecMult(LiTLi,temp_nz,z,nz,nz);				// xk = LiTLi*(WT*lambda-g)
		
	}
}
// check constraints for the reduced problem

void QPSolver::checkConstraints(){
	
	viol_idx = 0;
	real_t max_error = -INFVAL;
	for (int i=0;i<nc;++i){
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
	viol = (max_error>TOL);

}

void QPSolver::updateProblem()
{
	// update the active set for warm start
	activeCons->updateActiveSet();
}



void QPSolver::activeSetIterations(const int_t extra_idx){
	/* This is only called in two cases: 
	 * 1. When the lagrange multiplier for an active cons>0
	 * 2. When there are nz active constraints and we would like to add one more
	 * 
	 * In case 1, extra_idx is 0. In case 2 its the new constraint to be added
	 * Only the active constraints and extra constraint are considered for this function
	 */
		
	// inactive contains the indices of the inactive constraints considered for iterations
	ConstraintSet inactive;
	
	if(extra_idx){
		inactive.incrementSet(extra_idx);
	}
	
	// find index of positive lambda in active set
	int_t index;
	real_t max_val;
	Utils::max_value_idx(lambda,max_val,index,activeCons->getActiveSetSize());

	// delete the positive lambda constraint
	inactive.incrementSet(activeCons->getActiveIndex(index));		
	activeCons->removeConstraint(index);

	// update solution
	Utils::VectorCopy(z,z_sol,nz);
	
	int_t ac_iter =1;
	while (ac_iter<MAXITER){
		calcLambda();	
		calc_z();
		Utils::VectorSubstract(z_sol,z,delta,nz);
		
		if(Utils::VectorInfNorm(delta,nz)<TOL){
			// solution hasn't changed
			Utils::max_value_idx(lambda,max_val,index,activeCons->getActiveSetSize());
			if(max_val<=0){
				// all lagrange multipliers are negative ==> SOLVED!
				break;
			}else{
				// positive lambda exists: delete constraint
				inactive.incrementSet(activeCons->getActiveIndex(index));		
				activeCons->removeConstraint(index);
			}
		}else{
			// z_sol is not the solution for EP
			int_t n_ind = 0;
			Utils::VectorSubstract(z,z_sol,delta,nz);
			for(int i = 0; i<inactive.getSize();++i ){			
				if(inactive.getIndex(i)>0){
					// upper bound
					Utils::DotProduct(&AiZ[(inactive.getIndex(i)-1)*nz],delta,nz,a_del[i]);
				}else{ 
					// lower bound
					Utils::DotProduct(&AiZ[(-inactive.getIndex(i)-1)*nz],delta,nz,a_del[i]);
					a_del[i] = -a_del[i];
				}	
			}
			
			// find all indices which are not active and move in direction of constraint (a_del[i]>TOL)
			for(int_t i=0;i<inactive.getSize();++i){
				if(a_del[i]>TOL){			
					indices[n_ind] = i;
					n_ind++;
				}
			}
			
			// alpha_k = min((w-W[i]*xsol)/AiZ[i]*(zk-zsol))
			real_t alpha_k = 1, temp;
			int_t idx = 0;
			for(int_t i=0; i<n_ind; ++i){
				calculateError(inactive.getIndex(indices[i]),z_sol,&temp); // temp = W[i]*zsol-w;				
				temp = -temp/a_del[indices[i]];					// temp = ((w-W[i]*zsol)/AiZ[i]*(zk-zsol))
				if(temp<alpha_k){
					alpha_k = temp;								// alpha_k = min((w-W[i]*zsol)/AiZ[i]*(zk-zsol))
					idx = indices[i];
				}
			}
			// assert(alpha_k>0 && "Alpha must be positive!");
			
			if(alpha_k<1){
				// blocking constraint exists: add to active set
				activeCons->addConstraint(inactive.getIndex(idx));
				inactive.decrementSet(idx);
				// z_sol = z_sol + alpha_k*delta; (to account for when alpha<1)
				for(int_t i=0;i<nz;++i){
					z_sol[i] += alpha_k*delta[i];	
				}

			}else{
				Utils::VectorCopy(z,z_sol,nz);
			}
		}

		++ac_iter;
		
	}
#ifdef TIME
	it_count+=ac_iter;
#endif		
	if (ac_iter == MAXITER) {	
		// Max Iterations Reached within primal active set
		exitFlag = -2;
	}
}

void QPSolver::addConstraint(const int_t viol_idx){
	int_t	t_nac = activeCons->getActiveSetSize();					// get the initial active set size
	
	
	if (t_nac < nz) { // active set is not full
		activeCons->addConstraint(viol_idx);						// add constraint
	}
	
	// kick out constraints to find new active set
	int_t	temp_out_idx,											// index of constraint kicked out
		temp_idx;												// temporary storage of index to be inserted
	int add_iter = 1;
	
	if (t_nac == nz || activeCons->getLD_Flag()){

		// initialize kicking out constraints
		if (t_nac == nz) {	// for systems with full active set
			temp_out_idx = activeCons->getActiveIndex(t_nac - 1); ;
			activeCons->removeConstraint(t_nac - 1);
			
			// add the violation
			activeCons->addConstraint(viol_idx);
		}else{	// for systems with linearly dependent active set
			// viol is added already: just kickout the last constraint (not viol)
			temp_out_idx = activeCons->getActiveIndex(t_nac - 1); ;
			activeCons->removeConstraint(t_nac - 1);
		}


		while(add_iter<t_nac){

			if (activeCons->getLD_Flag()) {	// check flag again because constraint removed
				// skip current iteration if the active set was not updated
				// kick out the next constraint
				temp_idx = temp_out_idx;
				temp_out_idx = activeCons->getActiveIndex(t_nac - 1 - add_iter);
				activeCons->removeConstraint(t_nac - 1 - add_iter);
				activeCons->addConstraint(temp_idx);
				++add_iter;
				continue;
			}

			calcLambda();
			calc_z();
			real_t err;
			calculateError(temp_out_idx, z, &err);
			if(err>TOL){	// kick out the next constraint(starting from last)
			
				temp_idx = temp_out_idx;
				temp_out_idx =  activeCons->getActiveIndex(t_nac-1-add_iter);
				activeCons->removeConstraint(t_nac-1-add_iter);
				activeCons->addConstraint(temp_idx);

			}else if(Utils::anyPositive(lambda,t_nac)){
				// Perform iterations with one extra constraint!!
#ifdef TIME
	titStart = Utils::getCPUtime();
#endif
				activeSetIterations(temp_out_idx);
#ifdef TIME
	titTot += Utils::getCPUtime()-titStart;
#endif
				break;
			}else{
				// removed constraint successfully
				break;
			}

			
			++add_iter;
		}

	// assert(add_iter<t_nac && "No constraint can be removed.");
		if (add_iter == t_nac) {	// infeasible initial conditions
			exitFlag = -1;		
		}
	}
}



void QPSolver::calculateError(const int_t idx,const real_t *const x, real_t *const err) const{
	if(idx>0){
		// upper bound
		Utils::DotProduct(&AiZ[(idx-1)*nz],x,nz,*err);
		*err -= ubineq[idx-1];
	}else{ 
		// lower bound
		Utils::DotProduct(&AiZ[(-idx-1)*nz],x,nz,*err);
		*err = lbineq[-idx-1]-*err;
	};

}