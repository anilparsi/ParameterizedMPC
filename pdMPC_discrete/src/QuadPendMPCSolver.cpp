#include "QuadPendMPCSolver.h"
#include "QuadPendActiveCons.h"

QuadPendMPCSolver::QuadPendMPCSolver(std::string dir) : QuadMPCSolver(dir) {
	std::string tmp;
	int_t tmp2;
	
	{
		tmp = dir + "/params";
		real_t *tmpvec;

		Utils::LoadVec(tmp.c_str(), &tmpvec, tmp2);

		mu = tmpvec[8];

		delete[] tmpvec;
	}

	tmp = dir + "/soft_indices";
	Utils::LoadVec(tmp.c_str(), &soft_indices, tmp2);			// tmp2 = nc*n

	w_b		= new real_t[nz]();
	temp_wb = new real_t[nz]();

	// create activeCons with one extra variable (slack)
	delete activeCons;
	activeCons = new QuadPendActiveCons(AiZ, lbineq, ubineq, Li, nz+1,soft_indices,mu);
}

QuadPendMPCSolver::~QuadPendMPCSolver(){
	delete[] w_b;
	delete[] temp_wb;

	delete[] soft_indices;
}

void QuadPendMPCSolver::calc_z()
{
	QPSolver::calc_z();
	
	// calculate ep = -w_b'*lambda /mu;
	activeCons->dotProductFunction(lambda, ep);
	ep = -ep / mu;

	if (ep != ep) {
		// ep is nan
		int a = 1;
		// must set exitflag off
	}
}

void QuadPendMPCSolver::checkConstraints() {
	if (s < nz) {	// implement QPSolver::checkConstraints() with soft bounds
		viol_idx = 0;
		real_t max_error = -INFVAL;
		for (int i = 0; i<nc; ++i) {
			real_t prod = 0.0;

			for (int j = 0; j<nz; ++j) {
				prod += AiZ[i*nz + j] * z[j];
			}

			// errors
			real_t	e1 = prod - ubineq[i];
			// subtract slack from error
			if (ubineq[i] > 0)		{e1 -= ubineq[i] * soft_indices[i] * ep;}
			else					{e1 += ubineq[i] * soft_indices[i] * ep;}				
			
			if (e1>max_error)
			{	// update max error
				viol_idx = i + 1;
				max_error = e1;
			}

			e1 = lbineq[i] - prod;
			// subtract slack from error
			if (lbineq[i] > 0) { e1 -= lbineq[i] * soft_indices[i] * ep; }
			else { e1 += lbineq[i] * soft_indices[i] * ep; }

			if (e1>max_error)
			{	// update max error
				viol_idx = -i - 1;
				max_error = e1;
			}
		}
		viol = (max_error>TOL);
	}
	else 
	{
		checkConstraints_skip();
	}

}

void QuadPendMPCSolver::checkConstraints_skip() {
	//eta_w = C0*x0 + C1*z
	Utils::MatVecMult(C0, x0, temp_nw, m_nw, n);
	Utils::MatVecMult(C1, z, eta_w, m_nw, nz);
	Utils::VectorAdd(eta_w, temp_nw, eta_w, m_nw);

	for (int i = 0; i<m_np; ++i) {
		norm_w[i] = Utils::VectorNorm(&eta_w[i*s], s);
	}

	// val = tau0d'*eta_w;

#ifdef DEBUG_SKIP
	eval = 0;
#endif

	viol_idx = 0;
	real_t max_error = -INFVAL;
	// calculate the upper and lower bound errors for omega_rate
	// no slack for these
	for (int_t i = 0; i<m; ++i) {
		real_t prod = 0.0;

		for (int j = 0; j<nz; ++j) {
			prod += AiZ[i*nz + j] * z[j];
		}

		// errors
		real_t	e1 = prod - ubineq[i];

		if (e1>max_error)
		{	// update max error
			viol_idx = i + 1;
			max_error = e1;
		}

		e1 = lbineq[i] - prod;
		if (e1>max_error)
		{	// update max error
			viol_idx = -i - 1;
			max_error = e1;
		}
	}

	real_t val = 0.0;
	int_t idx1 =  m;							// index in constraint matrix

	// calculate the upper and lower bound errors for input constraints at t=0 
	for (int_t i = 0; i<m; ++i) {
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
	// calculate the upper and lower bound errors for state constraints at t=0 
	for (int_t i = m; i < m_np; ++i) {	// do not increment index because these are not in inequality matrix
		est_ubErr[i] = INFVAL;
		est_lbErr[i] = INFVAL;
	}

	// skip constraints loop
	for (int i = 0; i < t_star; ++i) {
		for (int k = 0; k<m_np; ++k) {

			if (time_indices[(i + 1)*m_np + k]>0) { // constraint is in non-redundant set
				++idx1;

				// update estimates
				est_ubErr[k] += norms[i] * norm_w[k];
				est_lbErr[k] += norms[i] * norm_w[k];

				if (est_ubErr[k] > max_error || est_lbErr[k] > max_error) {
					// estimate crosses bound: find exact value
					Utils::DotProduct(&tauk[(i + 1)*s], &eta_w[k*s], s, val);
					est_ubErr[k] = val - b_u[k];
					if (ubineq[idx1] > 0)	{ est_ubErr[k] -= ubineq[idx1] * soft_indices[idx1] * ep; }
					else					{ est_ubErr[k] += ubineq[idx1] * soft_indices[idx1] * ep; }

					est_lbErr[k] = b_l[k] - val;
					if (lbineq[idx1] > 0) { est_lbErr[k] -= lbineq[idx1] * soft_indices[idx1] * ep; }
					else { est_lbErr[k] += lbineq[idx1] * soft_indices[idx1] * ep; }

					if (est_ubErr[k]>max_error) {
						max_error = est_ubErr[k];
						viol_idx = idx1;
					}

					if (est_lbErr[k]>max_error) {
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

void QuadPendMPCSolver::updateMPCProblem(const real_t * const x_IC, const real_t * const u_old)
{
	QuadMPCSolver::updateMPCProblem(x_IC,u_old);
	activeCons->updateActiveSet();
}
