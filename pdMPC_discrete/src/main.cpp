#pragma once
#include "QPSolver.h"
#include "MPCSolver.h"
#include "DefineSettings.h"
#include "Utils.h"

#include <string>

#ifdef TIME
	extern real_t tqpTot;
	extern real_t tqpStart;
	extern real_t tupTot;
	extern real_t tupStart;
	extern real_t tdwnTot;
	extern real_t tdwnStart;
	extern real_t tdwnTot;
	extern real_t tdwnStart;
	extern real_t tchkTot;
	extern real_t tchkStart;
	extern real_t titTot;
	extern real_t titStart;
	extern int_t it_count;
#endif

int main(){
	
#ifdef TIME
	real_t tstart,tqpStart,tqpTot=0.0;
	tstart=Utils::getCPUtime();
#endif
	real_t *A,*B;
	int_t n,m,temp;
	std::string dir = "Data/MPCmat";

	// load matrices A and B
	std::string tmp=dir+"/A";
	Utils::LoadVec(tmp.c_str(),&A,temp);
	n = (int_t)Utils::getSqrt(temp);

	tmp=dir+"/B";
	Utils::LoadVec(tmp.c_str(),&B,temp);
	m = temp/n;
	real_t *Ax = new real_t[n],*Bu = new real_t [n];

	// loads MPC problem data into pmpc object
	MPCSolver pmpc(dir);

	int_t tmax = 600;						// simulation duration
	real_t x0[] = {0.5, 0.5, 0.5, 0.5};
	
	
	real_t *x = new real_t [tmax*n],
		   *u = new real_t [tmax*m];
	
	Utils::VectorCopy(x0,x,n);

	// for the QuadMPC, solve one time step outside loop	(due to u(-1))
	pmpc.solve(&x[0]);	
	pmpc.getSolutionCopy(&u[0]);
	Utils::MatrixMult(A,&x[0],Ax,n,n,1);
	Utils::MatrixMult(B,&u[0],Bu,n,m,1);
	Utils::VectorAdd(Ax,Bu,&x[n],n);

	// simulate for tmax;
	for(int i=1; i<tmax-1;++i){	
#ifdef TIME
	tqpStart = Utils::getCPUtime();
#endif
		pmpc.solve(&x[n*i]);	
		pmpc.getControlInputs(&u[m*i]);
#ifdef TIME
	tqpTot += Utils::getCPUtime()-tqpStart;
#endif
		Utils::MatrixMult(A,&x[n*i],Ax,n,n,1);
		Utils::MatrixMult(B,&u[m*i],Bu,n,m,1);
		Utils::VectorAdd(Ax,Bu,&x[n*(i+1)],n);
	}
	
#ifdef TIME
	tstart=Utils::getCPUtime()-tstart;
	printf("\nexec time: %.5f",tstart);
	printf("\nQP solve time: %.5f",tqpTot);
	printf("\niteration time: %.5f, iteration count %d",titTot,it_count);
	
	printf("\n\ncons check time: %.5f",tchkTot);
#endif


	delete[] x;
	delete[] u;
	delete[] Ax;
	delete[] Bu;
	return 0;
}
