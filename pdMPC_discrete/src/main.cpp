#include "QPSolver.h"
#include "MPCSolver.h"
#include "DefineSettings.h"
#include "Utils.h"
#include <cmath>
#include <string>


// Check implementation of parameterized MPC 

int main(){
	real_t *A,*B;
	int_t n,m,temp;
	std::string dir = "Data/MPCmat";

	// load matrices A and B
	std::string tmp=dir+"/A";
	Utils::LoadVec(tmp.c_str(),&A,temp);
	n = (int_t)sqrt(temp);

	tmp=dir+"/B";
	Utils::LoadVec(tmp.c_str(),&B,temp);
	m = temp/n;
	real_t *Ax = new real_t[n],*Bu = new real_t [n];

	// loads MPC problem data into pmpc object
	MPCSolver pmpc(dir);

	int_t tmax = 1000;						// simulation duration
	real_t x0[] = {0.5, 0.5, 0.5, 0.0};
	
	
	real_t *x = new real_t [tmax*n],
		   *u = new real_t [tmax*m];
	
	Utils::VectorCopy(x0,x,n);

	// for the QuadMPC, solve one time step outside loop	(due to u(-1))
	pmpc.solve(&x[0]);	
	pmpc.getControlInputs(&u[0]);
	Utils::MatrixMult(A,&x[0],Ax,n,n,1);
	Utils::MatrixMult(B,&u[0],Bu,n,m,1);
	Utils::VectorAdd(Ax,Bu,&x[n],n);
	// simulate for tmax;
	for(int i=1; i<tmax-1;++i){	
		pmpc.solve(&x[n*i]);	
		pmpc.getControlInputs(&u[m*i]);
		Utils::MatrixMult(A,&x[n*i],Ax,n,n,1);
		Utils::MatrixMult(B,&u[m*i],Bu,n,m,1);
		Utils::VectorAdd(Ax,Bu,&x[n*(i+1)],n);
	// printf("state is %f, %f, %f, %f.\n",x[n*i-4],x[n*i-3],x[n*i-2],x[n*i-1]);	
	// printf("input is %f.\n",u[0]);	
	}


	delete[] x;
	delete[] u;
	delete[] Ax;
	delete[] Bu;
	return 0;
}


/*
// Check the QPSolver constructor with matrices given as pointers
int main(){

real_t *Li,*g,*Aineq,*lbineq,*ubineq;
int_t nz,nc,temp;
std::string dir = "Data/MPCmat";

// load matrices
std::string tmp = dir+"/g";
Utils::LoadVec(tmp.c_str(), &g, nz);

tmp = dir + "/Li";
Utils::LoadVec(tmp.c_str(), &Li, temp);

tmp = dir + "/AiZ";
Utils::LoadVec(tmp.c_str(), &Aineq, temp);

tmp = dir + "/lbineq";
Utils::LoadVec(tmp.c_str(), &lbineq, nc);

tmp = dir + "/ubineq";
Utils::LoadVec(tmp.c_str(), &ubineq, temp);

// Call constructor from matrices
QPSolver solver(Li,g,Aineq,lbineq,ubineq,nz,nc,1e-9,1e-9);

solver.solve();

real_t *z = new real_t[nz];
solver.getSolutionCopy(z);

for (int_t i = 0; i < nz; i++)
{
	printf("\n%f", z[i]);
}

delete[] z;
return 0;
}
*/