#include "QPSolver.h"
#include "MPCSolver.h"
#include "DefineSettings.h"
#include "Utils.h"
#include <cmath>
#include <string>


// Implementation of parameterized MPC for a given system

int main(){
	real_t *A,*B;
	int_t n,m,temp;
	std::string dir = "Data/MPCmat";

	// load matrices A and B
	std::string tmp = dir+"/A";
	Utils::LoadVec(tmp.c_str(),&A,temp);
	n = (int_t)sqrt(temp);

	tmp = dir+"/B";
	Utils::LoadVec(tmp.c_str(),&B,temp);
	m = temp/n;
	
	real_t *Ax = new real_t[n];
	real_t *Bu = new real_t [n];

	// define a parameterized MPC solver
	MPCSolver pmpc(dir);

	int_t tmax = 1000;						// simulation duration
	real_t x0[] = {0.3, 0.3, 0.3, 0.3};	
	
	real_t *x = new real_t [tmax*n],
		   *u = new real_t [tmax*m];
	
	Utils::VectorCopy(x0,x,n);

	// simulate for tmax;
	for(int i=0; i<tmax-1;++i){	
		pmpc.solve(&x[n*i]);	
		pmpc.getControlInputs(&u[m*i]);
		Utils::MatrixMult(A,&x[n*i],Ax,n,n,1);
		Utils::MatrixMult(B,&u[m*i],Bu,n,m,1);
		Utils::VectorAdd(Ax,Bu,&x[n*(i+1)],n);
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