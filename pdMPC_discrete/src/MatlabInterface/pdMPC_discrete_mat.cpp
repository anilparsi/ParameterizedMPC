#pragma once

#include "mex.h"
#include "QPSolver.cpp"
#include "MPCSolver.cpp"
#include "ActiveConstraints.cpp"
#include "Rmatrix.cpp"
#include "Utils.cpp"
#include <string>
#include <vector>

/* Work-around for settings where mexErrMsgTxt causes unexpected behaviour. */
#ifdef __AVOID_MEXERRMSGTXT__
	#define myMexErrMsgTxt( TEXT ) mexPrintf( "%s\n\n",(TEXT) );
#else
	#define myMexErrMsgTxt mexErrMsgTxt
#endif



// global pointer to ParametrizedMPC objects
static std::vector<MPCSolver*> MPC_instances;


int_t allocateMPCProblem(std::string dir){
	MPC_instances.push_back(new MPCSolver(dir));
	return (int_t)(MPC_instances.size()-1);
}

void deleteMPCProblem(int_t handle){
	delete MPC_instances.at(handle);
}


void mexFunction( int_t nlhs, mxArray* plhs[], int_t nrhs, const mxArray* prhs[]){
	if(nrhs>3){
		myMexErrMsgTxt( "ERROR (parametrizedMPC): Invalid number of input arguments!\n" );
		return;
	}

	if (mxIsChar( prhs[0] ) != 1){
		myMexErrMsgTxt( "ERROR (parametrizedMPC): First input argument must be a string!\n" );
		return;
	}
	char typeString[2];

	mxGetString( prhs[0], typeString, 2);

	if(strcmp( typeString,"i")==0){
		// initialize
		if(mxIsChar( prhs[1] ) != 1){
			myMexErrMsgTxt( "ERROR (parametrizedMPC): Second input argument must be a string when initializing!\n" );
		}
		char dirString[100];
		mxGetString( prhs[1],dirString,100);

		int_t handle=allocateMPCProblem(std::string(dirString));
		nlhs=1;
		plhs[0]=mxCreateDoubleScalar(handle);
		
	}else if(strcmp( typeString,"s")==0){
		// get instance
		int_t handle = (uint_t)mxGetScalar( prhs[1] );

		// get x0
		double* x0_vec=(double*)mxGetPr( prhs[2] );
		
		
		// solve mpc problem
		if(MPC_instances.size()>handle){
			MPC_instances[handle]->solve(x0_vec);
			/*
			if(ret!=1){
				// do real_t max_error handling
				if(ret==-1){
					myMexErrMsgTxt( "ERROR (parametrizedMPC): Problem is infeasible!\n");
				}else if(ret==-2){
					myMexErrMsgTxt( "ERROR (parametrizedMPC): Unknown problem!\n");
				}
			}
			*/
		}else{
			// MPC instance has not been allocated
			myMexErrMsgTxt( "ERROR (parametrizedMPC): MPC problem has not yet been allocated!\n" );
			return;		
		}

		// retrieve solution
		nlhs=3;
		plhs[0]=mxCreateDoubleScalar(handle);
		plhs[1]= mxCreateDoubleMatrix(MPC_instances[handle]->getNumberOfOutputs(),1,mxREAL);
		plhs[2]=mxCreateDoubleScalar((real_t)MPC_instances[handle]->getIterNumber());
		double *u= mxGetPr( plhs[1] );
		MPC_instances[handle]->getControlInputs(u);
		//printf("u is %f\n",*u);

	}
	else if(strcmp( typeString,"d")==0){
		int_t handle = (uint_t)mxGetScalar( prhs[1] );
		deleteMPCProblem(handle);
		nlhs=1;
		plhs[0]=mxCreateDoubleScalar(-1);
	}
	else{
		myMexErrMsgTxt( "ERROR (parametrizedMPC): Command string is not correct!\n" );
	}
}
