#include "mex.h"
#include "../QPSolver.cpp"
#include "../MPCSolver.cpp"
#include "../ActiveConstraints.cpp"
#include "../Rmatrix.cpp"
#include "../Utils.cpp"
#include <string>
#include <vector>

/* Work-around for settings where mexErrMsgTxt causes unexpected behaviour. */
#ifdef __AVOID_MEXERRMSGTXT__
    #define myMexErrMsgIdAndTxt( ID, TEXT ) printf( "%s\n",(TEXT) );
#else
    #define myMexErrMsgIdAndTxt mexErrMsgIdAndTxt
#endif


// global pointer to ParametrizedMPC objects
static std::vector<MPCSolver*> MPC_instances;


int_t allocateMPCProblem(std::string dir){
    MPC_instances.push_back(new MPCSolver(dir));
    printf("Initialized handle %i.\n",(MPC_instances.size()-1));
    return (int_t)(MPC_instances.size()-1);     
}

void deleteMPCProblem(int_t handle){
    if(MPC_instances.at(handle)){
    	delete MPC_instances.at(handle);  
    	printf("Deleted handle %i.\n",handle);
    	MPC_instances[handle] = NULL;  	
	}
	else{
		myMexErrMsgIdAndTxt("pMPC:NoHandle", "This handle does not exist.");
	}

}


void mexFunction( int_t nlhs, mxArray* plhs[], int_t nrhs, const mxArray* prhs[]){
    if(nrhs>3){
        myMexErrMsgIdAndTxt( "pMPC:InvalidInputs", "ERROR (parametrizedMPC): Invalid number of input arguments!\n" );
        return;
    }

    if (mxIsChar( prhs[0] ) != 1){
        myMexErrMsgIdAndTxt( "pMPC:FirstInputType","ERROR (parametrizedMPC): First input argument must be a string!\n" );
        return;
    }
    char typeString[2];

    mxGetString( prhs[0], typeString, 2);

    if(strcmp( typeString,"i")==0){
        // initialize
        if(mxIsChar( prhs[1] ) != 1){
            myMexErrMsgIdAndTxt( "pMPC:SecondInputType","ERROR (parametrizedMPC): Second input argument must be a string when initializing!\n" );
        }
        char dirString[100];
        mxGetString( prhs[1],dirString,100);

        int_t handle = allocateMPCProblem(std::string(dirString));
        nlhs = 1;
        plhs[0] = mxCreateDoubleScalar(handle);
        
    }else if(strcmp( typeString,"s")==0){
        // get instance
        int_t handle = (uint_t)mxGetScalar( prhs[1] );
        
        // check handle 
        if(!MPC_instances.at(handle)){
    		myMexErrMsgIdAndTxt("pMPC:NoHandle", "This handle does not exist.");
            return;     
    	}    

        // get x0
        double* x0_vec = (double*)mxGetPr( prhs[2] );
                
        // solve mpc problem        
        MPC_instances[handle]->solve(x0_vec);	

        // retrieve solution
        nlhs = 3;
        plhs[0] = mxCreateDoubleScalar(handle);
        plhs[1] = mxCreateDoubleMatrix(MPC_instances[handle]->getNumberOfOutputs(),1,mxREAL);
        plhs[2] = mxCreateDoubleScalar((real_t)MPC_instances[handle]->getIterNumber());
        double *u = mxGetPr( plhs[1] );
        MPC_instances[handle]->getControlInputs(u);     
    }
    else if(strcmp( typeString,"d")==0){
        int_t handle = (uint_t)mxGetScalar( prhs[1] );
        deleteMPCProblem(handle);
        nlhs = 1;
        plhs[0] = mxCreateDoubleScalar(-1);
    }
    else{
        myMexErrMsgIdAndTxt("pMPC:InvalidString", "ERROR (parametrizedMPC): Command string is not correct!\n" );
    }
}
