function make(vargin)
% Generate the matlab interface of parameterized MPC.
%
% Example:
% make(clean)   cleans the mex files in the directory
% make          creates a mex file from the C++ code

if nargin==1
    if strcmp(vargin,'clean')
        eval('delete *.o;');
        eval('delete *.mexw64*;');
        disp('successfully cleaned');
    end
else
    GFDGPATH= '..\';
    
    IFLAGS = [' -I',GFDGPATH,' '];
    
%     DEBUGFLAGS=  ' -g CXXDEBUGFLAGS=''$CXXDEBUGFLAGS -Wall -pedantic -Wshadow'' ';    
%     OPTIMIZEFLAGS= '';
    
    DEBUGFLAGS= '';
	OPTIMIZEFLAGS= 'OPTIMFLAGS="/Ox" ';

    cmd = [ 'mex -output pdMPC_discrete ',IFLAGS, DEBUGFLAGS, OPTIMIZEFLAGS, '-D__SINGLE_OBJECT__ ', '-D__MATLAB__ ','pdMPC_discrete_mat.cpp'];
    
    eval(cmd);
    
    path( path,pwd );

end