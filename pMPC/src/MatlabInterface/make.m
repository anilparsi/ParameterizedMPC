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
    
    % Uncomment the following flags for DEBUG mode
%     DEBUGFLAGS=  ' -g CXXDEBUGFLAGS=''$CXXDEBUGFLAGS -Wall -pedantic -Wshadow'' ';    
%     OPTIMIZEFLAGS= '';
    
    % Uncomment the following flags for OPTIMIZED code
    DEBUGFLAGS= '';
	OPTIMIZEFLAGS= 'OPTIMFLAGS="/Ox" ';

    cmd = [ 'mex -output pMPC ',IFLAGS, DEBUGFLAGS, OPTIMIZEFLAGS, '-D__SINGLE_OBJECT__ ', '-D__MATLAB__ ','pMPC_matlab.cpp'];
    
    eval(cmd);
    
    path( path,pwd );

end