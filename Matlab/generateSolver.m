function handle = generateSolver(moas)
% This function generates all the matrices required for the C++ implementation 
% of the problem and initializes the parameterized MPC solver.
%% User defined settings for parameterized MPC

% path to the .txt files: inputs for C++
path = fullfile('..','pMPC','build','Data','MPCmat');

tolMin  = 1e-9;     % minimum tolerance for constraint check
tolMax  = 1e-5;     % maximum tolerance for constraint check
MAXITER = 50;       % maximum number of iterations of active-set method

%% Derived variables

% Objective function
H = blkdiag(kron(moas.sys.Q,eye(moas.apx.s)), kron(moas.sys.R,eye(moas.apx.s)));

Z = moas.Z;

G  = moas.Z'*H*moas.Z;
F  = moas.Z'*H*moas.C;
g  = F*moas.sys.x0;
L  = chol(G,'lower');
Li = inv(L);

AiC = moas.Aineq * moas.C;
AiZ = moas.Aineq * moas.Z;


%% Generate matrices

% combine all the data in the parameter vector
params=[tolMin;                     
        tolMax;
        MAXITER;
        moas.sys.m * moas.apx.s - moas.sys.n;   % nz := number of variables in reduced problem
        length(moas.lbineq);                    % nc := number of inequality constraints
        moas.sys.n;
        moas.sys.m;
        moas.apx.s;
        ];
    
% print the parameters and problem data to txt files that can the be read
% by the c++ implementation
mkdir(path);

% (:) converts matrices into a vector, coloumn wise!
% Using transpose to correct this and store row-wise.
tmp=params';
vec2dense(tmp(:),fullfile(path,'params'));


tmp=AiC';
vec2dense(tmp(:),fullfile(path,'AiC'));

tmp=moas.C';
vec2dense(tmp(:),fullfile(path,'C'));

tmp=kron(eye(moas.sys.m),moas.apx.tau0d(:,1)')';
vec2dense(tmp(:),fullfile(path,'eta2u'));

tmp=F';
vec2dense(tmp(:),fullfile(path,'F'));

tmp=g';
vec2dense(tmp(:),fullfile(path,'g'));

tmp=moas.Z';
vec2dense(tmp(:),fullfile(path,'Z'));

tmp=Li';
vec2dense(tmp(:),fullfile(path,'Li'));

tmp = AiZ';
vec2dense(tmp(:),fullfile(path,'AiZ'));

tmp = (moas.Cs*moas.C)';
vec2dense(tmp(:),fullfile(path,'C0'));

tmp = (moas.Cs*Z)';
vec2dense(tmp(:),fullfile(path,'C1'));

tmp = int32(moas.time_indices');
vec2dense(tmp(:),fullfile(path,'time_indices'));

tmp = moas.norms; 
vec2dense(tmp(:),fullfile(path,'norms'));

% store tauk column wise, because each column has tau(k);
tmp = moas.tauk; 
vec2dense(tmp(:),fullfile(path,'tauk'));

tmp = moas.sys.b_l';
vec2dense(tmp(:),fullfile(path,'b_l'));

tmp = moas.sys.b_u';
vec2dense(tmp(:),fullfile(path,'b_u'));

tmp = moas.lbineq';
vec2dense(tmp(:),fullfile(path,'lbineq'));

tmp = moas.ubineq';
vec2dense(tmp(:),fullfile(path,'ubineq'));

tmp=moas.sys.A';
vec2dense(tmp(:),fullfile(path,'A'));

tmp=moas.sys.B';
vec2dense(tmp(:),fullfile(path,'B'));

%% Initialize parameterized MPC solver

% generate mex function
run(fullfile('..','pMPC','src','MatlabInterface','make.m'))
% initialize solver
handle = pMPC('i',path);
