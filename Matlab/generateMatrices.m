%% Modified for QR implementation
% 24.05.18

%% 
% This function is intended to be a matlab interface to the C++
% implementation. It collects and calculates all the necessary problem data
% for the C++ implementation.
%
% Input:    H   Hessian of the QP.
%           A,B equality constraint of the QP (A*eta=B*x0, with x0 the
%               specific initial condition)
%           T   matrix that defines the semi-infinite input and state 
%               constraints. The constraint reads as u_l <= T [x;u] <= u_u, 
%               where x is the state trajectory, u the input trajectory,
%               that are both parametrized with the basis functions.
%       Md,tau0 defines the basis functions \tau(k+1)=Md \tau(k),
%               \tau(0)=tau0
%
% Output:  handle  handle to the generated qp-instance, obtain the solution
%                  via pdMPC('s',handle,x0)
%
function generateMatrices(moas,path)
%% User defined settings for parameterized MPC

tolMin  = 1e-9;     % minimum tolerance for constraint check
tolMax  = 1e-5;     % maximum tolerance for constraint check
MAXITER = 50;       % maximum number of iterations of active-set method

%% Derived variables

% add mex-compiled-function to the path.
addpath('..\pdMPC_discrete\src\MatlabInterface');

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
vec2dense(tmp(:),strcat(path,'\params'));


tmp=AiC';
vec2dense(tmp(:),strcat(path,'\AiC'));

tmp=moas.C';
vec2dense(tmp(:),strcat(path,'\C'));

tmp=kron(eye(moas.sys.m),moas.apx.tau0d(:,1)')';
vec2dense(tmp(:),strcat(path,'\eta2u'));

tmp=F';
vec2dense(tmp(:),strcat(path,'\F'));

tmp=g';
vec2dense(tmp(:),strcat(path,'\g'));

tmp=moas.Z';
vec2dense(tmp(:),strcat(path,'\Z'));

tmp=Li';
vec2dense(tmp(:),strcat(path,'\Li'));

tmp = (Li'*Li)';
vec2dense(tmp(:),strcat(path,'\LiTLi'));

tmp = AiZ';
vec2dense(tmp(:),strcat(path,'\AiZ'));

tmp = (moas.Cs*moas.C)';
vec2dense(tmp(:),strcat(path,'\C0'));

tmp = (moas.Cs*Z)';
vec2dense(tmp(:),strcat(path,'\C1'));

tmp = int32(moas.time_indices');
vec2dense(tmp(:),strcat(path,'\time_indices'));

tmp = moas.norms; 
vec2dense(tmp(:),strcat(path,'\norms'));

% store tauk column wise, because each column has tau(k);
tmp = moas.tauk; 
vec2dense(tmp(:),strcat(path,'\tauk'));

tmp = moas.sys.b_l';
vec2dense(tmp(:),strcat(path,'\b_l'));

tmp = moas.sys.b_u';
vec2dense(tmp(:),strcat(path,'\b_u'));

tmp = moas.lbineq';
vec2dense(tmp(:),strcat(path,'\lbineq'));

tmp = moas.ubineq';
vec2dense(tmp(:),strcat(path,'\ubineq'));

tmp=moas.sys.A';
vec2dense(tmp(:),strcat(path,'\A'));

tmp=moas.sys.B';
vec2dense(tmp(:),strcat(path,'\B'));

