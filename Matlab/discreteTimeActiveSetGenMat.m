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
function discreteTimeActiveSetGenMat(A,B,Q_cost,R_cost,Md,tau0d,x0,datafile)
    
% need to make sure that the mex-compiled-function pdMPC is in the path.
addpath('..\..\..\pdMPC_discrete\src\MatlabInterface');
%% NEW MATRICES
s = length(tau0d);  % number of basis functions
n = length(A);      % number of states 
m = length(B(1,:)); % number of inputs

% Objective function
H = blkdiag(kron(Q_cost,eye(s)), kron(R_cost,eye(s)));
h = zeros(n*s+m*s,1);

% load inequality constraints from MOAS
meas = load(datafile);

Z = meas.Z;
C = meas.C; 

G  = Z'*H*Z;
F  = Z'*H*C;
g  = F*x0 + Z'*h;  
L  = chol(G,'lower');
Li = inv(L);

AiC = meas.AiC;
AiZ = meas.AiZ;

%%

% set up tolerances 
tolMin  = 1e-9;     % minimum tolerance for constraint check
tolMax  = 1e-5;     % maximum tolerance for constraint check
MAXITER = 50;      % maximum number of iterations

%%

% combine all the data in the parameter vector
params=[tolMin;
        tolMax;
        MAXITER;
        m*s-n;
        length(meas.lbineq);
        n;
        m;
        s;
        ];
% print the parameters and problem data to txt files that can the be read
% by the c++ implementation
path = '..\..\..\pdMPC_discrete\build\Data\MPCmat';
mkdir(path);


% (:) converts matrices into a vector, coloumn wise!
% Using transpose to correct this and store row-wise.
tmp=params';
vec2dense(tmp(:),strcat(path,'\params'));

tmp=g';
vec2dense(tmp(:),strcat(path,'\g'));

tmp=AiC';
vec2dense(tmp(:),strcat(path,'\AiC'));

tmp=C';
vec2dense(tmp(:),strcat(path,'\C'));

tmp=kron(eye(m),tau0d(:,1)')';
vec2dense(tmp(:),strcat(path,'\eta2u'));

tmp=F';
vec2dense(tmp(:),strcat(path,'\F'));

tmp=Z';
vec2dense(tmp(:),strcat(path,'\Z'));

tmp=Li';
vec2dense(tmp(:),strcat(path,'\Li'));

tmp = (Li'*Li)';
vec2dense(tmp(:),strcat(path,'\LiTLi'));

tmp=AiZ';
vec2dense(tmp(:),strcat(path,'\AiZ'));

tmp = (meas.Cs*C)';
vec2dense(tmp(:),strcat(path,'\C0'));

tmp = (meas.Cs*Z)';
vec2dense(tmp(:),strcat(path,'\C1'));

tmp = int32(meas.time_indices');
vec2dense(tmp(:),strcat(path,'\time_indices'));

tmp = meas.norms; 
vec2dense(tmp(:),strcat(path,'\norms'));

% store tauk column wise, because each column has tau(k);
tmp = meas.tauk; 
vec2dense(tmp(:),strcat(path,'\tauk'));

tmp = meas.b_l';
vec2dense(tmp(:),strcat(path,'\b_l'));

tmp = meas.b_u';
vec2dense(tmp(:),strcat(path,'\b_u'));

tmp = meas.lbineq';
vec2dense(tmp(:),strcat(path,'\lbineq'));

tmp = meas.ubineq';
vec2dense(tmp(:),strcat(path,'\ubineq'));

tmp=A';
vec2dense(tmp(:),strcat(path,'\A'));

tmp=B';
vec2dense(tmp(:),strcat(path,'\B'));

