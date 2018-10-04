%% Script to define system matrices
% Ac    : continuous-time state matrix
% Bc    : continuous-time input to state matrix
% A     : discrete-time state matrix
% B     : discrete-time input to state matrix
% Cxu   : constraint matrix on the state and input variables
% b_l   : lower bounds on constraints
% b_u   : upper bounds on constraints
% Q     : penalty on state variables
% R     : penalty on input variables
% n     : number of states
% m     : number of inputs
% p     : number of constraints at each time step

%% User-defined variables
% continuous-time description \dot{x}=Ax + Bu
Ac = [0,1,0,0;
      0,0,1,0;
      0,0,0,1;
      0,0,0,0];
Bc=[0;0;0;1;];

% constraints: b_l <= Cxu * [x(k);u(k)] <= b_u
Cxu = [1 0 0 0 0
       0 0 0 0 1];   
b_l = [-.1; -0.5];
b_u = [10; 0.5,];

% cost matrices: cost(k) =  x(k)^T Q x(k) + u(k)^T R u(k)
Q = diag([1 1 1 1]);
R = 0.05;

% sampling time for discretization
Ts=.02;

% Give initial condition
x0=5*[.1;.1;.1;.1];

%% Derived variables

n = size(Ac,1); 
m = size(Bc,2); 
p = size(Cxu,1);

% define discrete-time system
sysd=c2d(ss(Ac,Bc,eye(n),zeros(n,m)),Ts);
A = sysd.a;
B = sysd.b;

% normalize constraints
for i = 1:p
    temp = norm(Cxu(i,:));
    Cxu(i,:) = Cxu(i,:)/temp;
    b_l(i) = b_l(i)/temp;
    b_u(i) = b_u(i)/temp;
end
