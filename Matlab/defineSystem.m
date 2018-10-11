function sys = defineSystem()
%% Function to define system matrices
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
% stack all state constraints at the top
sys.Cxu = [1 0 0 0 0
           0 0 0 0 1];   
sys.b_l = [-.1; -0.5];
sys.b_u = [10; 0.5,];


% cost matrices: cost(k) =  x(k)^T Q x(k) + u(k)^T R u(k)
sys.Q = diag([1 1 1 1]);
sys.R = 0.05;

% sampling time for discretization
sys.Ts=.02;

% Give initial condition
sys.x0=5*[.1;.1;.1;.1];

%% Derived variables

sys.n = size(Ac,1); 
sys.m = size(Bc,2); 
sys.p = size(sys.Cxu,1);

% find number of state constraints
for i = 1:sys.n
    if any(sys.Cxu(i,sys.n+1:end))
        sys.px = i-1;
        break; 
    end    
end

% define discrete-time system
sysd=c2d(ss(Ac,Bc,eye(sys.n),zeros(sys.n,sys.m)),sys.Ts);
sys.A = sysd.a;
sys.B = sysd.b;

% normalize constraints
for i = 1:sys.p
    temp = norm(sys.Cxu(i,:));
    sys.Cxu(i,:) = sys.Cxu(i,:)/temp;
    sys.b_l(i) = sys.b_l(i)/temp;
    sys.b_u(i) = sys.b_u(i)/temp;
end

end