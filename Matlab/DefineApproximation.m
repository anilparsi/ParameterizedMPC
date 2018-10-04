%% Script to define approximations
% Default basis functions: Laguerre functions
% alpha : decay constant
% s     : number of  basis functions
% Md    : basis function evolution matrix
% tau0d : discrete-time basis functions at k = 0
% tau0c : continuous-time basis functions at k = 0

%% User-defined variables
alpha = getDefaultAlpha(A,B,Q,R);                
s     = floor(n/m)*2;                  % min s: n/m       

%% Derived variables
if s <= floor(n/m)
    error('Increase number of basis functions')
end

% Laguerre basis functions(default):
tau0c = sqrt(2*alpha)*ones(s,1);
Mc    = -2*alpha*tril(ones(s))+alpha*eye(s);

% obtain the discretized version by taking the matrix exponential
Md = expm(Mc*Ts); 

% orthonormalize to obtain discrete-time basis functions
T     = chol(dlyap(Md,tau0c*tau0c'),'lower');       % transformation matrix
Md    = T\(Md*T);
tau0d = T\tau0c;