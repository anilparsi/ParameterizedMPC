function apx = defineApproximation(sys)
%% Function to define approximations
% 
% Default basis functions: Laguerre functions
% alpha : decay constant
% s     : number of  basis functions
% Md    : basis function evolution matrix
% tau0d : discrete-time basis functions at k = 0
% tau0c : continuous-time basis functions at k = 0

%% User-defined variables
% apx.alpha = getDefaultAlpha(sys.A,sys.B,sys.Q,sys.R);         
apx.alpha = 0.7;               

% apx.s     = floor(sys.n/sys.m)*2;                  % min s: n/m    
apx.s     = 9;    

%% Derived variables
if apx.s <= floor(sys.n/sys.m)
    error('Increase number of basis functions')
end

% Laguerre basis functions(default):
tau0c = sqrt(2*apx.alpha)*ones(apx.s,1);
Mc    = -2*apx.alpha*tril(ones(apx.s))+apx.alpha*eye(apx.s);

% obtain the discretized version by taking the matrix exponential
apx.Md = expm(Mc*sys.Ts); 

% orthonormalize to obtain discrete-time basis functions
P     = dlyap(apx.Md,tau0c*tau0c');
T     = chol(P,'lower');       % transformation matrix
apx.Md    = T\(apx.Md*T);
apx.tau0d = T\tau0c;

end