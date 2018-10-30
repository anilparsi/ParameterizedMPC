function moas = generateMOAS(sys,apx)
%% To find maximum output admissible set(MOAS)
% The MOAS for an autonomous system can be given by t_star, which is the
% time index until which the constraints must be included. A variable
% time_indices can be used to give the time-based active constraints at
% each time step. 
% This function uses a branch and bound technique to find the value of
% t_star. The window in which t_star is estimated is given in the value
% window. 
%% User defined settings

% Window size to estimate t_star 
window = 50;

% Maximum value of t_star the algorithm can converge to
max_size = 1600;  

% Data storage variables
save_data = true;
datafile = 'moas.mat';

%% Algorithm to generate MOAS

CX = sys.Cxu(:,1:sys.n);
CU = sys.Cxu(:,sys.n+1:end);

options = optimoptions('linprog','display','off');

% Eliminate dynamics to improve speed        eta_x = Meq * eta_u
Meq = -inv(kron(sys.A,eye(apx.s))- ...
       kron(eye(sys.n),apx.Md'))*kron(sys.B,eye(apx.s));   

% Generate all constraints until max_size
tauk = [apx.tau0d];
all_cons = [kron(CX,tauk(:,end)')*Meq + kron(CU,tauk(:,end)')];
for k = 1:max_size
    tauk = [tauk apx.Md*tauk(:,end)];
    all_cons = [all_cons;         
            kron(CX,tauk(:,end)')*Meq + kron(CU,tauk(:,end)')];
end
all_cons_lb = repmat(sys.b_l,max_size+1,1);
all_cons_ub = repmat(sys.b_u,max_size+1,1);

% time_indices: contains the active constraints at each time step
% time_indices(i,j) =  1 if j-th constraint at i-th time step is active, 
%                   = -1 if j-th constraint at i-th time step is inactive
% initialize all constraints as active
time_indices = ones(max_size+1,sys.p); 
% State constriants at initial time must be inactive
time_indices(1,:) = [-1*ones(1,sys.px) ones(1,sys.p-sys.px)];      

t_min = 1;
t_max = max_size;
% Check t_star = max_size
[cons, cons_lb, cons_ub] = generateConsSet(all_cons, all_cons_lb, all_cons_ub, time_indices);

f1 = kron(CX,tauk(:,end)')*Meq + kron(CU,tauk(:,end)');      
[converged,exitflag] = solve_lp(f1,cons,cons_lb,cons_ub,sys.b_l,sys.b_u,options);
loopflag = true;
if exitflag == -2 || exitflag == -4 
    % infeasible problem
    error('Infeasibility: Please change the approximation');
    return;
elseif exitflag == -3
    % unbounded problem
    warning('Unboundedness: Using t_star = max_size');
    t_star = max_size;
    loopflag = false;
elseif ~all(converged)
    % t_star > max_size for the given problem
    fprintf('No convergence. t_star = max_size being used. \n');
    t_star = max_size;
    loopflag = false;
end

while t_max-t_min > window && loopflag 
    % current estimate for t_star
    t_guess = round((t_max+t_min)/2);
    
    % update indices
    time_indices = update_indices(time_indices,t_max,t_min,t_guess,converged);
    
    % generate constraint set
    [cons, cons_lb, cons_ub] = generateConsSet(all_cons, all_cons_lb, all_cons_ub, time_indices);
    
    % objective function
    f1 = kron(CX,tauk(:,t_guess+1)')*Meq + kron(CU,tauk(:,t_guess+1)');            
    [converged, exitflag] = solve_lp(f1,cons,cons_lb,cons_ub,sys.b_l,sys.b_u,options);
    
    if exitflag == -2 
        % infeasible problem
        error('Infeasibility: Please change the approximation');
    elseif all(converged)
        % all constraints converged
        t_max = t_guess; 
    else
        % not all constraints converged, or problem is unbounded
        % increase t_min
        t_min = t_guess;
    end    
end

if t_max-t_min <= window 
   t_star = t_guess; 
end

%% change data representation into [eta_x;eta_u] = eta_z

% equality constraints (IC, dynamics)
Aeq = [kron(eye(sys.n),apx.tau0d') zeros(sys.n,sys.m*apx.s);
       kron(sys.A,eye(apx.s))-kron(eye(sys.n),apx.Md'),kron(sys.B,eye(apx.s))];
D = [eye(sys.n); zeros(sys.n*apx.s,sys.n)];

peq = length(Aeq(:,1));   % # of equality constraints
[Qeq,Req] = qr(Aeq');

% Change of variables: eta_z = Y*y + Z*z 
% Aeq*Z = 0
Y = Qeq(:,1:peq);
Z = Qeq(:,peq+1:end);
Req = Req(1:peq,1:peq);

% inequality constraints
Aineq = [];
lbineq = [];
ubineq = [];
for i = 1:length(time_indices)    
    ind = time_indices(i,:)>0;    
    Aineq  = [Aineq; [kron(CX(ind,:),(apx.Md^(i-1)*apx.tau0d)'),...
                kron(CU(ind,:),(apx.Md^(i-1)*apx.tau0d)')]];
    lbineq = [lbineq; sys.b_l(ind)];
    ubineq = [ubineq; sys.b_u(ind)];
end

% eta_z = C*x0 + Z*z
C = Y/(Req')*D;

AiC = Aineq*C;
AiZ = Aineq*Z;

%% Data for parameterized constraint check
norms = zeros(t_star+1,1);

for i = 1:t_star+1
    norms(i) = norm(tauk(:,i)'* (apx.Md-eye(apx.s)));
end

Cs = kron(sys.Cxu,eye(apx.s));  % eta_w = Cs * eta_z

%% save data
moas.sys = sys;
moas.apx = apx;
moas.t_star = t_star;
moas.Aineq = Aineq;
moas.norms = norms;
moas.Cs = Cs;
moas.lbineq = lbineq;
moas.ubineq = ubineq;
moas.time_indices = time_indices;
moas.tauk = tauk;
moas.Z = Z;
moas.C = C;

if (save_data)    
    save(datafile,'moas')
end


end

%% Function to generate constraint sets 
function [cons, cons_lb, cons_ub] = generateConsSet(all_cons, all_cons_lb, all_cons_ub, time_indices)
    % get active constraint indices
    indices = reshape(time_indices,[],1);

    % build constraint set
    cons = all_cons(indices>0,:);
    cons_lb = all_cons_lb(indices>0,:);
    cons_ub = all_cons_ub(indices>0,:);
end

%% function to update indices of active constraints
function time_indices = update_indices(time_indices,t_max,t_min,t_guess,converged)
    time_indices(t_min+1:t_guess,:) = ones(t_guess-t_min,size(time_indices,2));
    time_indices(t_guess+1:t_max,:) = zeros(t_max-t_guess,size(time_indices,2));
end