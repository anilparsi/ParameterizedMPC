function MOAS_eta(skip_iter,MAXITER,datafile)
%% To find maximum output admissible set
%% User defined settings

% skip_iter;    % The algorithm includes all constraints until this iteration
                % Can be used to speed up the algorithm

% MAXITER;      % If MAXITER is NaN, skip_iter will be set to t_star

%% Gilbert algorithm 3.2

CX = Cxu(:,1:n);
CU = Cxu(:,n+1:end);
% a constraint is non-redundant if either lb or ub is unsatisfied

converged = zeros(length(b_l),1);

options = optimoptions('linprog','display','off');

% Eliminate dynamics to improve speed 
Meq = -inv(kron(A,eye(s))-kron(eye(n),Md'))*kron(B,eye(s));   % eta_x = Meq * eta_u

% Warm start for t*: initially not many redundant constraints
j = 1;
cons = kron(CX(px+1:end,:),tau0d')*Meq + kron(CU(px+1:end,:),tau0d');     
                    % state constraint at t=0 is not needed
cons_lb = b_l(px+1:end);
cons_ub = b_u(px+1:end);
time_indices = [-1 1];      
tauk = [tau0d Md*tau0d];
for i = 1:j
    cons = [cons;         
            kron(CX,tauk(:,end)')*Meq + kron(CU,tauk(:,end)')];
    cons_lb = [cons_lb; b_l];   
    cons_ub = [cons_ub; b_u];
    time_indices(i+1,:) = ones(1,p); 
    tauk = [tauk Md*tauk(:,end)];
end


while ~all(converged) && j<MAXITER
    
    f1 = kron(CX,tauk(:,end)')*Meq + kron(CU,tauk(:,end)');            
    converged = run_cplex(f1,cons,cons_lb,cons_ub,b_l,b_u,options);
    
    if all(converged)
        t_star = j;
        break
    else
        cons = [cons; kron(CX(converged==0,:),tauk(:,end)')*Meq ...
                        + kron(CU(converged==0,:),tauk(:,end)')];   
                    
        cons_lb = [cons_lb; b_l(converged==0,:)];   
        cons_ub = [cons_ub; b_u(converged==0,:)];
        
        j = j+1;       
        
        
        for i = 1:p
            if ~converged(i)
                time_indices(j+1,i) = 1;
            else
                % store -1 for all those constraints which are redundant
                time_indices(j+1,i) = -1;
            end
            
        end
    end
    
    tauk = [tauk Md*tauk(:,end)];
end

%% change data representation into [eta_x;eta_u] = eta_z

% equality constraints (IC, dynamics)
Aeq = [kron(eye(n),tau0d') zeros(n,m*s);
       kron(A,eye(s))-kron(eye(n),Md'),kron(B,eye(s))];
D = [eye(n); zeros(n*s,n)];

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
    
    Aineq  = [Aineq; [kron(CX(ind,:),(Md^(i-1)*tau0d)'),    kron(CU(ind,:),(Md^(i-1)*tau0d)')]*Z];
    lbineq = [lbineq; b_l(ind)];
    ubineq = [ubineq; b_u(ind)];
end

% eta_z = C*x0 + Z*z
C = Y/(Req')*D;

AiC = Aineq*C;
AiZ = Aineq*Z;

%% Data for cons_skip
norms = zeros(t_star+1,1);

for i = 1:t_star+1
    norms(i) = norm(tauk(:,i)'* (Md-eye(s)));
end

Cs = kron(Cxu,eye(s));  % eta_w = Cs * eta_z

%% save data
save(datafile,'alpha','Ts','s','t_star','Cxu','b_l','b_u','Aineq','norms',...
      'Cs','lbineq','ubineq','time_indices','Aeq','tauk','AiC','AiZ','Z','C')
end