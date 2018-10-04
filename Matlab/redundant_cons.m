function [cons_n,c_lb_n,c_ub_n,redundant,count,time_indices,cons_indices] = redundant_cons(cons,cons_lb,cons_ub,time_indices,cons_indices);
%% Find redundant constraints
% clhs are the LHS of constraints, crhs the RHS
% bz is used to set the first constraint set to not redundant
% tic
redundant = zeros(length(cons),2);
options = optimoptions('linprog','display','off');

% logical indexing for selecting non-redundant constraints
time_log = time_indices;
cons_log = cons_indices;
l1 = cellfun('length',time_indices)';
ind1 = 2; 
ind2 = 1;
time_log{1} = ones(1,l1(1));
for i = 1:length(cons_log)    
    cons_log{i}(1) = 1;
end

for i = l1(1)+1:length(cons_lb)
    cons_t = cons;
    c_lb_t = cons_lb;
    c_ub_t = cons_ub;

    cons_t(i,:) = [];
    c_lb_t(i) = [];
    c_ub_t(i) = [];
    
    f1 = cons(i,:);        
    
    [red,marg] = run_cplex(f1,cons_t,c_lb_t,c_ub_t,cons_lb(i),cons_ub(i),options);
    
    if red
        % constraint is redundant
        redundant(i,:) = [1,marg];        
        cons_log{time_indices{ind1}(ind2)}(ind1) = 0;
        time_log{ind1}(ind2) = 0;
    else
        % constraint is non-redundant
        redundant(i,:) = [0,marg];
        cons_log{time_indices{ind1}(ind2)}(ind1) = 1;
        time_log{ind1}(ind2) = 1;
    end
    
    ind2 = ind2 + 1;
    
    if ind2>l1(ind1)
        ind1 = ind1+1;
        ind2 = 1;
    end
end

% update indices
for i = 1:length(time_indices)
    time_indices{i} = time_indices{i}(time_log{i}==1);
end

for i = 1:length(cons_indices)
    cons_indices{i} = cons_indices{i}(cons_log{i}==1);
end

% output non-redundant constraints
cons_n = cons;
c_lb_n = cons_lb;
c_ub_n = cons_ub;

cons_n(redundant(:,1)==1,:) = [];
c_lb_n(redundant(:,1)==1,:) = [];
c_ub_n(redundant(:,1)==1,:) = [];

count = sum(redundant(:,1));
% toc
end

%% Change representation of constraints and run cplex
function [red,marg] = run_cplex(f1,cons,cons_lb,cons_ub,lb,ub,options)

tol = 1e6*eps;

% bounds because the problem is unbounded sometimes
LB = -1e6*ones(length(f1),1); 
UB = -LB;

% convert constraints into Ax<=b form
clhs = [cons;-cons];
crhs = [cons_ub;-cons_lb];


% LP for lb (minimisation problem)
[~,fval1,exitflag1,~] = cplexlp(f1',clhs,crhs,[],[],LB,UB,[],options);

% LP for ub (minimisation problem)    
[~,fval2,exitflag2,~] = cplexlp(-f1',clhs,crhs,[],[],LB,UB,[],options);

marg = max((-fval1 +lb),(-fval2-ub));

if exitflag1<0 || exitflag2<0
    keyboard
elseif marg<=-tol
    red = 1;
else
    red = 0;
end

end