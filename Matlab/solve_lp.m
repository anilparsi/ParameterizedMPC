%% Change representation of constraints and run cplex
function [converged, exitflag] = solve_lp(f,cons,cons_lb,cons_ub,lb,ub,options)
% This function solves linear programs encountered while finding a maximum
% output admissible set for an autonomous system. The function solves a
% series of problems which have each row of f as the objective. These
% problems are solved as minimization and maximization problems, to find if
% the value of f is within the bounds (lb,ub). The variable converged 
% returns 1 for the LPs which have objective functions
% within the bounds (lb,ub). Otherwise, it returns 0. MOAS is found when
% all(converged) is true.
%
% Inputs:
% f: objective functions
% cons: constraint matrix
% cons_lb: lower bounds
% cons_ub: upper bounds
% lb: lower bounds on f
% ub: upper bounds on f
%
% Outputs:
% converged: vector containing status of each sub-problem 
% exitflag: vector containing exitflag from each sub-problem
%% User defined variables

% tolerance
tol = 1e6*eps;

% bounds on variables
LB = -1e6*ones(length(f),1); 
UB = -LB;

% convert constraints into Ax<=b form
clhs = [cons;-cons];
crhs = [cons_ub;-cons_lb];

ef = zeros(size(lb));
converged = zeros(size(lb));
for i = 1:length(lb)
    
    % LP for lb (minimisation problem)
    [~,fval1,exitflag1,~] = linprog(f(i,:)',clhs,crhs,[],[],LB,UB,[],options);

    % LP for ub (minimisation problem)    
    [~,fval2,exitflag2,~] = linprog(-f(i,:)',clhs,crhs,[],[],LB,UB,[],options);

    if exitflag1==-2 || exitflag1== -5 || exitflag2== -2 || exitflag2== -5
        % infeasibile
        ef(i) = -2;        
    elseif exitflag1==-3 || exitflag2==-3
        % unbounded
        ef(i) = -3;
    elseif exitflag1<0 || exitflag2<0
        % ill conditioned problem
        ef(i) = -4;
    elseif (-fval1 +lb(i))<=-tol && (-fval2-ub(i))<=-tol
        converged(i) = 1;
        ef(i) = 1;
    end
    
    if all(ef>0)
       exitflag = 1;
    elseif any(ef == -2)  || any( ef == -4)
       exitflag = -2;
    elseif any(ef == -3)
       exitflag = -3;
    else
       exitflag = 0;
    end
    
end

end