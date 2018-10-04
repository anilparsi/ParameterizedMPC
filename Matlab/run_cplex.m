%% Change representation of constraints and run cplex
function [converged] = run_cplex(f1,cons,cons_lb,cons_ub,lb,ub,options)

tol = 1e6*eps;

% bounds because the problem is unbounded sometimes
LB = -1e6*ones(length(f1),1); 
UB = -LB;

% convert constraints into Ax<=b form
clhs = [cons;-cons];
crhs = [cons_ub;-cons_lb];

    
for i = 1:length(lb)
    % LP for lb (minimisation problem)
    [~,fval1,exitflag1,~] = cplexlp(f1(i,:)',clhs,crhs,[],[],LB,UB,[],options);

    % LP for ub (minimisation problem)    
    [~,fval2,exitflag2,~] = cplexlp(-f1(i,:)',clhs,crhs,[],[],LB,UB,[],options);

    if exitflag1<0 || exitflag2<0
        error('Infeasible problem. Please change approximation.')
    elseif (-fval1 +lb(i))<=-tol && (-fval2-ub(i))<=-tol
        converged(i) = 1;
    else
        converged(i) = 0;
    end
end

end