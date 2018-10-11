% simulates the system with the MPC controller in closed-loop.
% input:    handle   handle to the MPC solver
%           x0       initial condition
%           A,B      continuous time system dynamics \dot{x}=Ax+Bu
%           Q,R      matrices that define the cost function (continuous
%                    time)
%           tau0d    basis functions evaluate at time index 0
%           Ts       sampling time
%           tend     simulation horizon
%
% output:   time     time vector
%           xtraj    state trajectory (closed-loop)
%           utraj    input trajectory (closed-loop)
%           iter     number of iterations needed for convergence of MPC
%                    solver
%           exectime execution time in seconds of MPC solver
%           J        closed-loop cost
%           opt      exitflag of the MPC solver (not implemented)
%
function [time,xtraj,utraj,iter,exectime,J] = simulatepdMPC(handle,moas)
%% User defined variables

% simulation time
tend = 15;

%%
% number of simulation steps
nsteps = tend/moas.sys.Ts+1;

% set up vectors to store the input and state trajectories
time = [0:moas.sys.Ts:tend]';
xtraj = zeros(nsteps,moas.sys.n);
utraj = zeros(nsteps,moas.sys.m);
iter = zeros(nsteps-1,1);         % store the number of iterations
exectime = zeros(nsteps-1,1);     % store the exection times


xtraj(1,:) = moas.sys.x0;
J = 0;

% just measure the execution time for warm-starts -> solve the optimization
% a first time
[handle,utmp,itertmp] = pdMPC_discrete('s',handle,xtraj(1,:)');
utraj(1,:) = utmp;

        
for k = 2:nsteps
    if(isnan(utmp))
        a = 0;
        error('Infeasible problem')
    end
    
    % solve optimization problem
    tic;
    [handle,utmp,iter(k-1)] = pdMPC_discrete('s',handle,xtraj(k-1,:)');
    exectime(k-1) = toc;
    utraj(k) = utmp;    
            
    % evaluate the solution, apply the input to the system and 
    % calculate the state at the next time instant
    xtraj(k,:) = xtraj(k-1,:)*moas.sys.A'+utraj(k)'*moas.sys.B';
    
    % calculate closed-loop cost
    J = J + 1/2*xtraj(k-1,:)*moas.sys.Q*xtraj(k-1,:)'*moas.sys.Ts ...
          + 1/2*utraj(k-1,:)*moas.sys.R*utraj(k-1,:)'*moas.sys.Ts;
    
end
