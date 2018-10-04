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
function [time,xtraj,utraj,iter,exectime,J,opt]=simulatepdMPC(handle,x0,A,B,Q,R,tau0d,Ts,tend)

% retrieve problem dimensions
n=size(A,1);
m=size(B,2);
N=length(tau0d)-1;

% number of simulation steps
nsteps=tend/Ts+1;

% set up vectors to store the input and state trajectories
time=[0:Ts:tend]';
xtraj=zeros(nsteps,n);
utraj=zeros(nsteps,m);
iter=zeros(nsteps-1,1);     % store the iteration numbers
exectime=zeros(nsteps-1,1); % store the exection times
opt=zeros(nsteps-1,1);      % is the problem solved to full optimality?

sysd=c2d(ss(A,B,eye(n),zeros(n,m)),Ts);
A=sysd.a;
B=sysd.b;

xtraj(1,:)=x0;
J=0;

% just measure the execution time for warm-starts -> solve the optimization
% a first time
[handle,utmp,itertmp]=pdMPC_discrete('s',handle,xtraj(1,:)');
utraj(1,:) = utmp;

        
for k=2:nsteps
    if(isnan(utmp))
        a = 0;
        keyboard;
    end
    % print out initial condition and interation number
%     fprintf('k: %i, x0: %f ',k,xtraj(k-1,1));
%     for kk=2:length(x0)
%         fprintf('%f ',xtraj(k-1,kk));
%     end
%     fprintf('\n');
    
    % solve optimization problem
    tic;
    [handle,utmp,iter(k-1)]=pdMPC_discrete('s',handle,xtraj(k-1,:)');
    exectime(k-1)=toc;
    utraj(k) = utmp;
    

            
    % evaluate the solution, apply the input to the system and 
    % calculate the state at the next time instant
    opt(k-1)=1;
    xtraj(k,:)=xtraj(k-1,:)*A'+utraj(k)'*B';
    
    % calculate closed-loop cost
    J=J+1/2*xtraj(k-1,:)*Q*xtraj(k-1,:)'*Ts+1/2*utraj(k-1,:)*R*utraj(k-1,:)'*Ts;
    
end
