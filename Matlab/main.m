clc;
clear;

%% Define system
sys = defineSystem();

%% Define approximation
apx = defineApproximation(sys);

%% Generate maximum output admissible set (MOAS)
moas = generateMOAS(sys,apx);

%% set-up the parameterized MPC solver
handle = generateSolver(moas);

%% Closed Loop Simulation

% simulate in closed-loop
[t,xtraj,utraj,iter,exectime,J] = simulatepdMPC(handle,moas);

%% Plot trajectories
figure(20); clf;

subplot(211); plot(t,xtraj(:,1)); title('state') 
line([t(1),t(end)],-[0.1,0.1],'color','k','linestyle','--')
axis tight
subplot(212); plot(t,utraj); title('input')
line([t(1),t(end)],[0.5,0.5],'color','k','linestyle','--')
line([t(1),t(end)],-[0.5,0.5],'color','k','linestyle','--')
axis tight
xlabel('time [s]')
legend('trajectory','bound')