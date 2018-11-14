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

subplot(211); plot(t,xtraj); title('state') 
subplot(212); plot(t,utraj); title('input')