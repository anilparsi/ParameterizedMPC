clc;
clear;

%% Add CPLEX to matlab path
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');

%% Define system
sys = defineSystem();

%% Define approximation
apx = defineApproximation(sys);

%% Generate maximum output admissible set (MOAS)
moas = generateMOAS(sys,apx);

%% set-up the parameterized MPC solver
handle = generateMatrices(moas);

%% Closed Loop Simulation

% simulate in closed-loop
[t,xtraj,utraj,iter,exectime,J] = simulatepdMPC(handle,moas);

%% Plot trajectories
figure(20); clf;

subplot(211); plot(t,xtraj); title('state') 
subplot(212); plot(t,utraj); title('input')