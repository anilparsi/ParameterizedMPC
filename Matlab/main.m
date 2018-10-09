clc;
clear;

%% Add CPLEX to matlab path
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');

%% Define system
sys = defineSystem();

%% Define approximation
apx = defineApproximation();

%% Generate maximum output admissible set (MOAS)
moas = generateMOAS(sys,apx);

%% set-up the MPC solver

% run make file 
run('..\pdMPC_discrete\src\MatlabInterface\make.m')

% path to the .txt files: inputs for C++
path  = '..\pdMPC_discrete\build\Data\MPCmat';

% Generate matrices and save into .txt files
generateMatrices(moas,path);

% initialize solver
handle = pdMPC_discrete('i',path);
%% Closed Loop Simulation

% simulate in closed-loop
[t,xtraj,utraj,iter,exectime,J,opt] = simulatepdMPC(handle,moas);
