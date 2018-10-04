clc;
clear;
if(ispc)
    addpath('C:\Program Files\qpOASES-3.2.1\interfaces\matlab');
    addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\matlab\x64_win64');
end
%% Define System
DefineSystem;

%% Define approximation
DefineApproximation;

%% Generate maximum output admissible set (MOAS)
% give CPLEX path/solver path
% give t* directly:

%% set-up the MPC solver

% run make file 
run('..\pdMPC_discrete\src\MatlabInterface\make.m')

% Generate matrices
discreteTimeActiveSetGenMat(A,B,Q,R,Md,tau0d,x0,datafile);

% initialize solver

% set up the MPC problem
% handle = 0;
path2  = '..\..\..\pdMPC_discrete\build\Data\MPCmat';
handle=pdMPC_discrete('i',path2);
%% Closed Loop Simulation

% simulate in closed-loop
[t,xtraj,utraj,iter,exectime,J,opt]=simulatepdMPC(handle,x0,Ac,Bc,Q,R,tau0d,Ts,15);
