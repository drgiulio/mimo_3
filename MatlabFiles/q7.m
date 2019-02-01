%clear
close all
clc
%% Q7

% Load channels and MAC transmit covariance matrices
load('exampleMIMOBCs.mat')

% Get number of users K
Ptx_dB = 10;
Ptx = 10^(Ptx_dB/10);

% Compute optimal dual MAC transmit covariance matrices Q_k & sum capacity
[Q,Csum] = DualMACSumRateMaximization(H,Ptx);
