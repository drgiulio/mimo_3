function [S] = ptpTransform(Q,Heff)
% Function
% [C] = ptpTransform(Q,H)
%
% The effective transmit covariance matrix of
% the flipped channel H' is calculated to 
% preserve the achievable rate of the system.
%
% Inputs
% Q: M x M covariance matrix Q
% Heff: effective channel matrix H
% Output
% C: N x N transmit covariance matrix

% Perform SVD of channel matrix
[U,~,V] = svd(Heff,'econ');

% Compute transmit covariance matrix S according to p-t-p duality principle
S = V*U'*Q*U*V';

