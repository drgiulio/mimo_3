function [ Q, Csum ] = DualMACSumRateMaximization( H, Ptx )
% function [ Q, Csum ] = DualMACSumRateMaximization( H, Ptx )
%
% The function calculates the sum capacity and the corresponding
% transmit covariance matrices of the K user MIMO MAC with a joint sum
% transmit power constraint.
%
% Inputs
% H: K x 1 cell array with channel Hk per cell
% Ptx: joint transmit power
%
% Outputs
% Q: K x 1 cell array of users' optimal transmit covariance Qk per cell
% Csum: sum capacity of the MIMO MAC

% parameters
[M,N] = size(H{1});
K = length(H);

% Set YALMIP options
options = sdpsettings('solver','sdpt3','verbose',0);

% Initialize optimization variables
Q = cell(K,1);
for k = 1:K
    Q{k} = sdpvar(M,M,'hermitian','complex');
end

% Define constraint set
powers = 0;
for k = 1:K
    powers = powers + trace(Q{k});
end
Constraints = [powers<=Ptx];

% Define objective function
Z = eye(N);
for k = 1:K
    Z = Z + H{k}'*Q{k}*H{k};
end
Objective = logdet(Z);

% Solve
optimize(Constraints, Objective, options);

% Retrieve solution
for k = 1:K
    Q{k} = value(Q{k});
end

% Compute sum capacity
det_arg = eye(N);
for k = 1:K
    det_arg = det_arg + H{k}'*Q{k}*H{k};
end
Csum = log2(det(det_arg));

end

