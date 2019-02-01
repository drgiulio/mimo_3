function [ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order )
% function [ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order )
%
% The function calculates the achievable rates of the MAC and
% BC for given channels Hk and transmit covariance matrices Qk and Sk and
% a given decoding and encoding order.
%
% Input Specifications:
% H: Kx1 cell array with channel Hk per cell
% Q: Kx1 cell array with MxM matrix Qk per cell
% S: Kx1 cell array with NxN matrix Sk per cell
% order: 1xK vector in {1,...,K} with the BC encoding order
%
% Output Specifications:
% R_BC:  Kx1 array with achievable rates of the BC
% R_MAC: Kx1 array with achievable rates of the MAC

[M,N] = size(H{1});
K = length(H);

R_BC = zeros(K,1);
R_MAC = zeros(K,1);

%% Compute BC rates
for i = 1:K
    
    % Get user k
    k = order(i);
    
    % Sum of S_ks
    sum_S = 0;
    for j = i+1:K
        sum_S = sum_S + S{order(j)};
    end
    
    % Interference matrix
    C_in_k = eye(M) + H{k}*sum_S*H{k}';
    
    % BC Rate for user k
    R_BC(k) = real(log2(det(eye(N)+H{k}'*inv(C_in_k)*H{k}*S{k})));

end

%% Compute MAC rates
for i = K:-1:1
    
    % Get user k (reversed order!)
    k = order(i);
    
    % Interference matrix
    X_k = eye(N);
    for j = i-1:-1:1
        X_k = X_k + H{order(j)}'*Q{order(j)}*H{order(j)};
    end
    
    % BC Rate for user k
    R_MAC(k) = real(log2(det(eye(M)+H{k}*inv(X_k)*H{k}'*Q{k})));

end

