function [S] = MACtoBCtransform(Q,H,order)
% Function [S] = MACtoBCtransform(Q,H,order)
% 
% Rate based duality transformation from the
% multiple access channel to the broadcast channel.
%
% Inputs
% Q: K x 1 cell array with matrix Qk per cell
% H: M x N cell array with channel Hk per cell
% order: 1 x K vector with BC encoding order
% Output
% S: K x 1 cell array of covariance matrix Sk per cell

% parameters
[M,N] = size(H{1});
K = length(H);

% initialization
S = cell(K,1);

% Compute BC covariance matrices from Q_ks (starting with last encoded user
% due to dependence of F_k on S_i, i > k!
for i = K:-1:1
    
    % Get user k
    k = order(i);
    
    % Compute interference matrix X_k
    X_k = eye(N);
    for j = i-1:-1:1
        X_k = X_k + H{order(j)}'*Q{order(j)}*H{order(j)};
    end
    
    % Compute hermitian matrix F_k
    sum_S = 0;
    for j = i+1:K
        sum_S = sum_S + S{order(j)};
    end
    F_k = eye(M) + H{k}*sum_S*H{k}';
    
    % Compute k-th effective MAC transmit covariance matrix
    Qk_eff = sqrtm(F_k)*Q{k}*sqrtm(F_k)';
    
    % Compute k-th effective channel
    H_keff = inv(sqrtm(F_k))'*H{k}*inv(sqrtm(X_k))';
    
    % Compute effective BC transmit covariance matrix
    S_keff = ptpTransform(Qk_eff,H_keff);
    
    % Compute S_k
    S{k} = inv(sqrtm(X_k))'*S_keff*inv(sqrtm(X_k));
    
end

