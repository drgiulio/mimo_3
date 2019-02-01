clear
close all
clc
%% Q5

% Load channels and MAC transmit covariance matrices
load('exampleMIMOBCs.mat')

% Get number of users K
K = length(H);

% Define BC encoding orders
orders = [1:K; K:-1:1];
n_orders = size(orders,1);

% Compute BC rates by transformation of MAC channels
S = cell(n_orders,1);
R_BC = cell(n_orders,1);
R_MAC = cell(n_orders,1);
BC_powers = cell(n_orders,1);
MAC_powers = cell(n_orders,1);
for i = 1:n_orders
    
    % Get order
    order = orders(i,:);
    
    % Transformation to dual BC transmit matrices
    S{i} = MACtoBCtransform(Q,H,order);
    
    % Compute rates
    [R_BC{i}, R_MAC{i}] = MAC_BC_rates(H,Q,S{i},order);
    
    % Compute transmit powers for comparison
    BC_powers{i} = zeros(K,1);
    MAC_powers{i} = zeros(K,1);
    for k = 1:K
        BC_powers{i}(k) = real(trace(S{i,1}{k,1}));
        MAC_powers{i}(k) = real(trace(Q{i}));
    end
    
end