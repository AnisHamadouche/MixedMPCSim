%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Authoer: Yun Wu

%% 
function x = prox_l1(v, lambda)
% PROX_L1    The proximal operator of the l1 norm.
%
%   prox_l1(v,lambda) is the proximal operator of the l1 norm
%   with parameter lambda.

    x = max(0, v - lambda) - max(0, -v - lambda);
end