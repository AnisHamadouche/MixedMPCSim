%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Author: Anis Hamadouche

%% 
function x = prox_l0(v, lambda)
    thr = sqrt(2*lambda);
    x = wthresh(v,'h',thr);
end