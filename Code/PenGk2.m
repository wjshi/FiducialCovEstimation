function [ penk ] = PenGk2(G, n, penc1, penc2)
% Compute log-transformed MDL penalty for the clique model
%   penk = c1 * sum(gis.^2) * log(n) + c2 * (p + 1) * log(k)
%   G: list of indices for each clique
%   n: number of observations
%   penc1 & penc2: penalty parameters
%       Nopen: penc1 = penc2 = 0; 
%       oldpen: penc1 = 1/2, penc2 = 1; 
%       newpen: penc1 = 1/4, penc2 = 1.


    k = length(G);
    gis = cellfun(@(x) size(x, 2), G);
    p = sum(gis);
    penk = penc1 * sum(gis.^2) * log(n) + penc2 * (p + 1) * log(k);
end