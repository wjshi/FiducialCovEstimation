function [ penk ] = PenGk(G, n)
% Compute log-transformed MDL penalty for the clique model
%   G: list of indices for each clique
%   n: number of observations
%     penk = 0;
    k = length(G);
    gis = cellfun(@(x) size(x, 2), G);
    p = sum(gis);
%     penk = sum(gis.^2) / 2 * log(n) + (p + 1) * log(k);
    penk = sum(gis.^2) / 4 * log(n) + (p + 1) * log(k);
end