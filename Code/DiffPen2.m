function [ diffpen ] = DiffPen2(oldclique_size, newclique_size, ...
    newclique, numcliques, dimp, samplesize, penc1, penc2)
% Compute log-transformed MDL penalty difference for the clique model
%   One index move from clique s to clique q
%   diffpen = PenGk2(newG, n, c1, c2) - PenGk2(oldG, n, c1, c2)
%
%   oldclique_size: size of clique s
%   newclique_size: size of clique q
%   newclique: new clique assignment 
%   numcliques: current number of cliques
%   dimp: row number of the covariate matrix
%   samplesize: number of observations
%   penc1 & penc2: penalty parameters
%       e.g. Nopen: penc1 = penc2 = 0; 
%            oldpen: penc1 = 1/2, penc2 = 1; 
%            newpen: penc1 = 1/4, penc2 = 1.

    
    if oldclique_size == 1
        diffpen = penc1 * 2 * newclique_size * log(samplesize) + ...
            penc2 * (dimp + 1) * log((numcliques - 1)/numcliques);
    elseif newclique == numcliques+1
        diffpen = penc1 * 2 * (1 - oldclique_size) * log(samplesize) + ...
            penc2 * (dimp + 1) * log((numcliques + 1)/numcliques);
    else 
        diffpen = penc1 * 2 * (1 + newclique_size - oldclique_size) * ...
            log(samplesize);
    end
end