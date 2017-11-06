function [newlogris, newGk, newGInd]...
= Clique1step(V, presentlogris, presentGk, presentGInd, numblocksample,...
penc1, penc2)
%One step Gibbs sampling for the clique model 
%
% Randomly choose a tuple from 1:p, and update its clique index from s to
% another group q. 
%   If the chosen tuple is the only member in clique s, propose to move it
%   to an existing clique; otherwise, propose to move it to a new clique 
%   or an existing clique.
%
%
%Input: (V, Gk, logrs, GInd, numblocksample)
%   V: stacked observations. V=[Y1,...,Yn]'.
%   presentlogris: a vector of log-fiducial likelihood without constant
%   for each tuple
%   presentGk: list of tuples for each clique
%   presentGInd: corresponding clique indicies for tuples 1:p
%   numblocksample: number of block determinant calculated. It is only used
%       when nchoosek(n, gi) is large
%
%Output: [newGk, newlogrs, newGInd]
%   newlogris: updated vector of log-fiducial likelihood without constant 
%   newGk: updated list of tuples for each clique
%   newGInd: updated corresponding clique indicies for tuples 1:p


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    newGk = presentGk; 
    newlogris = presentlogris;
    newGInd = presentGInd;
    
    [n, p] = size(V);
    k = length(presentGk);
    g = cellfun(@(x) size(x,2), presentGk);
    tuple = randsample(p, 1); %randomly choose a tuple from 1:p
    s = presentGInd(tuple); %clique index for the chosen tuple
    logrs = presentlogris(s);
    gs = g(s);
    
    if gs == 1 %the chosen tuple is the only element in clique s 
        rM_ratio = ones(k, 1);
        newlogrqs = ones(k, 1);
        for ind = setdiff(1:k, s)            
            logrq = presentlogris(ind);
            gq = g(ind); 
            newlogrqs(ind) = Logri( V, [presentGk{ind} tuple],...
            numblocksample );
%             rM_ratio(ind) = exp(newlogrqs(ind) - logrs - logrq) /...
%                 (n^gq * ((k - 1)/k)^(p + 1)); 
            rM_ratio(ind) = exp(newlogrqs(ind) - logrs - logrq -...
                DiffPen2(gs, gq, ind, k, p, n, penc1, penc2) );            
        end
        rM_cum = cumsum(rM_ratio/sum(rM_ratio));
        ra = rand(1);
        q = sum(ra > rM_cum) + 1;
        
        if q ~= s
            newGInd(tuple) = q; %assign clique index q to the tuple
            newGk{q} = [presentGk{q} tuple];
            newGk(s,:) = []; %remove clique s index list
            Ind_adj = find(newGInd > s); 
            if ~isempty(Ind_adj) %adjust clique indices
                newGInd(Ind_adj) = newGInd(Ind_adj) - 1;
            end
            newlogris(q) = newlogrqs(q); %update logrq
            newlogris(s) = []; %remove logrs
        end  
        
    else %clique s has multiple members
        rM_ratio = ones(k + 1, 1);
        newlogrqs = ones(k + 1, 1);
        newlogrs = Logri(V, setdiff(presentGk{s}, tuple), numblocksample); 
        %update logrs
        for ind = setdiff(1:k, s)            
            logrq = presentlogris(ind);
            gq = g(ind);
            newlogrqs(ind) = Logri(V, [presentGk{ind} tuple],...
                numblocksample);
%             rM_ratio(ind) = exp(newlogrqs(ind) + newlogrs - logrs -...
%                 logrq)/(n^(gq - gs + 1)); 
            rM_ratio(ind) = exp(newlogrqs(ind) + newlogrs - logrs -...
                logrq - DiffPen2(gs, gq, ind, k, p, n, penc1, penc2) );             
        end
        newlogrqs(k + 1) = Logri(V, tuple, numblocksample);
%         rM_ratio(k + 1) = exp(newlogrqs(k + 1) + newlogrs - logrs)/...
%             (n^(1 - gs) * (1 + 1/k)^(p + 1)); %add a clique
        rM_ratio(k + 1) = exp(newlogrqs(k + 1) + newlogrs - logrs -...
            DiffPen2(gs, 0, k+1 , k, p, n, penc1, penc2) ); %add a clique        
        rM_cum = cumsum(rM_ratio/sum(rM_ratio));
        ra = rand(1);
        q = sum(ra > rM_cum) + 1;
        
        if q ~= s
            newGInd(tuple) = q; %assign clique index q to the tuple
            newlogris(q) = newlogrqs(q); %update logrq
            newlogris(s) = newlogrs;
            newGk{s} = setdiff(presentGk{s}, tuple);
            if q == k+1
                newGk{q} = tuple;
            else 
                newGk{q} = [presentGk{q} tuple];
            end
        end         
        
    end
end