function[ Gk ] = getGk( GInd )
%Convert the clique indices to list of tuples contained in each clique
%
%   GInd: p by 1 vector of clique indices

k = max(GInd);
Gk = cell(k, 1);
for ind = 1:k
    Gk{ind} = find(GInd == ind)';
end
