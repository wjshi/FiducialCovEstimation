function[ CliqueM ] = getM( GInd )
%Convert the clique indices to a matrix form with entries indicating if the
%column index and row index belong to the same clique.
%   GInd: p by 1 vector of clique indices

p = length(GInd);
k = max(GInd);
CliqueM = diag(ones(1, p));
for ind1 = 1:k
    G_ind1 = find(GInd == ind1);
    g_ind1 = length(G_ind1);
    if g_ind1 > 1
        for ind2 = 1:g_ind1
            tuple = G_ind1(ind2);
            g_other = setdiff(G_ind1, tuple);
            CliqueM(tuple, g_other) = 1;
            CliqueM(g_other, tuple) = 1;
        end
    end
end
