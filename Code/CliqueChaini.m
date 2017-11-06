function[ rM, LogrM ] = CliqueChaini(dataFile, saveFile, chaini,...
    w, burnin, numblocksample, penc1, penc2)
%Run one chain of Gibbs sampler for the clique model
%
%Input ( dataFile, saveFile, chaini, w, burn, numblocksample ):
%   dataFile: data file that includes V, n, p
%   saveFile: file name for saved results
%   chaini: chain index 
%   w: window size for saved iterations
%   burnin: number of burn-in iterations
%   numblocksample: number of subsamplers used to approximate average
%   determinant of submatrices in LogDi.m
%Output [ rM, LogrM ]
%   rM: fiducial likelihood matrix for cliques
%   LogrM: a w by 1 matrix that records penalized fidicual likelihood of
%   each clique model without normalizing constant on the log scale

    load(dataFile); 
%     load(saveFile);
%     rng(w + p * 3 + chaini * 5)
    rng(w + p * 10 + chaini * 5)
    
    inds = randi(p, p, 1);
    inds_unique = unique(inds);
    k0 = length(inds_unique);
    GInd0 = ones(p, 1);
    for ind = 1:k0
        Gi = inds == inds_unique(ind);
        GInd0(Gi) = ind;
    end
    Gk0 = getGk(GInd0);
    Lris0 = ones(k0,1);
    for ind = 1:k0
        Lris0(ind) = Logri(V, Gk0{ind}, numblocksample);
    end

    t = 0;
    T = burnin + w;
    Lris_t = Lris0; Gk_t = Gk0; GInd_t = GInd0; 
    LogrM = ones(w, 1);
    GInds = ones(w, p);
    while t < T
        [Lris_t, Gk_t, GInd_t] = Clique1step(V, Lris_t, Gk_t, GInd_t,...
            numblocksample, penc1, penc2);
        t = t + 1;
        wi = t - burnin;
        if wi > 0
            LogrM(wi) = sum(Lris_t) - PenGk2(Gk_t, n, penc1, penc2);
            GInds(wi, :) = GInd_t;
        end
    end

    save(saveFile, 'GInd0', 'Gk0', 'Lris0', 'LogrM', 'burnin', 'w',...
        'Lris_t', 'Gk_t', 'GInds', 'penc1', 'penc2')

    rM = zeros(p,p);
    for wi = (w/2+1) : w
        GInd_wi = GInds(wi,:);
        rM = rM + getM(GInd_wi);
    end
    rM = rM/(w/2);

    save(saveFile, 'GInd0', 'Gk0', 'Lris0', 'rM', 'LogrM', 'burnin',...
        'w', 'Lris_t', 'Gk_t', 'GInds', 'penc1', 'penc2')
end
    