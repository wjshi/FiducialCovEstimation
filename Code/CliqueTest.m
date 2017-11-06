%Create a 100 by 100 covariance matrix Sigma with 10 equivariant cliques.
%Diagonal of Sigma is 1, for an off-diagonal position, it is 0.5 if the
%column and row indices belong to the same clique, 0 other wise. 


%Block diagonal format
Cor = 0.5; Var = 1; p = 10; n = 50; k = 3;
pndir = sprintf('SimData/CliqueM/p%dn%d', p, n);
mkdir(pndir);
datafile = sprintf('%s/Data.mat', pndir);
rng(k*11 + p*3)

G_end = [sort(randperm(p - 1, k - 1)) p];
Gk = cell(k, 1);
Gk{1} = 1:G_end(1);
for ind = 2:k
    Gk{ind} = (G_end(ind - 1) + 1) : G_end(ind);
end
g = cellfun(@(x) size(x,2), Gk);
Bs = cell(k,1);
for ind = 1:k
    gi = g(ind);
    Bs{ind} = repmat(Cor, gi) + diag(repmat(Var - Cor, [gi 1]));
end
Sigma = blkdiag(Bs{:});
G_ind = ones(p, 1); 
for ind = 1:k
    tuples = Gk{ind};
    G_ind(tuples) = ind;
end

rng(p * 11 + n * 7)
V = mvnrnd(zeros(1, p), Sigma, n);
Sn = V' * V / n;

save(datafile, 'Cor', 'Var', 'Gk', 'g', 'k', 'n', 'p',...
    'Sigma', 'G_ind', 'V', 'Sn');

figure
imagesc(Sigma)
colormap(flipud(gray))
colorbar
title('Sigma')

figure
imagesc(Sn)
colormap(flipud(gray))
colorbar
title('Sn')
%%
p = 100; n = 200;
pndir = sprintf('SimData/CliqueM/p%dn%d_newpen0', p, n);
datafile = sprintf('%s/Data.mat', pndir);
load(datafile);
numblocksample = 100;
w = 10000; burn = 10000; chain_num = 10;
c1 = 1/4; c2 = 1;

% parfor c = 1:chain_num
%     savefile = sprintf('%s/Chain%dw%d.mat', pndir, c, w);
%     CliqueChaini(datafile, savefile, c, w, burn, numblocksample, ...
%         c1, c2);
% end


subplot(2,2,1)
M = zeros(p);
for c = 1:chain_num
    savefile = sprintf('%s/Chain%dw%d.mat', pndir, c, w);
    load(savefile)
    M = rM + M;
    plot(LogrM)
    hold on
%     display(GInds(w, :))
end
title(sprintf('logrM trace plot (burnin = %d)', burn))
M = M/chain_num;


subplot(2,2,2)
imagesc(Sigma)
colormap(flipud(gray))
colorbar
title('Sigma')

subplot(2,2,3)
imagesc(Sn)
colormap(flipud(gray))
colorbar
title('Sn')

subplot(2,2,4)
imagesc(M)
colormap(flipud(gray))
colorbar
title(sprintf('%d chains', chain_num))

hold off
  


    